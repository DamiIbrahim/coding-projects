#Citations: Vidhata Singh & Savi Dhoat

import array
import sys
import numpy as np
import pysam
import argparse
import logging 
import time
logger = logging.getLogger()


class MinimizerIndexer(object):
    """ Simple minimizer based substring-indexer. 
    
    """
    def __init__(self, targetString, w, k, t):
        """ The target string is a string/array of form "[ACGT]*".
        
        Stores the lexicographically smallest k-mer in each window of length w, such that w >= k positions. This
        smallest k-mer is termed a minmer. 
        
        If a minmer occurs in the target sequence more than t times as a minmer then it is omitted from the index, i.e. if the given minmer (kmer) is a minmer
        in more than t different locations in the target string. Note, a minmer may be the minmer for more than t distinct windows
        and not be pruned, we remove minmers only if they have more than t distinct occurrences as minmers in the sequence.
        """
        
        self.targetString = targetString
        self.w = w
        self.k = k
        self.t = t 
        self.minimizerMap = {}

        for i in range(len(self.targetString) - self.w + 1) : #iterates over all possible substrings within the self.w window substring
            length = self.targetString[i:i + self.w]
            minmerLst = []
            for j in range(len(length) - self.k + 1) :
                minmerLst.append((length[j:j + self.k], i + j)) #add substring and its starting index within the self.w window
            minmer = min(minmerLst) #find the smallest self.k window substring and its starting index
            if minmer[0] not in self.minimizerMap : #Add the minimum length substring of self.k as a key and the starting index as a value
                self.minimizerMap[minmer[0]] = (minmer[1],)
            elif minmer[1] not in self.minimizerMap[minmer[0]] :#if the starting index for the minimum substring is not present then add it as a value to the dictionary
                self.minimizerMap[minmer[0]] = self.minimizerMap[minmer[0]] + (minmer[1],)

    def getMatches(self, searchString):
        """ Iterates through search string finding minmers in searchString and
        yields their list of minmer occurrences in targetString, each as a pair of (x, (y,)*N),
        where x is the index in searchString and y is an occurrence in targetString.
        """
        for i in range(len(searchString) - self.k + 1):
            if searchString[i:i + self.k] in self.minimizerMap: #checks if the window substring is in the dictionary, if it is then yield the starting indices of the substring
                yield(i, self.minimizerMap[searchString[i:i + self.k]])

class SeedCluster:
    """ Represents a set of seeds between two strings.
    """
    def __init__(self, seeds):
        """ Seeds is a list of pairs [ (x_1, y_1), (x_2, y_2), ..., ], each is an instance of a seed 
        (see static cluster seeds method below: static methods: https://realpython.com/blog/python/instance-class-and-static-methods-demystified/)
        """
        seeds = list(seeds)
        seeds.sort()
        self.seeds = seeds
        # Gather the minimum and maximum x and y coordinates
        self.minX = seeds[0][0]
        self.maxX = seeds[-1][0]
        ys = [y for x,y in seeds] # python3: map(lambda...) changed to list comprehension
        self.minY = min(ys)
        self.maxY = max(ys)

    @staticmethod
    def clusterSeeds(seeds, l):
        """ Cluster seeds (k-mer instances) in two strings. This is a static constructor method that creates a set
        of SeedCluster instances.
        
        Here seeds is a list of tuples, each tuple has the form (x, (y_1, y_2, ... )), where x is the coordinate
        in the first string and y_1, y_2, ... are coordinates in the second string. Each pair of x and y_i
        is an occurence of a shared k-mer in both strings, termed a *seed*, such that the k-mer 
        occurrence starts at position x in the first string and starts at position y_i in the second string.
        The input seeds list contains no duplicates and is sorted in ascending order, 
        first by x coordinate (so each successive tuple will have a greater  
        x coordinate), and then each in tuple the y coordinates are sorted in ascending order.
        
        Two seeds (x_1, y_1), (x_2, y_2) are *close* if the absolute distances | x_2 - x_1 | and | y_2 - y_1 |
        are both less than or equal to l.   
        
        Consider a *seed graph* in which the nodes are the seeds, and there is an edge between two seeds if they
        are close. clusterSeeds returns the connected components of this graph
        (https://en.wikipedia.org/wiki/Connected_component_(graph_theory)).
        
        The return value is a Python set of SeedCluster object, each representing a connected component of seeds in the 
        seed graph.

        """ 

        cluster = set()
        seedLst = []
        for i, j in seeds: #add each individual seed to the list
            for item in j:
                seedLst.append((i, item))

        if len(seedLst) == 0:
            return cluster
        visited = [seedLst[0]] #Seeds that have been considered
        notVisited = [seedLst[0]] #Seeds that have not been considered
        seedCluster = [seedLst[0]]

        while len(visited) < len(seedLst):
            if len(notVisited) > 0: #Checks if there are seeds that have not been visited
                next = notVisited.pop(0)
            else:
                final = SeedCluster(seedCluster) #add cluster to the final cluster and reset the list
                cluster.add(final)
                seedCluster = []

                for seed in seedLst:
                    if seed not in visited:
                        visited.append(seed)
                        notVisited.append(seed) #ensures that the seed will be visited in the next iteration
                        seedCluster.append(seed)
                        break

            x1, y1 = next[0], next[1]
            for seed in seedLst: #iterate through the seeds that have not been considered yet
                if seed not in visited:
                    x2, y2 = seed[0], seed[1]
                    if abs(x1 - x2) <= l and abs(y1 - y2) <= l: #Check if the current seed is connected to the current cluster, if so, then add seed to all the lists
                        visited.append(seed)
                        notVisited.append(seed)
                        seedCluster.append(seed)
            if len(visited) == len(seedLst): #check if all seeds have been considered and add to final cluster
                final = SeedCluster(seedCluster)
                cluster.add(final)
        return cluster


class SmithWaterman(object) :
    def __init__(self, string1, string2, gapScore=-2, matchScore=3, mismatchScore=-3) :
        """ Finds an optimal local alignment of two strings.

        Implements the Smith-Waterman algorithm:
        https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm

        """

        self.matrix = np.zeros(shape=[len(string1) + 1, len(string2) + 1], dtype=int) #initialize a matrix of zeros to store the scores
        self.string1 = string1
        self.string2 = string2
        self.gapScore = gapScore
        self.matchScore = matchScore
        self.mismatchScore = mismatchScore

        for i in range(1, len(string1) + 1) : #iterates over the rows and columns of the matrix starting from the second row and column
            for j in range(1, len(string2) + 1) :
                indel1 = self.matrix[i][j - 1] + self.gapScore #compute the score for adding a gap
                indel2 = self.matrix[i - 1][j] + self.gapScore
                if string1[i - 1] == string2[j - 1] : #check if its a match or mismatch
                    alignment = self.matrix[i - 1][j - 1] + self.matchScore
                else :
                    alignment = self.matrix[i - 1][j - 1] + self.mismatchScore
                self.matrix[i][j] = max(indel1, indel2, alignment, 0) #set the score to the max value of the three possible scores

    def getAlignment(self) :
        """ Returns an optimal local alignment of two strings. Alignment
        is returned as an ordered list of aligned pairs.

        e.g. For the two strings GATTACA and CTACC an optimal local alignment
        is (GAT)TAC(A)
             (C)TAC(C)
        where the characters in brackets are unaligned. This alignment would be returned as
        [ (3, 1), (4, 2), (5, 3) ]

        In cases where there is a tie between optimal sub-alignments use the following rule:
        Let (i, j) be a point in the edit matrix, if there is a tie between possible sub-alignments
        (e.g. you could choose equally between different possibilities), choose the (i, j) to (i-1, j-1)
        (match) in preference, then the (i, j) to (i-1, j) (insert in string1) in preference and
        then (i, j) to (i, j-1) (insert in string2).
        """

        maxScore = np.where(self.matrix == self.getMaxAlignmentScore()) #returns where the max value in the matrix occurs


        x = int(maxScore[0])
        y = int(maxScore[1])

        pairs = []
        while x and y > 0: #iterates back through the matrix until it reaches the top left square to find the optimal path
            leftSquare = self.matrix[x][y - 1] + self.gapScore
            topSquare = self.matrix[x - 1][y] + self.gapScore
            if self.string1[x - 1] == self.string2[y - 1] :
                alignment = self.matrix[x - 1][y - 1] + self.matchScore
            else :
                alignment = self.matrix[x - 1][y - 1] + self.mismatchScore

            if max(leftSquare, topSquare, alignment) == alignment : #if the diagonal score is the highest move diagonal towards the top left and add the current position to the list
                x -= 1
                y -= 1
                pairs.append((x, y))
            elif max(leftSquare, topSquare, alignment) == leftSquare : #if the left score is the max value then move left.
                y -= 1
            else : #if the top score is the max value move up.
                x -= 1

        pairs = pairs[::-1] #Reverse the pairs since we're moving backwards

        return pairs

    def getMaxAlignmentScore(self) :
        """ Returns the maximum alignment score
        """
        return np.max(self.matrix) #returns max value in the matrix
    
def simpleMap(targetString, minimizerIndex, queryString, config):
    """ Function takes a target string with precomputed minimizer index and a query string
    and returns the best alignment it finds between target and query, using the given options specified in config.
    
    Maps the string in both its forward and reverse complement orientations.

    """
    bestAlignment = [None]
    
    def mapForwards(queryString):
        """ Maps the query string forwards
        """
        # Find seed matches, aka "aligned kmers"
        seeds = list(minimizerIndex.getMatches(queryString))
        
        # For each cluster of seeds
        for seedCluster in SeedCluster.clusterSeeds(list(seeds), l=config.l):
            
            # Get substring of query and target to align
            queryStringStart = max(0, seedCluster.minX - config.c) # Inclusive coordinate
            queryStringEnd = min(len(queryString), seedCluster.maxX + config.k + config.c) # Exclusive coordinate
            querySubstring = queryString[queryStringStart:queryStringEnd]
            
            targetStringStart = max(0, seedCluster.minY - config.c) # Inclusive coordinate
            targetStringEnd = min(len(targetString), seedCluster.maxY + config.k + config.c) # Exclusive coordinate
            targetSubstring = targetString[targetStringStart:targetStringEnd]
            
            #print "target_aligning", targetStringStart, targetStringEnd, targetSubstring
            #print "query_aligning", queryStringStart, queryStringEnd, querySubstring
            
            # Align the genome and read substring
            alignment = SmithWaterman(targetSubstring, querySubstring, 
                                      gapScore=config.gapScore, 
                                      matchScore=config.matchScore,
                                      mismatchScore=config.mismatchScore)
            
            # Update best alignment if needed
            if bestAlignment[0] == None or alignment.getMaxAlignmentScore() > bestAlignment[0].getMaxAlignmentScore():
                bestAlignment[0] = alignment
        
        return bestAlignment
    
    def reverseComplement(string):
        """Computes the reverse complement of a string
        """
        rMap = { "A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
        return "".join(rMap[i] for i in string[::-1])
                
    # Run mapping forwards and reverse
    mapForwards(queryString)
    mapForwards(reverseComplement(queryString))
    
    return bestAlignment[0]

class Config():
    """ Minimal configuration class for handing around parameters
    """
    def __init__(self):
        self.w = 30
        self.k = 20
        self.t = 10
        self.l = 30
        self.c = 100
        self.gapScore=-2
        self.matchScore=3
        self.mismatchScore=-3
        self.logLevel = "INFO"
        
def main():
    # Read parameters
    config = Config()
    
    #Parse the inputs args/options
    parser = argparse.ArgumentParser(usage="target_fasta query_fastq [options]") # , version="%prog 0.1")

    parser.add_argument("target_fasta", type=str,
                        help="The target genome fasta file.")
    parser.add_argument("query_fastq", type=str,
                        help="The query sequences.")
    
    parser.add_argument("--w", dest="w", help="Length of minimizer window. Default=%s" % config.w, default=config.w)
    parser.add_argument("--k", dest="k", help="Length of k-mer. Default=%s" % config.k, default=config.k)
    parser.add_argument("--t", dest="t", help="Discard minmers that occur more frequently " 
                                            "in the target than t. Default=%s" % config.w, default=config.w)
    parser.add_argument("--l", dest="l", help="Cluster two minmers into the same cluster if within l bases of"
                                            " each other in both target and query. Default=%s" % config.l, default=config.l)
    parser.add_argument("--c", dest="c", help="Add this many bases to the prefix and suffix of a seed cluster in the"
                                            " target and query sequence. Default=%s" % config.c, default=config.c)
    parser.add_argument("--gapScore", dest="gapScore", help="Smith-Waterman gap-score. Default=%s" % 
                      config.gapScore, default=config.gapScore)
    parser.add_argument("--matchScore", dest="matchScore", help="Smith-Waterman match-score. Default=%s" % 
                      config.gapScore, default=config.gapScore)
    parser.add_argument("--mismatchScore", dest="mismatchScore", help="Smith-Waterman mismatch-score. Default=%s" % 
                      config.mismatchScore, default=config.mismatchScore)
    parser.add_argument("--log", dest="logLevel", help="Logging level. Default=%s" % 
                      config.logLevel, default=config.logLevel)
    
    options = parser.parse_args()
    
    # Parse the log level
    numeric_level = getattr(logging, options.logLevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.logLevel)
    
    # Setup a logger
    logger.setLevel(numeric_level)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(numeric_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.debug("Established logger")
    
    startTime = time.time()
    
    # Parse the target sequence and read the first sequence
    with pysam.FastaFile(options.target_fasta) as targetFasta:
        targetString = targetFasta.fetch(targetFasta.references[0])
    logger.info("Parsed target string. Length: %s" % len(targetString))
    
    # Build minimizer index
    minimizerIndex = MinimizerIndexer(targetString.upper(), w=options.w, k=options.k, t=options.t)
    minmerInstances = sum(map(len, minimizerIndex.minimizerMap.values()))
    logger.info("Built minimizer index in %s seconds. #minmers: %s, #minmer instances: %s" %
                 ((time.time()-startTime), len(minimizerIndex.minimizerMap), minmerInstances))
    
    # Open the query files
    alignmentScores = [] # Array storing the alignment scores found
    with pysam.FastqFile(options.query_fastq) as queryFastq:
        # For each query string build alignment
        for query, queryIndex in zip(queryFastq, range(sys.maxsize)):  # python3: xrange to range, maxint to maxsize
            print (queryIndex) # python3
            alignment = simpleMap(targetString, minimizerIndex, query.sequence.upper(), config)
            alignmentScore = 0 if alignment is None else alignment.getMaxAlignmentScore()
            alignmentScores.append(alignmentScore)
            logger.debug("Mapped query sequence #%i, length: %s alignment_found?: %s "
                         "max_alignment_score: %s" % 
                         (queryIndex, len(query.sequence), alignment is not None, alignmentScore)) 
            # Comment this out to test on a subset
            #if queryIndex > 100:
            #    break
    
    # Print some stats
    logger.critical("Finished alignments in %s total seconds, average alignment score: %s" % 
                    (time.time()-startTime, float(sum(alignmentScores))/len(alignmentScores)))
    
if __name__ == '__main__':
    main()
