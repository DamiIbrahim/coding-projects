#!/usr/bin/env python3
# Name: Dami Ibrahim (oaibrahi)
# Group Members: None
'''Program docstring: This program reads a file from STDIN and returns the unique and essential subsequence of each tRNA. It also minimizes all 22 sets of
 tRNA such that no string in a unique subsequence is a substring of another unique set'''
import sys
sys.stdout.reconfigure(encoding='utf-8') #Prints out special characters
class FastAreader :
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
    def doOpen (self):
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)
    def readFasta (self):
        header = ''
        sequence = ''
        with self.doOpen() as fileH:
            header = ''
            sequence = ''
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                if not line: # we are at EOF
                    return header, sequence
                line = fileH.readline()
            header = line[1:].rstrip()
            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
        yield header,sequence
########################################################################
# Main
# Here is the main program
#
########################################################################
class tRNAs:
    '''
    Defines the init, buildUnique, buildEssential, and printEssentials methods and contains the tRNA list.
    '''
    tRNAlist = []
    def __init__(self, head, seq):
        ''' Captures and defines self.seq, self.head, self.unique, self.essential, and self.powerSet that will be in the later methods including main.
        Also builds the powerset.'''
        self.head = head.replace(' ', '') #removes any spaces from the header line
        self.seq = seq.upper()#Makes the sequence uppercase
        self.unique =  set()
        self.essential = set()
        self.powerSet = set()
        for start in range(0,len(seq)): #builds powerset
            for end in range(start+1, len(seq) +1):
                self.powerSet.add(seq[start:end])
        tRNAs.tRNAlist.append(self)
    def buildUnique(self):
        ''' Finds all the unique subsequences'''
        superSet = set()
        for tRNA in tRNAs.tRNAlist:
            if tRNA is not self:
                superSet = superSet | tRNA.powerSet #Computes the union of the tRNA sets
        self.unique = self.powerSet - superSet #Removes union from current tRNA set
        #print(self.unique)
    def _buildEssential(self):
        ''' Finds all the essential subsequences'''
        superSet = set()
        for seq in self.unique:
            if (seq[1:] in self.unique) or (seq[:-1] in self.unique): #cuts off the letters from each end of the sequence
                superSet.add(seq)
        self.essential = self.unique - superSet
        #print(self.essential)
    def printEssentials(self):
        ''' Prints out the header and sequence as well as the sorted essentials'''
        essentials = {}
        print(self.head)
        print(self.seq)
        for item in self.essential:
            start = 0
            while self.seq.find(item, start) != -1:
                essentials[self.seq.find(item, start)] = item
                start += 1
        for i in sorted(essentials.keys()): #sorts the essentials by key value
            print('.'* i + essentials[i])

def main(inCL=None):
    '''Sorts through the tRNAs in the list by header and calls upon all the methods'''
    myReader = FastAreader(inCL)
    for head, seq in myReader.readFasta() :
        tRNAs(head, seq)
    sortedList = sorted(tRNAs.tRNAlist, key = lambda item: item.head) #sorts by tRNA header
    for tRNA in sortedList:
        tRNA.buildUnique() #calls unique method
        tRNA._buildEssential() #calls essential method
        tRNA.printEssentials() #calls printEssential method
if __name__ == "__main__":
    main(inCL = 'bos-tRNA-7.fa')
