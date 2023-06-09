# coding-projects
A few coding projects I have worked on during my time as an undergraduate at the University of California, Santa Cruz. 

The simpleMap python program creates a basic seed and extend mapping algorithm for aligning short query strings
to a larger target string. Additionally, it clusters seeds to determine probable alignments and finds the optimal
local alignment between two strings. 
Run the program with: python simpleMap.py hg19.chr3.9mb.fa NA12878.ihs.chr3.100kb.1.fastq.tiny --log=DEBUG

The pairwise python program is a program that my lab partners and I (Shivani Rao & Arshia Kapil) 
developed in our research lab that returns a pairwise correlation matrix of all cell adhesion cadherin/protocadherin genes 
based on their pearson correlation coefficient. In addition, it also returns the pairs that reject the null hypothesis
which is when the pearson coefficient is greater than or equal to the critical r value based on the level of significance. 
A user can input the level of significance of their choice and also change the tsv file to a dataset of their choice.
In our research, we split the dataset found from the Human Protein Atlas into three different cell types-Neural, Enodthelial, and Epithetlial-and used the 0.05 level of significance. 

The DAG python notebook finds the longest path between two nodes in an edge-weighted directed acyclic graph given an adjaceny list
This program is based off problem 15 Rosalind "Find the Longest Path in a DAG". 

The mt.tRNA folder contains a python file and two fasta files. Defects in mt.tRNA can lead to diseases and in order to identify the mutations, it is important to identify the abduances of these molecules in a population. This can be done by finding the unique subsequences in the 22 mt.tRNA that can be used as tags to assess the abundance. The program does that by reading in files of fasta sequences from
STDIN, finding the unique subsequences in each single tRNA that is not present amongst any other tRNA sets. The program outputs the tRNA name, tRNA sequence, and the each unique element to STDOUT. 
