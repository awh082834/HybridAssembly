#!/usr/bin/python

import sys
import os.path
#Assembly Filename
assemFile = sys.argv[1]
#Name of output file
outfileName = "IsolatedContigs.fasta"
headerFound = False
secondHeaderFound = False
try:
    with open(assemFile, 'r') as infile:
        lines = infile.readlines()
        #If the user wants all contigs other than the first, presumably
        #chromosomal contig, this flag will grab all of them and 
        #add them to a new file.
        for line in lines:
            if secondHeaderFound:
                outfile.write(line)
            #Looks for the second contig as that is the second
            #smallest, excluding the largest, the chromosomal contig
            if ">2" in line:
                secondHeaderFound = True
                outfile = open(outfileName, 'a')
                outfile.write(line)

except:
    print("Assembly file does not exist.")