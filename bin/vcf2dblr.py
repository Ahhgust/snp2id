#!/usr/bin/env python3
# Written by August Woerner
# this creates a DBLR-compliant "weights" input file.
# presented as the normalized likelihoods from the PL tag (VCF file format).
# this ONLY support biallelic genotypes.
# AND
# it "lies" about what alleles are observed (IMO, for low pass WGS not the best criterion)
# it does not lie about what the corresponding genotype probabilities are, though.
# in my thinking, it is the probabilities that matter, not the observations.

import argparse
import os
import sys


def die(args):
    print(args, file=sys.stderr)
    exit(1)

# extract (and ignore) metadata
for line in sys.stdin:
    if len(line)==0 or line[0] != "#":
        die("Invalid file format from stdin")

    if len(line)> 1 and line[1] != '#':
        break

cols = line.rstrip().split("\t") # column headers

print("File Version: 3")

locusCounter=0
# file body...
for line in sys.stdin:
    sp = line.rstrip().split("\t")
    meta=sp[8].split(":")
    geno=sp[9].split(":")
    
    # skip no calls
    # it's either that or emit flat likelihoods...
    if geno[0][0] == '.':
        continue
    
    plindex=1
    if len(meta) < 3:
        continue # per the spec, index 1 must be the gt.
    elif meta[plindex] != 'PL':
        for plindex in range(3, len(meta)):
            if meta[plindex]=='PL':
                break
    if meta[plindex] != "PL":
        continue
    
    pls = [int(i) for i in geno[plindex].split(",")]
    # biallelic only...
    if len(pls) != 3:
        continue
      
    probs = [ 10**(i/-10.) for i in pls ] # Qscores to likelihood (ratios)
    s = sum(probs)  # 
    probs = [ i/s for i in probs ] # then normalized into "probabilities"
    
    # locus name: chr1:400_A_C
    lname = sp[0] + "_" + sp[1] + "_" + sp[3] + "_" + sp[4]
    locusCounter +=1
    print("Locus ", locusCounter, " (", lname , ")", sep="")
    print("Observed: ", sp[3] , ",", sp[4], sep="")
    print(sp[3] , "\t", sp[3], "\t\t", probs[0])
    print(sp[3] , "\t", sp[3], "\t\t", probs[1])    
    print(sp[4] , "\t", sp[4], "\t\t", probs[2])
    #print(meta, geno, line, sep="\n")
    #exit(1)







