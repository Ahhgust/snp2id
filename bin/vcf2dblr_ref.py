#!/usr/bin/env python3
# Written by August Woerner
# this creates a DBLR-compliant ground truth (taken to be accurate single source genotypes)

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
if len(cols) < 10:
    print("Incomplete vcf detected...", cols, file=sys.stderr)
    exit(1)
    
sample=cols[9]

print("Sample Name", "Locus", "Allele 1", "Allele 2", sep="\t")


for line in sys.stdin:
    sp = line.rstrip().split("\t")
    geno=sp[9].split(":")
    
    # skip no calls
    if geno[0][0] == '.':
        continue

    locus=sp[0] + "_" + sp[1] + "_" + sp[3] + "_" + sp[4]
    
    alleles=[sp[3],  sp[4]]
    
    call = geno[0].split("/")
    if len(call)==1:
        call = geno[0].split("|")
        
        if len(call)==1:
            print("Diploid calls only...", line, sep="\n", file=sys.stderr)
            continue
    
    acall = [alleles[ int(c) ] for c in call]    
    
    print(sample, locus, acall[0], acall[1], sep="\t")







