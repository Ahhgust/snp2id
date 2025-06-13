#!/usr/bin/env python3
# this takes in a VCF file (stdin)
# and creates a dblr-compliant allele frequency file (population specific)

import sys
import os
import argparse

parser = argparse.ArgumentParser(description="Let's make some allele frequency tables!")
parser.add_argument('-p', '--population', dest='P', help="The population in question. AF_{pop} and AC_{pop} tags are needed in the VCF", default="NFE")
parser.add_argument('-s', '--sex-chromosomes', dest='S', help="Should the analysis only consider the autosomes (default), or also consider X/Y markers?", action='store_true')
 
 
(results, xtra) = parser.parse_known_args(sys.argv[1:])

if len(xtra):
    print("Extra arguments detected...", xtra, file=sys.stderr)
    exit(1)
    

# gnomad population label.
pop = results.P.lower()
# tag presents as AF_afe=0.990"
afTag = "AF_" + pop + "="
anTag = "AN_" + pop + "="


print("Locus,Allele,Allele Seq,Count,Frequency")

for line in sys.stdin:
    if line.startswith("#"):
        continue
    sp= line.rstrip().split("\t")
    if not results.S:
        # non-numeric chromosome...
        if sp[0][-1] < '0' or sp[0][-1] > '9':
            break

    
    refAllele= sp[3]
    altAlleles=sp[4].split(",")
    
    tags = sp[7].split(";")
    an=[-1]
    af = [-1];
    for t in tags:
        if t.startswith(afTag): # 1+ alt alleles
            af = t.split("=")[1]
            af = [float(f) for f in af.split(",")] # alt allele frequencies in the order given
        elif t.startswith(anTag):
            an = t.split("=")[1]
            an = [int(f) for f in an.split(",")]
    
    if sum(af) < 0 or sum(an) < 0:
        print("No allele frequency data for\n", line, file=sys.stderr)
        continue
        
    # as w/ the genotype weights file, chrom, position, ref, alt
    locus=sp[0] + "_" + sp[1] + "_" + sp[3] + "_" + sp[4]
    
    altFreq=0.
    nalleles=len(af)
    for i in range(nalleles):
        print(locus, ord(altAlleles[i])-64, altAlleles[i], round(an[i]*af[i]), af[i], sep=",")
        altFreq += af[i]
    
    anMean =round( sum(an)/len(an) ) # should be singular, but just in case.
    refFreq = 1.0-altFreq
    
    print(locus, ord(refAllele)-64, refAllele, round(refFreq*anMean), refFreq, sep=",")

    
    
    