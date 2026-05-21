#!/usr/bin/env python3
# Written by August Woerner
# this takes in a VCF file 
# and applies some basic thinning (approximating independence)
# thinning is applied using separation in the genetic map (ie, in units of cM)
# TODO: double-check code. and double check the AF tag
import sys
import os
import argparse
import gzip

parser = argparse.ArgumentParser(description="Let's thin some sites")
parser.add_argument('-m', '--map-vcf', dest='M', help="A VCF file with cM positions", default="")
parser.add_argument('-t', '--tag-in-vcf', dest='T', help="The FILTER field in the vcf file with the cM values", default="GPOS")
parser.add_argument('-c', '--min-cm-apart', dest='C', help="The minimum distanct in cM between markers", default=1.0, type=float)
parser.add_argument('-f', '--min-maf', dest='F', help="The minimum minor allele frequency", default=0.2, type=float)
parser.add_argument('-F', '--min-maf-tag', dest='AF', help="The minimum minor allele frequency TAG in the vcf file (same as cM file)", default="AF")
parser.add_argument('-v', '--transversions_only', dest='V', default=False, help="Only outputs transversions, ie, sites more resilient to DNA damage", action='store_true')


def getGmap(filename, tag,aftag, minAF):
    """
    Parses a VCF file decorated with cM positions
    produces a dictionary that ties the chr,pos -> cM values
    NOTE
    only markers with a minimum allele frequency (>=minMAF) are included
    filename: name of the file
    tag: the genetic map FILTER tag
    aftag: the genetic map allele frequency tag
    minAF: the smallest allowable minor allele frequency
    """
    if filename.endswith(".gz"):
        fh = gzip.open(filename, "rt")
    else:
        fh = open(filename, 'r')
    
    gmap = {}
    tag += "="
    taglen=len(tag)


    aftag += "="
    aftaglen=len(aftag)
    

    for line in fh:
        if line.startswith("#"):
            continue
        
        # column indexes come from the VCF standard
        sp = line.rstrip().split("\t")
        if len(sp) < 8:
            print("Bad file format for ", filename, file=sys.stderr)
            exit(1)
        fields = sp[7].split(";")
        cm = -1
        gotAf=False
        for f in fields:
            if f.startswith(tag):
                cm = float(f[taglen:])
                
                if sp[0] not in gmap:
                    gmap[sp[0]] = {}
                # dictionary of dictionaries
                # chr->pos = cm pos
                
                
            elif f.startswith(aftag):
                # get the alt allele frequency
                af = float( f[aftaglen:])
                gotAf=True
                # and ignore markers with *minor* allele frequency < than what is suggested
                if af < minAF:
                    cm = -1
                    break
                elif af > 1.0-minAF:
                    cm = -1
                    break
                    
        if cm > -1 and gotAf:
            gmap[sp[0]][sp[1]] = cm

            
    fh.close()
    return gmap

(results, xtra) = parser.parse_known_args(sys.argv[1:])

def thinSnpsSimple(gmap, cmsep, transversionsOnly):
    """
    takes in a gmap dictionary of dictionaries
    parses stdin, prints to standard out a VCF file that has 
    a separation between all snps of at least cmsep
    """

    prevcm = 0
    prevchrom="chr0"
    for line in sys.stdin:
        if line.startswith("##"):
            print(line, end="")
            continue
        elif line.startswith("#"):
            print('##INFO=<ID=GPOS,Number=A,Type=Float,Description="cM position in genetic map">')
            print(line, end="")
            continue
        sp = line.rstrip().split("\t")
        chrom = sp[0]
        pos = sp[1]
        # assumes vcf file in sane
        if chrom not in gmap:
            break
        
        # note, this assumes that the alt allele is annotated, regardless of the genotype called (e.g., it's not . even when the call is 0/0)
        # it's how this package genotypes, so it's OK
        # but this logic may not port to other callers
        # in practice this removes a lot of markers. this is a bit heavy handed, too.
        if transversionsOnly:   
            # slightly simpler to just remove transitions
            if sp[3] == 'C' and sp[4] == 'T':
                continue
            elif sp[3] == 'T' and sp[4] == 'C':
                continue           
            elif sp[3] == 'G' and sp[4] == 'A':
                continue 
            elif sp[3] == 'A' and sp[4] == 'G':
                continue 
            if len(sp[3]) !=1 or len(sp[4]) != 1: # additional check for snp-ness
                continue
            
        cms = gmap[chrom]
        if pos not in cms:
            continue
        cm = cms[pos]
        
        # ensure that cMs are increasing.
        #if prevcm>cm and prevchrom==chrom:
        #    print("genetic map error", chrom, pos, cm, prevcm, file=sys.stderr)
        #    exit(1)
        
        if sp[7] == "" or sp[7] == ".":
            sp[7] = "GPOS=" + str(cm)
        else:
            sp[7] += ";GPOS=" + str(cm)
            
        if prevchrom != chrom:           
            print("\t".join(sp))
            prevchrom=chrom
            prevcm=cm
        elif cm - prevcm >= cmsep:
            print("\t".join(sp))
            prevcm=cm
        


if len(xtra):
    print("Extra arguments detected...", xtra, file=sys.stderr)
    exit(1)

if results.M == "" or not os.path.exists(results.M):
    print("Failed to read the genetic map VCF file. Looks like you mis-typed it", file=sys.stderr)
    exit(1)


gmap = getGmap(results.M, results.T, results.AF, results.F)

thinSnpsSimple(gmap, results.C, results.V)

