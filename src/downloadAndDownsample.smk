import os, sys, random

VERSION=0.02
project_dir = os.path.join(workflow.current_basedir, "..")

# select 5 individuals from 4 populations
# at present, these are just from 30x recall of the 1000 Genomes Project
# (we have metadata for hgdp as well)
nsamps=5
kpops=["CEU", "MXL", "ASW", "JPT"]
# and downsample to the specified coverages
# modified to include 30x (high pass)
coverages=[0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.50, 1.0, 5., 10., 30.]

# and stuff needed to relate the above to
# actually downloading and downsampling
url_lut={} # 2d; pop -> sample -> url
url_flat={} # 1D pop/sampleid -> url
kgfile=os.path.join(project_dir, "pop_metadata", "1kg.urls")

URLCOL=0
POPCOL=10
SCOL=9

random.seed(10)

# parse the 1KG data
with open(kgfile) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        sp = line.rstrip().split("\t")
        #annoyingly, some identifiers have leading/trailing whitespace
        (u, p, s) = (sp[URLCOL].strip(), sp[POPCOL].strip(), sp[SCOL].strip())
        #print("'", u, "' '", p, "' '",s, "',", sep="")
        if p not in kpops:
            continue
        if p not in url_lut:
            url_lut[p]=dict()
        url_lut[p][s] = u

# randomly select nsamps individuals from each pop

for pop in kpops:
    if pop not in url_lut:
        print("Missing population: ", pop, list(url_lut.keys()), file=sys.stderr)
        exit(1)
        
    inner = url_lut[ pop ]
    if len(inner) < nsamps:
        print("Not enough samples available for", pop, file=sys.stderr)
        exit(1)
    

    toremove=list(inner.keys() )
    random.shuffle( toremove )
    toremove = toremove[nsamps:]
    
    for s in toremove:
        del inner[s]
    
    for s, u in inner.items():
        url_flat[ pop + "." + s  ] = u


rule whatcha_want:
    input: 
        expand("original_crams/{fname}.cram", fname= list(url_flat.keys())),
        expand("original_crams/{fname}.cram.cov", fname= list(url_flat.keys())),
        expand("logs/{fname}.ds.complete", fname= list(url_flat.keys()) )

# neat trick; treat the url as a parameter;
# that way, snakemake treats this as a generator function (no input, but makes output)
rule download_crams:
    output:
        "original_crams/{fname}.cram"  
    params:
        url = lambda wildcards: url_flat[ wildcards.fname ]        
    threads:
        1
    log:
        "logs/dl.{fname}.log"
    shell:
        """
        curl -s -o {output} {params.url}  &>{log}
        curl -s -o {output}.crai {params.url}.crai &>>{log}
        """

rule estimate_coverage:
    input: 
        "original_crams/{fname}.cram"
    output:
        samstats="original_crams/{fname}.cram.samstats",
        cov="original_crams/{fname}.cram.cov"
    params:
        ss_py= os.path.join(project_dir, "bin", "samStats.py"),
        cov_R= os.path.join(project_dir, "bin", "samstats2cov.R"),
        bed=os.path.join(project_dir, "beds", "gsa_24v3-0_A2.GRCh38.rand10k.bed"),
        ref=os.path.join(project_dir, "hg38", "GRCh38_full_analysis_set_plus_decoy_hla.fa")
    log:
        "logs/cov.{fname}.log"
    shell:
        """
        python3 {params.ss_py} -r {params.ref} {input} {params.bed} > {output.samstats} 2> {log}
        Rscript {params.cov_R} {output.samstats} > {output.cov} 2>> {log}
        """

def getCov(file):
    """
        Opens a "cov" file, and returns, well, the (estimated) coverage.
    """
    # apparently file is a set of strings...
    with open(file.pop()) as fh:
        line=fh.readline()
        line=fh.readline()
        s = line.rstrip().split("\t")
        return float(s[1]) # return the mean read depth...
    print("Truncated file", file, file=sys.stderr)
    exit(1)

# a bit of laziness; snakemake isn't tracking the downsampling. I'm okay with that.
rule downsample_cram:
    input:
        cram="original_crams/{fname}.cram",
        cov="original_crams/{fname}.cram.cov"
    output:
        "logs/{fname}.ds.complete"
    params:
        fname=lambda wildcards: wildcards.fname,
        ref=os.path.join(project_dir, "hg38", "GRCh38_full_analysis_set_plus_decoy_hla.fa")
    log:
        "logs/ds.{fname}.log"
    run:
        cov = getCov({input.cov})
        i=1 # random number seed for samtools view -s
        for ds in coverages:
            # added max; still makes the correct output when the required coverage > available (as in, it still gives you all of the data).
            shell("samtools view -q 20 -T {params.ref} -F 3844 -s %f -o downsampled/{params.fname}.%f.cram {input.cram} " % ( i + min(ds/cov,1), ds) )
            shell("samtools index downsampled/{params.fname}.%f.cram" % (ds))
            i += 1
        shell("touch {output}")
        # -F: excludes reads that are: unmapped, NOT primary alignment, fails QC, is dup, is supplemental alignment
        