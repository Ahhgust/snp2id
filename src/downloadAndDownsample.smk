import os, sys, random

VERSION=0.02
project_dir = os.path.join(workflow.current_basedir, "..")

# select 5 kids from 4 populations
# at present, these are just from 30x recall of the 1000 Genomes Project
# (we have metadata for hgdp as well)
nsamps=5
kpops=["CEU", "MXL", "ASW", "YRI"] # Modified; includes YRI instead of JPT (for trios)
# and downsample to the specified coverages
# modified to include 30x (high pass)
coverages=[0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.50, 1.0, 5., 10., 30.]

# and stuff needed to relate the above to
# actually downloading and downsampling
url_lut={} # 2d; pop -> sample -> url
url_flat={} # 1D pop/sampleid -> url
# okay. fun. relatives can be either in the unrelated set, or the related set from the 1000 Genomes.
# (this makes sense, as you'd want the (hopefully unrealted) parents but not the kid in a trio)
kgfile1=os.path.join(project_dir, "pop_metadata", "1000G_698_related_high_coverage.sequence.index")
kgfile2 = os.path.join(project_dir, "pop_metadata", "1kg.urls.all")
URLCOL=0
POPCOL=10
SCOL=9

# let's forcibly include some samples
#
forcible = {}
# two half siblings, followed by father of the second and shared mother
forcible["YRI"] = ["NA18913", "NA19240", "NA19239", "NA19238"]
#  full siblings (first two) and parents (second two)
forcible["MXL"] = ["NA19662", "NA19685", "NA19661", "NA19660"]
# and another quartet; full siblings (first two) and parents (second two)
forcible["MXL"].extend( ["NA19675", "NA19680", "NA19679", "NA19678"] )


random.seed(10)
samp2pop={}
# parse the 1KG data
for kgfile in (kgfile1, kgfile2):
    with open(kgfile) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            sp = line.rstrip().split("\t")
            #annoyingly, some identifiers have leading/trailing whitespace
            (u, p, s) = (sp[URLCOL].strip(), sp[POPCOL].strip(), sp[SCOL].strip())
            #print("'", u, "' '", p, "' '",s, "',", sep="")
            samp2pop[s]=p
            if p not in kpops:
                continue
            if p not in url_lut:
                url_lut[p]=dict()
            url_lut[p][s] = u

# child id -> dad id, mom id
kids2parents= {}
# population -> list of kids
popkids = {}
pedfile = os.path.join(project_dir, "pop_metadata", "1kGP.3202_samples.pedigree_info.txt")

with open(pedfile) as fh:
    line = fh.readline() # byebye header
    for line in fh:
        s = line.rstrip().split(" ")
        
        if len(s) < 4 or s[0] not in samp2pop:
            continue
            
        samp = s[0]
        pop = samp2pop[samp]
        if pop not in kpops:
            continue    
        if s[1] == '0' or s[2] == '0': # need both parents
            continue
        
        if pop not in popkids:
            popkids[pop] = []
           
        kids2parents[samp] = [ s[1], s[2] ]
        popkids[pop].append(samp)  
    

for pop in kpops:
    if pop not in url_lut:
        print("Missing population: ", pop, list(url_lut.keys()), file=sys.stderr)
        exit(1)
        
    inner = url_lut[ pop ]
    
    if pop not in popkids:
        print("Pop: ", pop, " has no trios...", file=sys.stderr)
        exit(1)
    
    # list of kid sample ids from this pop
    kids = list(popkids[pop])
    #print(len(kids), pop)
    
    if len(kids) < nsamps:
        print("Not enough samples available for", pop, len(kids), kids, file=sys.stderr)
        exit(1)
    
    
    random.shuffle(kids)
    toinclude = set()
 
    for i in range(nsamps):
        toinclude.add(kids[i])
        toinclude.add(kids2parents[kids[i]][0])
        toinclude.add(kids2parents[kids[i]][1])
    
    
    for s, u in inner.items():
        if s in toinclude:
            #1
            url_flat[ pop + "." + s  ] = u


for p, samps in forcible.items():
    for s in samps:
        if s not in url_lut[p]:
            print("Unexpected; missing ", s , " from ", p, file=sys.stderr)
            continue
        url_flat[ p + "." + s ] = url_lut[p][s]


# manually add pedigree 1463 as well.
# See: https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Datasets
# note: need awscli (to download) AND these are bam files. because of course they are.
# also note that they used ... noncanonical names for the bam files. because of course.
url_flat["CEU.NA12889"] = "s3://platinum-pedigree-data/data/illumina/mapped/GRCh38/2281.v4.2.4.grc38.bam"
url_flat["CEU.NA12890"] = "s3://platinum-pedigree-data/data/illumina/mapped/GRCh38/2280.v4.2.4.grc38.bam" # current MIA. issue raised on github repo
url_flat["CEU.NA12891"] = "s3://platinum-pedigree-data/data/illumina/mapped/GRCh38/2214.v4.2.4.grc38.bam"
url_flat["CEU.NA12892"] = "s3://platinum-pedigree-data/data/illumina/mapped/GRCh38/2213.v4.2.4.grc38.bam"
url_flat["CEU.NA12877"] = "s3://platinum-pedigree-data/data/illumina/mapped/GRCh38/2209.v4.2.4.grc38.bam"
url_flat["CEU.NA12878"] = "s3://platinum-pedigree-data/data/illumina/mapped/GRCh38/2188.v4.2.4.grc38.bam"
url_flat["CEU.NA12879"] = "s3://platinum-pedigree-data/data/illumina/mapped/GRCh38/2216.v4.2.4.grc38.bam"
url_flat["CEU.NA12881"] = "s3://platinum-pedigree-data/data/illumina/mapped/GRCh38/2211.v4.2.4.grc38.bam"
url_flat["CEU.NA12882"] = "s3://platinum-pedigree-data/data/illumina/mapped/GRCh38/2212.v4.2.4.grc38.bam"
url_flat["CEU.NA12885"] = "s3://platinum-pedigree-data/data/illumina/mapped/GRCh38/2217.v4.2.4.grc38.bam"
url_flat["CEU.NA12886"] = "s3://platinum-pedigree-data/data/illumina/mapped/GRCh38/2189.v4.2.4.grc38.bam"

platinumSamples = set( ["CEU.NA12889", "CEU.NA12890",  "CEU.NA12891", "CEU.NA12892", \
"CEU.NA12877", "CEU.NA12878", "CEU.NA12879", "CEU.NA12881", "CEU.NA12882", "CEU.NA12885", "CEU.NA12886"] ) 

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
        url = lambda wildcards: url_flat[ wildcards.fname ],
        ref=os.path.join(project_dir, "hg38", "GRCh38_full_analysis_set_plus_decoy_hla.fa")        
    threads:
        1
    log:
        "logs/dl.{fname}.log"
    run:
        # download BAM from platinum pedigree website; convert to cram. delete bam
        if wildcards.fname in platinumSamples:
            shell("""
            aws s3 cp --no-sign-request {params.url} {output}.bam 
            samtools view -T {params.ref} --write-index -C -o {output} {output}.bam
            rm {output}.bam
            """)
        else: # download CRAM from the 1kGP website
            shell("curl -s -o {output} {params.url}  &>{log}")
            shell("curl -s -o {output}.crai {params.url}.crai &>>{log}")

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
        