import os, sys, glob

VERSION=0.01
project_dir = os.path.normpath( os.path.join(workflow.current_basedir, "..") )
ref = os.path.join(project_dir, "hg38", "GRCh38_full_analysis_set_plus_decoy_hla.fa")

filesToRead= glob.glob( os.path.join(project_dir, "genotypes", "bcftools", "*30.00000*.vcf.gz") )
highcovs = []
for f in filesToRead:
    b = os.path.basename(f)
    # CEU.NA12414.30.0000* gets converted to just CEU.NA12414
    samp = b.split(".")[0:2]
    highcovs.append(".".join(samp) )
    
pnames="kintelligence"
 
files= glob.glob( os.path.join(project_dir, "genotypes", "bcftools", "*.vcf.gz") )
# ignore the 30x ones.
lowcovs= [ file for file in files if file.find("30.00000") ==-1 ]
# and just grab the filenames (no extension)
samps = [ os.path.basename( file[:-7] ) for file in lowcovs ]

# TODO: this can be dynamic.
panels=["kintelligence.hg38.autossnps_except6"]
callers=["bcftools"]
kpops=["CEU", "MXL", "ASW", "JPT"]
rpops=['NFE', 'AFR', "EAS", "AMR"] # gnomad population labels
#nonfinnish euro, african/aa, east asian, admixed american

rule all:
    input:
        expand("dblr/panels/{panel}.{rpop}.csv" , panel=panels, rpop=rpops),
        expand("dblr/lowpass/{prefix}/Weights_Export.xml", prefix=samps), # TODO: set to all of samps ONCE WE KNOW IT WORKS
        expand("dblr/hipass/{prefix}.{panel}.txt", prefix=highcovs, panel=panels)

def sn2s(sobj):
    """
    converts a snakemake.io object to a string. only necessary for "run" directives.
    honestly, that I have to write this function is maddening.
    it just removes {'...'} crap when you cast it to a string.
    obviously, the object must be singlular (ie, refer to just 1 file).
    had I used a shell() directive, I wouldn't need this...
    """
    return str(sobj)[2:-2]


def getHicovs(wildcards):
    return glob.glob("genotypes/bcftools/*30.0000*" + wildcards.panel + ".vcf.gz")

rule make_ibdgem_panels:
    input:
        panel="panels/{panel}_no_unanalyzed_AND_chrXY.afannos.anannos.bcf"
    output:
        "dblr/panels/{panel}.{rpop}.csv"
    params:
        rpop="{rpop}",
        binary="python3 " + os.path.join(project_dir, "bin", "gnomad2allelefreqs_long.py")
    log:
        "dblr/logs/{panel}.{rpop}.log"	
    shell: 
        """
        bcftools view -H {input.panel} | {params.binary} -p {params.rpop} > {output} 2> {log}
        """

rule make_lowpass:
    input:
        vcf="genotypes/bcftools/{prefix}.vcf.gz"
    output:
        file="dblr/lowpass/{prefix}/weights",
        xml="dblr/lowpass/{prefix}/Weights_Export.xml"
    params:
        binary="python3 " + os.path.join(project_dir, "bin", "vcf2dblr.py"),
        xbinary="python3 " + os.path.join(project_dir, "bin", "makeDblrWeights.py") 
    log:
        "dblr/logs/lowpass.{prefix}.log"	
    shell: # output is a pair; the weights file has the genotype probabilities, the XML file simply says where to find the weights file.
        """
        bcftools view {input.vcf} | {params.binary} > {output.file} 2> {log}
        {params.xbinary} `basename {input.vcf}` weights > {output.xml} 2>> {log}
        """

rule make_hipass:
    input:
        vcf="genotypes/bcftools/{prefix}.30.000000_{panel}.vcf.gz"
    output:
        "dblr/hipass/{prefix}.{panel}.txt"
    wildcard_constraints:
        prefix="|".join(highcovs) 
    params:
        binary="python3 " + os.path.join(project_dir, "bin", "vcf2dblr_ref.py")
    log:
        "dblr/logs/hipass.{prefix}.{panel}.log"	
    shell: 
        """
        bcftools view {input.vcf} | {params.binary} > {output} 2> {log}
        """

# TODO;  add dblr_long2wide

