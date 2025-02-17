import os, sys, glob

VERSION=0.01
project_dir = os.path.join(workflow.current_basedir, "..")

panels=glob.glob( os.path.join(project_dir, "panels", "*.tsv.gz"))

filesToRead= glob.glob( os.path.join(project_dir, "downsampled", "*.cram") )

fnames= []

for f in filesToRead:
    b = os.path.basename(f)[0:-5] # trims off the .cram; note the file format is hard coded
    fnames.append( b )

 
pnames = []
for p in panels:
    pnames.append( os.path.basename(p)[0:-7] ) # and the .tsv.gz

rule main:
    input: 
        expand("genotypes/bcftools/{fname}_{panel}.vcf.gz", fname= fnames, panel=pnames)


rule genotype_bcftools:
    input:
        cram="downsampled/{fname}.cram",
        panel="panels/{panel}.tsv.gz"
    output:
        "genotypes/bcftools/{fname}_{panel}.vcf.gz"  
    params:
        pileup_params="-d 512 -q 20 -Q 10 -I -E -a FORMAT/DP,FORMAT/AD,FORMAT/SP",
        call_params="--ploidy 2 -Am -P 0. -C alleles -Oz9",
        ref=os.path.join(project_dir, "hg38", "GRCh38_full_analysis_set_plus_decoy_hla.fa")
    wildcard_constraints:
        panel= "|".join(pnames)
    threads: 4
    log:
        "logs/bcftools_call.{fname}_{panel}.log"
    shell:
        """
        (bcftools mpileup {params.pileup_params} --fasta-ref {params.ref} -T {input.panel} {input.cram} |\
        bcftools call {params.call_params} -T {input.panel} --threads {threads} -o {output}) 2> {log} &&\
        bcftools index --threads {threads} {output} 2>> {log}
        """

