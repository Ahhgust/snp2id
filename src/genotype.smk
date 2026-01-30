import os, sys, glob
# Written by August Woerner
# 1/30/2026; performs genotyping either using bcftools (no extra dependencies needed)
# or glimpse (with subsampling loci down to much smaller panels)
# the glimpse routines have external dependencies. in our application, not just one reference panel
# but one reference panel for each population (holding out, for example, the YRI when a YRI
# sample is imputed/refined.
# Please reach out to me if you want these panels. They're... big. too big to easily share.

VERSION=0.01
project_dir = os.path.join(workflow.current_basedir, "..")

panels=glob.glob( os.path.join(project_dir, "panels", "*.tsv.gz"))

filesToRead= glob.glob( os.path.join(project_dir, "downsampled", "*.cram") )
#filesToRead = [filesToRead[0]] # for debugging purposes...

fnames= []

for f in filesToRead:
    b = os.path.basename(f)[0:-5] # trims off the .cram; note the file format is hard coded
    fnames.append( b )

 
pnames = []
for p in panels:
    pnames.append( os.path.basename(p)[0:-7] ) # and the .tsv.gz

rule main:
    input: 
        expand("genotypes/bcftools/{fname}_{panel}.vcf.gz", fname= fnames, panel=pnames),
        expand("genotypes/glimpse2/{fname}.vcf.gz", fname= fnames),
        expand("genotypes/glimpse2_subsample/{fname}_{panel}.vcf.gz", fname= fnames, panel=pnames)
        

# runs bcftools on a specified panel of snps
rule genotype_bcftools:
    input:
        cram="downsampled/{fname}.cram",
        panel="panels/{panel}.tsv.gz"
    output:
        "genotypes/bcftools/{fname}_{panel}.vcf.gz"  
    params:
        pileup_params="-d 512 -q 20 -Q 10 -I -E -a FORMAT/DP,FORMAT/AD,FORMAT/SP",
        call_params="--ploidy 2 -Am -P 0. -C alleles -Oz9",
        bcftoolsbinary="bcftools", # pull from mamba/conda env,
        ref=os.path.join(project_dir, "hg38", "GRCh38_full_analysis_set_plus_decoy_hla.fa")
    wildcard_constraints:
        panel= "|".join(pnames)
    threads: 4
    log:
        "logs/bcftools_call.{fname}_{panel}.log"
    shell:
        """
        ({params.bcftoolsbinary} mpileup {params.pileup_params} --fasta-ref {params.ref} -T {input.panel} {input.cram} |\
        {params.bcftoolsbinary} call {params.call_params} -T {input.panel} --threads {threads} -o {output}) 2> {log} &&\
        {params.bcftoolsbinary} index --threads {threads} {output} 2>> {log}
        """

# intersects the tapir imputation panel with the panels considered. (at present, kintell and GSA; only autosomes)
rule subsample_glimpse:
    input:
        vcf="genotypes/glimpse2/{fname}.vcf.gz",
        panel="panels/{panel}.bed"
    output:
        "genotypes/glimpse2_subsample/{fname}_{panel}.vcf.gz"  
    params:
        bcftoolsbinary="bcftools" # pull from mamba/conda env,
    wildcard_constraints:
        panel= "|".join(pnames)
    log:
        "logs/glimpse2_subsample.{fname}_{panel}.log"
    shell:
        """
        {params.bcftoolsbinary} view -R {input.panel} -Oz9 -o {output} {input.vcf} &&\
        {params.bcftoolsbinary} index {output}
        """
        
# runs glimpse2 on our predefined panel of snps (TAPIR)
# note that this holds out the individual's population from the imputation panel.
rule glimpse2_holdout_population:
    input: 
        bam="downsampled/{population}.{fname}.cram"
    output:
        "genotypes/glimpse2/{population}.{fname}.vcf.gz"
    log:
        "logs/glimpse.{population}.{fname}.log"
    threads: 40 # a white lie. controls the degree of multi-processing (not threads per se).
    wildcard_constraints:
        population= "[^.]+" # no periods allowed in the population name...
    params:
        phasebinary=os.path.join(project_dir, "bin", "GLIMPSE2_phase_static"),
        phaseopt="--mapq 20",
        ligatebinary=os.path.join(project_dir, "bin", "GLIMPSE2_ligate_static"),
        ref=os.path.join(project_dir, "hg38", "GRCh38_full_analysis_set_plus_decoy_hla.fa"),
        bcftoolsbinary="bcftools", # pull from mamba/conda env,
        compressionthreads=4,
        paneldir="ImputationPanels/HoldingOutPopulations/{population}.out/reference_panel/split", # very large directory. at present, not shared; one per 1kg population
        # note that ${population}, which indicates that that population has been held out of the imputation panel.
        tmpdir="/tmp"     # notable performance boost if you make this a nvme drive...   
    shell: 
        """
        if ! outdir=`mktemp -d --tmpdir={params.tmpdir} XXXXXXglimpse2`; then
            echo "Failed to make tempdir!"
            exit 1
        fi
        
        for reference in {params.paneldir}/2.0.chr*bin; do
            bn=`basename $reference`
            chrom=`echo $bn | cut -f2 -d '_'`
            start=`echo $bn | cut -f3 -d '_'`
            stop=`echo $bn | cut -f4 -d '_'`
            echo "{params.phasebinary} {params.phaseopt} --bam-file {input.bam} -F {params.ref} --reference $reference --output $outdir/$chrom.$start.$stop.bcf "
        done | parallel -j {threads} &> {log}
        mkdir -p $outdir/ligate
        # note that curly braces are a total pain. hence the 1..22 (explicitly) below
        for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
            echo "{params.ligatebinary} --input <(ls -1v $outdir/chr$chrom.*bcf) --output $outdir/ligate/chr$chrom.bcf"
        done | parallel -j {threads} &>> {log}
        # and concatenate to make an autosomal vcf
        {params.bcftoolsbinary} concat -Oz9 --threads {params.compressionthreads} -o {output} $outdir/ligate/chr?.bcf $outdir/ligate/chr??.bcf &>> {log}
        {params.bcftoolsbinary} index --threads {params.compressionthreads} {output} &>> {log}
        rm -rf $outdir  
        """
