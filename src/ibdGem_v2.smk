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
    
print(highcovs)

# BCF dir path.
# Taken from Koenig et al 2024 (1KGP+ hgdp, relatives and unrelated removed)
#/eva/edatums/reference_materials/imputation/GnomadPanelV2/Raw/1kg_hgdp_phased_haplotypes_v2_unrelateds_cleanpcs/hgdp1kgp_chr22.filtered.SNV_INDEL.phased.shapeit5.unrelateds.bcf
kgpath="/eva/edatums/reference_materials/imputation/GnomadPanelV2/Raw/1kg_hgdp_phased_haplotypes_v2_unrelateds_cleanpcs"

# UCSC chromosome naming.
chroms= [ "chr" + str(i) for i in range(1,23) ]
pnames="kintelligence"
 
# input crams
crams=glob.glob("downsampled/*cram")
# ignore the 30x ones.
lowcovs= [ file for file in crams if file.find("30.00000") ==-1 ]
# and just grab the filenames (no extension)
samps = [ os.path.basename( file[:-5] ) for file in lowcovs ]

# TODO: this can be dynamic.
panels=["kintelligence.hg38.autossnps_except6"]
callers=["bcftools"]
kpops=["CEU", "MXL", "ASW", "JPT"]


rule all:
    input:
        expand("ibdgem_v2/panels/{chrom}_{panel}.impute.hap" , chrom=chroms, panel=panels),
        expand("ibdgem_v2/pops/{population}" , population=kpops),
        expand("ibdgem_v2/{population}/{lowcov}_{panel}/{chrom}/out.txt", population=kpops, chrom=chroms, panel=panels, lowcov=samps)
        
        

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

rule merge_hicovs:
    input:
        hicovs=getHicovs,
        panel="panels/{panel}.tsv.gz"
    output:
        "genotypes/bcftools/hicovs_{panel}.vcf.gz"
    log:
        "ibdgem_v2/logs/hicov.{panel}.log"	
    shell:
        """      
        bcftools merge -R {input.panel} -Oz9 --threads 4 -o {output} {input.hicovs} && bcftools index {output}
        """

rule make_ibdgem_panels:
    input:
        panelexclude="panels/exclude",
        panel="panels/{panel}.tsv.gz",
        hicov="genotypes/bcftools/hicovs_{panel}.vcf.gz"
    output:
        "ibdgem_v2/panels/{chrom}_{panel}.impute.hap"
    wildcard_constraints:
        chrom = "|".join(chroms)
    params:
        filt="--min-alleles 2 --max-missing 1 --IMPUTE",
        kgpath=kgpath
    log:
        "ibdgem_v2/logs/panel.{chrom}.{panel}.log"	
    shell: #notes: the "hicov" sample IDs will be fine; there will be "2:sampleid" versions of them as well. see notes on bcftools merge.
           # I consider the 2: versions to be a "feature" as we can separate the effects of how you estimate genotypes from the samples themselves.
        """
        #removeme=`cat {input.panelexclude} | while read indv; do echo -n "--remove-indv $indv "; done` # indv to remove from the pop dataset
        #bcf="{params.kgpath}/hgdp1kgp_{wildcards.chrom}.filtered.SNV_INDEL.phased.shapeit5.unrelateds.bcf" # the pop dataset  
        bcftools view -r {wildcards.chrom} {input.hicov} | sed "/^##/! s/\//|/g" | bcftools view -o ibdgem_v2/{wildcards.chrom}.vcf.gz 2> {log} # convert to quasi-phased genotypes, as per ibdgem manual
        bcftools index ibdgem_v2/{wildcards.chrom}.vcf.gz 2>> {log}
        bcftools merge --force-samples ibdgem_v2/{wildcards.chrom}.vcf.gz {params.kgpath}/hgdp1kgp_{wildcards.chrom}.filtered.SNV_INDEL.phased.shapeit5.unrelateds.bcf  |\
        vcftools {params.filt} --vcf - --positions <(zgrep -w {wildcards.chrom} {input.panel}| cut -f1,2) --out ibdgem_v2/panels/{wildcards.chrom}_{wildcards.panel} 2>> {log}
        rm ibdgem_v2/{wildcards.chrom}.vcf.gz ibdgem_v2/{wildcards.chrom}.vcf.gz.???
        """

# generates pileup "files" (plain text) 
# of data associated w/ a given snp panel for a particular individual.
# the whole panel is generated. subsample to a particular chrom as needed (use -c in ibdgem to constrain things)
rule make_pileups:
    input:
        cram="downsampled/{sprefix}.cram",
        panel="panels/{pprefix}.bed"
    output:
        "ibdgem/pileups/{sprefix}.{pprefix}.pile"
    wildcard_constraints:
        pprefix = "|".join(panels)
    log:
        "ibdgem/logs/{sprefix}_{pprefix}.log"
    params: 
        popt="--output-MQ -q 20 -Q 10", # add mapq, min map quality of 20, min bq of 10,
        reference=ref
    shell:
        "samtools mpileup {params.popt} --reference {params.reference} -l {input.panel} {input.cram} > {output} 2> {log}"


# generates a file that has all of the sample IDs associated with some population
# less those sample IDs we're using for HID.
rule make_pop_panels:
    input:
        panelexclude="panels/exclude",
        meta="pop_metadata/gnomad_meta_updated.tsv",
        cleanpcs="pop_metadata/kghgdp_cleanpc_samples.txt"
    output:
        "ibdgem_v2/pops/{population}"
    log:
        "ibdgem_v2/logs/pops.{population}"	
    wildcard_constraints:
        population = "|".join(kpops)    
    shell:
        """
        cut -f1,188 {input.meta} | fgrep -w {wildcards.population} | cut -f1 | fgrep -w -v -f {input.panelexclude} | fgrep -w -f {input.cleanpcs} > {output} 2> {log}
        """
        
rule ibdgem:
    input:
        panel="ibdgem_v2/panels/{chrom}_{panel}.impute.hap",
        population="ibdgem_v2/pops/{population}",
        pile="ibdgem/pileups/{lowcov}.{panel}.pile" # note the lack of v2 (this part is samsies)
    output:
        "ibdgem_v2/{population}/{lowcov}_{panel}/{chrom}/out.txt" # a little stickier than it looks
    wildcard_constraints:
        panel = "|".join(panels),
        population = "|".join(kpops) 
    params:
        ibdgem=os.path.join(project_dir, "bin", "ibdgem"),
        fcat="perl " + os.path.join(project_dir, "bin", "fcat.pl")
    log:
        "ibdgem_v2/logs/{lowcov}_{population}_{chrom}_{panel}.denom.log"
    run:
        #bin/ibdgem -P ibdgem/pileups/CEU.NA06994.0.001000.kintelligence.hg38.autossnps_except6.pile  -c chr2 -H ibdgem_v2/panels/chr2_kintelligence.hg38.autossnps_except6.impute.hap -L ibdgem_v2/panels/chr2_kintelligence.hg38.autossnps_except6.impute.legend -I ibdgem_v2/panels/chr2_kintelligence.hg38.autossnps_except6.impute.hap.indv \
        #-B ibdgem_v2/pops/CEU -S panels/exclude -r -O TEST
        i=input.panel[0:-4] #trims off .hap
        legend=i + ".legend" 
        indv=i+ ".hap.indv"
        outdir=str(output)[0:-8] # trims off /out.txt
        # ibdgem will crash iff there are 0 snps for a given chromosome. 
        # as this command is done per individual per chromosome (and not per individual, as per the numerator in the LR
        # this rule will sometimes fail. hence the tolerance for failure (pipefile command)
        shell("set +o pipefail && mkdir -p {outdir} && touch {output} && \
            {params.ibdgem} -H {input.panel} -L {legend} -I {indv} -c {wildcards.chrom} -P {input.pile} --LD -O {outdir} -B {input.population} -S panels/exclude -v 2> {log} && \
            {params.fcat} -h {outdir}/*summary.txt > {output} 2>> {log} || true")