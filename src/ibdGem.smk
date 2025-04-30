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
print(samps)

rule all:
    input:
        expand("ibdgem/panels/{chrom}_{panel}.impute.hap" , chrom=chroms, panel=panels),
        expand("ibdgem/pileups/{sprefix}.{pprefix}.pile", sprefix=samps, pprefix=panels),
        expand("ibdgem/hicov/{hicov}_{lowcov}_{panel}/out.txt", hicov=highcovs, panel=panels, lowcov=samps),
        expand("ibdgem/denom/{lowcov}_{panel}/{chrom}/out.txt", lowcov=samps,panel=panels, chrom=chroms)
        
        

def sn2s(sobj):
    """
    converts a snakemake.io object to a string. only necessary for "run" directives.
    honestly, that I have to write this function is maddening.
    it just removes {'...'} crap when you cast it to a string.
    obviously, the object must be singlular (ie, refer to just 1 file).
    had I used a shell() directive, I wouldn't need this...
    """
    return str(sobj)[2:-2]


rule make_ibdgem_panels:
    input:
        panelexclude="panels/exclude",
        panel="panels/{panel}.tsv.gz"
    output:
        "ibdgem/panels/{chrom}_{panel}.impute.hap"
    wildcard_constraints:
        chrom = "|".join(chroms)
    params:
        filt="--min-alleles 2 --max-missing 1 --IMPUTE",
        kgpath=kgpath
    log:
        "ibdgem/logs/panel.{chrom}.{panel}.log"	
    run:
        excl=""
        with open(input.panelexclude, "r") as infile:
            for line in infile:
                excl += " --remove-indv " + line.rstrip()
        out_prefix=str({output[0]})[2:-13] # okay, "output" is of type snakemake.io.OutputFiles; converting it to a string adds {} crap to it. even though it is represented just fine as a string naturally
        # we need to trim off the file suffix (impute.hap), and this is a correct (but ugly) way to do that.
        
        bcf=os.path.join(sn2s({params.kgpath}), "hgdp1kgp_" + sn2s({wildcards.chrom}) + ".filtered.SNV_INDEL.phased.shapeit5.unrelateds.bcf")
        shell("vcftools {params.filt} {excl} --bcf {bcf} --positions <(zgrep -w {wildcards.chrom} {input.panel}| cut -f1,2) --out {out_prefix} 2> {log}")


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


# files for the numerator in the LR
# they involve 1 high coverage sample (30x)
# compared to some data.
# (the contributor hypothesis may be true or false; 
# ie, the highcov sample may or may not be the same as the lowcov)
rule hi_cov_ibdgem:
    input:
        vcf="genotypes/bcftools/{hicov}.30.000000_{panel}.vcf.gz",
        pile="ibdgem/pileups/{lowcov}.{panel}.pile"
    output:
        "ibdgem/hicov/{hicov}_{lowcov}_{panel}/out.txt" # a little stickier than it looks
    wildcard_constraints:
        panel = "|".join(panels)
    params:
        ibdgem=os.path.join(project_dir, "bin", "ibdgem"),
        fcat="perl " + os.path.join(project_dir, "bin", "fcat.pl")
    log:
        "ibdgem/logs/hicov/{hicov}_{lowcov}_{panel}.log"
    shell:
        """
        if ! outdir=`mktemp -d --tmpdir=/tmp XXXXXXibdgem`; then
            echo "Failed to make tempdir!"
            exit 1
        fi
        rm -f {output}
        touch {output}
        
        dname=`dirname {output}`
        mkdir -p $dname
        for chrom in chr1  chr2  chr3  chr4  chr5  chr6  chr7  chr8  chr9  chr10  chr11  chr12  chr13  chr14  chr15  chr16  chr17  chr18  chr19  chr20  chr21  chr22;
        do
          if [[ $(fgrep -c -w $chrom {input.pile}) -ne 0 ]]; then #ibdgem crashes if there are no snps found on the chromosome
             mkdir -p $outdir/$chrom
             {params.ibdgem} -V {input.vcf} -P {input.pile} -c $chrom --LD -O $outdir/$chrom 2>> {log}
             {params.fcat} -h $outdir/$chrom/*summary.txt | tail -n +2 >> {output}
             cat $outdir/$chrom/*summary.txt > $dname/$chrom.summary
             fgrep -v '#' $outdir/$chrom/*tab.txt > $dname/$chrom.tab
          fi
        done
        rm -rf $outdir
        """
        
rule panel_ibdgem:
    input:
        panel="ibdgem/panels/{chrom}_{panel}.impute.hap",
        pile="ibdgem/pileups/{lowcov}.{panel}.pile"
    output:
        "ibdgem/denom/{lowcov}_{panel}/{chrom}/out.txt" # a little stickier than it looks
    wildcard_constraints:
        panel = "|".join(panels)
    params:
        ibdgem=os.path.join(project_dir, "bin", "ibdgem"),
        fcat="perl " + os.path.join(project_dir, "bin", "fcat.pl")
    log:
        "ibdgem/logs/{lowcov}_{chrom}_{panel}.denom.log"
    run:
        i=input.panel[0:-4] #trims off .hap
        legend=i + ".legend" 
        indv=i+ ".hap.indv"
        outdir=str(output)[0:-8] # trims off /out.txt


        # ibdgem will crash iff there are 0 snps for a given chromosome. 
        # as this command is done per individual per chromosome (and not per individual, as per the numerator in the LR
        # this rule will sometimes fail. hence the tolerance for failure (pipefile command)
        shell("set +o pipefail && mkdir -p {outdir} && touch {output} && \
            {params.ibdgem} -H {input.panel} -L {legend} -I {indv} -c {wildcards.chrom} -P {input.pile} --LD -O {outdir} 2> {log} && \
            {params.fcat} -h {outdir}/*summary.txt > {output} 2>> {log} || true")