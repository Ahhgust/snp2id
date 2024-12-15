

# Quick start

## Install dependencies

### Mamba
Install mamba. If you haven't already done so, I suggest following the directions found [here](https://www.usna.edu/Users/cs/fknoll/SD211/mamba.html)
<br>
(but ignoring the ```mamba create``` command. <br>
Install the dependencies
```
mamba install bioconda::pysam snakemake r-tidyverse r-Hmisc samtools bcftools
```

### Au Naturel
Ensure the following packages are installed:
- pysam
- snakemake
- R
  - tidyverse
  - Hmisc
as well as samtools and bcftools (any modern version should do)

## Install snp2id

###
clone this package. <br>
change directories to whereever you want your data to reside. **then **
```
git clone https://github.com/Ahhgust/snp2id.git
```

## Download the hg38 reference genome
Okay, so there are a lot of "reference genomes". We want the one used by the 1000 Genomes project.
```
cd snp2id
mkdir hg38 downsampled logs original_crams # not entirely necessary, but helpful to see
cd hg38
for suffix in fa dict fa.alt fa.amb fa.ann fa.bwt fa.fai fa.pac fa.sa
do
	wget "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.$suffix"
done
cd ..
```

**this will take a while**

alternatively, if you have the hla-masked hg38 reference genome, you can just add /symlink as appropriate.

### Downloading and downsampling

Download 5 (arbitrary) individuals from 4 populations ("CEU", "MXL", "ASW", "JPT"]). These are downloaded into ```original_crams``` <br>
(okay, not quite arbitrary; we only consider sample IDs from the "clean pcs" dataset of Koenig et al. 2024) <br>
And then downsample them to the following (arbitrary) coverages: [0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.50, 1.0, 5., 10.]
<br><br>
We want to be good people-- so let's limit the concurrency, both so the 1KG folks don't get mad at as, and so we don't blow out the IO when we downsample.
<br>
First, let's try a dry run
```
snakemake -n -s src/downloadAndDownsample.smk --core 5
```
and if we don't note any errors, do:
```
nohup snakemake -s src/downloadAndDownsample.smk --core 5 &>  dlAndds.outerr &
```

The directory structure (once complete) is:
* beds/ # genomic data in bed-file format
  - gsa_24v3-0_A2.GRCh38.rand10k.bed # 10,000 autosomal SNP positions; used to estimate coverage
* src/
  - various snakemake routines go here. For now, just a script to download and downsample is used
  - TODO: a snakemake script to genotype
* bin/
  - samStats.py # provides site-level summary statistics (verbose)
  - samstats2cov.R  # and summarizes these statistics
* logs/
  - *individual log-files go here*
* pop_metadata/
  - See README; information on what sample-information is available where.
* original_crams
  - the cram files, as downloaded from the 1KG project
* downsampled/
  - cram files that have been downsampled.
  - e.g., ASW.NA19900.0.005000.cram
  - format:  POP . SAMPLE . COVERAGE 



