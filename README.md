

# Quick start

## Install dependencies

### Mamba
Install mamba. If you haven't already done so, I suggest following the directions found [here](https://www.usna.edu/Users/cs/fknoll/SD211/mamba.html)
<br>
(but ignoring the ```mamba create``` command.) <br>
Install the dependencies
```
mamba install bioconda::pysam snakemake r-tidyverse r-Hmisc samtools bcftools vcftools awscli
```

### Au Naturel
Ensure the following packages are installed:
- pysam
- snakemake
- awscli (good luck)
- R
  - tidyverse
  - Hmisc
as well as samtools and bcftools (any modern version should do)
for supporting ibdgem, apparently we need VCFtools as well. boo.

## Install snp2id

clone this package. <br>
change directories to whereever you want your data to reside. **then**
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

### Take Deux!
Okay, so it was impressed on me that having relatives can be important.
As such, the original download snakemake file is now called: `src/downloadAndDownsampleOriginal.smk` (catchy, I know) <br>
I made some deep cuts to: `src/downloadAndDownsample.smk` as well.
Namely, it now grabs 4 trios in the: MXL, ASW, CEU and YRI populations. (trios meaning mom, dad, kid)
In addition, I manually added: <br>
two half siblings, followed by father of the second and shared mother <br>
`forcible["YRI"] = ["NA18913", "NA19240", "NA19239", "NA19238"]`
full siblings (first two) and parents (second two)
`forcible["MXL"] = ["NA19662", "NA19685", "NA19661", "NA19660"]`
and another quartet; full siblings (first two) and parents (second two)
`forcible["MXL"].extend( ["NA19675", "NA19680", "NA19679", "NA19678"])`

As well, I added CEPH pedigree 1463. This one was a total pain, as it's housed elsewhere. See:
```
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
```

Note that one of the samples is missing from AWS.
See [this repo](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Datasets) for a full description.
(and no, I haven't found gen 4 data at all!)

### Practicals

**Give this a day to complete**


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
  - coverage (and like information) goes here, too.
* downsampled/
  - cram files that have been downsampled.
  - e.g., ASW.NA19900.0.005000.cram
  - format:  POP . SAMPLE . COVERAGE 

## Genotyping (BCFtools)
Targeted genotyping (with bcftools) can be performed a la:
```
snakemake -c 256  -s src/genotypeBcftools.smk
```
(setting the number of cores, 256 in this case, to something sane for your system)

Targeted genotyping takes in a Panel (`panels/*.tsv.gz`) and a downsampled cram file (`downsampled/*.cram`) and creates a VCF file 
that characterizes every site in the panel. Files get written to: 
`genotypes/bcftools/{cramfile}_{panel}.vcf.gz`


### Panels
at present we just perform genotyping on autosomal loci from Kintelligence. `genotypeBcftools.smk` will create genotypes for every Panel (`panels/*.tsv.gz`) file present.
<br>
Do note that these files are in the format expected by `bcftools call`. Of note, use block gzip (`bgzip`) and index any file you wish to genotype (`tabix`).



