# Panels

Here are the genotyping panels. The input files can be a bit heterogeneous (some flavor of VCF is preferred).
At a minimum, they need to be converted to a BCFtools call-specific format.

## Kintelligence
See `kintelligence_readme` for notes on how this file was made. <br>
The following was run on the Kintelligence BCF file.
```
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' kintelligence.hg38.liftover.autossnps_except6_no_unanalyzed_AND_chrXY.afannos.bcf | grep -v 'chr[XY]'| bgzip -l 9 >  kintelligence.hg38.autossnps_except6.tsv.gz  # make an autosomal snp panel (block compressed)
tabix -s1 -b2 -e2 kintelligence.hg38.autossnps_except6.tsv.gz # and index it
tabix kintelligence.hg38.autossnps_except6.tsv.gz chr22 | head # and then do a sanity test.
```

## GSA
Here's the panel file:
`GSA-24v3-0_A2.hg38.gnomadannos.autos.sites2include.bcf`
In brief, this is in hg38 (UCSC style chromosome names), with gnomad annotations (v3.1.1); just those for allele frequency. additionally this is restricted to the autosomes (`autos`), and has the "callability" mask of Woerner 2022 (removing some nasty-to-call regions of the genome)
The final callset has 586624 sites; all biallelic.
The same query/tabix postprocessing was done (see ##Kintelligence section)

## FORCE Panel
VCF (from tsv) made a la:
```
(bcftools view -h kintelligence.hg38.autossnps_except6_no_unanalyzed_AND_chrXY.afannos.anannos.bcf ; grep -vw '[XY]' FORCE_panel.txt | fgrep -vi indel |  perl -e '$_=<>; while (<>){ chomp; @s=split /\t/; print "chr$s[2]\t$s[3]\t$s[0]\t$s[5]\t$s[6]\t.\tPASS\t.\n" }') | bcftools view -o force_hg38_autosnps.bcf -Oz9
bcftools index force_hg38_autosnps.bcf
```
In words; borrow the vcf header from the kintelligence sites; remove the XY markers (as well as those that are marked as eXcluded) and write a sites-only VCF

Make the SNP call set (bcftools compatible); borrowed from kintelligence processing steps
```
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' force_hg38_autosnps.bcf |  bgzip -l 9 > force_hg38_autosnps.tsv.gz
tabix -s1 -b2 -e2 force_hg38_autosnps.tsv.gz
tabix force_hg38_autosnps.tsv.gz chr22 | head
bcftools view -H force_hg38_autosnps.bcf | perl -e 'while (<>){ chomp; @s = split /\t/; print "$s[0]\t", $s[1]-1 , "\t$s[1]\n";}' > force_hg38_autosnps.bed
```
### Allele frequencies
gnomad (v3.1.1) AF_ and AN_ annotations were added; autosomes only, and only considers sites with 5%MAF <br>
In brief, some of the GSA sites will be unannotated. Maybe we should'nt be looking at them anyways?

### cM 
I repurposed a script from Amy Williams to add cM positions to a VCF file.
The "hapmap" genetic map was used for the annotations a la:
```
bcftools view GSA-24v3-0_A2.hg38.gnomadannos.autos.sites2include_no_unanalyzed_AND_chrXY.afannos.anannos.bcf | perl ../bin/add-map-vcf.pl /eva/edatums/reference_materials/genetic_maps/hg38/Hapmap_Liftover/chr*.geneticMap | bcftools view -Oz9 -o GSA-24v3-0_A2.hg38.gnomadannos.autos.sites2include_no_unanalyzed_AND_chrXY.afannos.anannos.hapmapgmap.vcf.gz
```

