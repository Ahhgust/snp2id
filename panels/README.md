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



