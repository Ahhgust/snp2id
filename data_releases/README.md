# Data releases
Just to make life a bit easier, we also provide some "canned" results from the pipeline.
Some data are small enough to share (directly), while others are given as dropbox links because they're big. <br>
Note that all results are **unfiltered**. If you're using BCFtools, and using the PL tag in particular, I wouldn't worry, but the called genotypes (in general) will not be so great (especially <=5x). For GLIMPSE the concerns are a bit more subtle. Namely, sites with high posterior probability (GP~1) may simply reflect a large prior. I doubt such sites are what you want, especially as coverage wanes. I (personally) would recommend filtering on a Bayes factor (BF); one perhaps appropriate way is to filter on the posterior odds / the prior odds (naive, assuming Hardy Weinberg equilibrium), though there certainly are other ways to do it.

## GSA
The GSA files are too big to (easily) post here. <br>
Please download the GSA files from: <br>
https://www.dropbox.com/scl/fi/9ej04ejxnzebe4veyg71r/Downsampling_Bcftools_GSA.tar?rlkey=72vm92dlifsn5io564sap3s4u&dl=1
<br>
Or click the [link](https://www.dropbox.com/scl/fi/9ej04ejxnzebe4veyg71r/Downsampling_Bcftools_GSA.tar?rlkey=72vm92dlifsn5io564sap3s4u&dl=1)

## GLIMPSE
And the GLIMPSE (v2.0.0) files are too big to share as well. <br>
You can click onthe following [dropbox link](https://www.dropbox.com/scl/fi/0loa78rytg0w7wpy4lky8/Glimpse_Tapirpanel.tar?rlkey=97rn55j2n1hxdo6w35bnfaca0&dl=1) <br>
to get the refined/imputed genotype calls. Note that the panel in this case is essentially the 1000 Genomes Project (+HGDP). It is ~30M sites that are also used in Tapir (doi: 10.1016/j.fsigen.2025.103387). <br>
I also provide just the GSA (as well as just kintelligence) calls [here](https://www.dropbox.com/scl/fi/wkjudodb9bf83j7gyc6uc/Glimpse_Tapirpanel_gsa_and_kintell.tar?rlkey=onwp8a9syiy0d77qvbjxtztxe&dl=1) <br>
strictly speaking, these files are the intersection of the GSA/kintell with the sites from Tapir. <br>

### Method
GLIMPSE was run as it is in Tapir with one small modification. Namely, the true population was held out of the imputation panel. The imputation panels themselves are essentially the same (e.g., the same "chunks" are used). I would add that key summary statistics will also change. For example, the RAF (alternative allele frequency from the reference panel) will now differ (for the same site) depending on the population considered. 

## Version
Alpha