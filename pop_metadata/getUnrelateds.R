# no shebang. sorry!
# written by August Woerner
# does some easy-peasy preprocessing to make two files; one of unrelated individuals (parents in trios, for examples)
# and another as kids
# this is a sanity check to ensure that every sample maps into one or the other


suppressPackageStartupMessages(library(tidyverse))

# Note this is just laziness. the routine only works w/ the 1000 Genomes data.
files <- Sys.glob("../genotypes/bcftools/*30*GSA*gz")
tibble(File=basename(files),
	   SampleID=substr(File, 5, 11)) |>
	   select(SampleID) |> distinct() -> sampsWeHave


plat <- read_tsv("platinumpedigree.info", col_types=cols())

# definitions are always squishy (relatives are defined on non-relatives and vice versa. eg, I could kick some samples out, and some relatives would become unrelateds)
# let's just be naive about things and take the grandparents
filter(plat, Generation=='G1') |> 
   select(SampleID=primary.id) |> 
   inner_join(sampsWeHave, by="SampleID")   -> platunrel

# note if we ever get data from the spouses of the G4 generation, this should be updated.   
filter(plat, Generation!='G1') |> 
   select(SampleID=primary.id) |>
   inner_join(sampsWeHave, by="SampleID")    -> platrel

ped <- read_delim("1kGP.3202_samples.pedigree_info.txt", delim=" ", col_types=cols())

rename(ped, SampleID=sampleID)-> ped

filter(ped, fatherID==0&motherID==0) -> parents

filter(ped, !(fatherID==0&motherID==0)) -> kiddos # includes half-sibs

inner_join(sampsWeHave, kiddos, by='SampleID') |>
    select(SampleID) |>
	bind_rows(platunrel) -> allkids

# some trios randomly selected also happened to hit on the plat pedigree. 
inner_join(sampsWeHave, parents, by='SampleID') |>
    select(SampleID) |>
	bind_rows(platrel) |>
	anti_join(allkids, by="SampleID") -> allparents

	
# should be 0 rows
filter(sampsWeHave, ! (SampleID %in% allkids$SampleID),
					! (SampleID %in% allparents$SampleID))

allparents |>
	write_tsv("Parents.tsv")
	
allkids|>
	write_tsv("Kids.tsv")