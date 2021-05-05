# Create 450k_CpGstoGene_min3CpGs.rds Data Set
# Gabriel Odom and Lissette Gomez
# 2021-05-02

library(parallel)
library(coMethDMR)
library(tidyverse)


######  Map Genes to CpGs  ####################################################
cpgGene_df <- 
	IlluminaHumanMethylation450kanno.ilmn12.hg19::Other %>% 
	subset(select = "UCSC_RefGene_Name") %>% 
	as.data.frame() %>% 
	rownames_to_column(var = "CpG") %>% 
	as_tibble() %>% 
	filter(UCSC_RefGene_Name != "") %>% 
	filter(str_detect(CpG, "cg")) %>% 
	separate_rows(UCSC_RefGene_Name, sep = ";") %>% 
	unique()

allGeneRegions_ls <- 
	cpgGene_df %>% 
	split(f = cpgGene_df$UCSC_RefGene_Name) %>% 
	map(pull, "CpG")



######  Find Regions with Clusters of CpGs  ###################################
geneRegionsMin3_ls <- allGeneRegions_ls[lengths(allGeneRegions_ls) >= 3]
# This is the "450k_CpGstoGene_min3CpGs.rds" data set.
# NOTE: for some stupid reason, we can create a data file that starts with a 
#   number, but we can't load it into the environment with that name. I've added
#   the "x" to the start of the name.
all.equal(
	geneRegionsMin3_ls[names(x450k_CpGstoGene_min3CpGs)],
	x450k_CpGstoGene_min3CpGs
)
# It is, only in a different gene order. We will keep it alphabetised.
saveRDS(geneRegionsMin3_ls, "data/450k_CpGstoGene_min3CpGs.rds")

