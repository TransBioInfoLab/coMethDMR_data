# Create EPIC_Gene_3_200.rds Data Set
# Gabriel Odom and Lissette Gomez
# 2021-05-02

library(parallel)
library(coMethDMR)
library(tidyverse)


######  Map Genes to CpGs  ####################################################
cpgGene_df <- 
	IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Other %>% 
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
# I think this is the "EPIC_CpGstoGene_min3CpGs.rds" data set.
all.equal(
	geneRegionsMin3_ls[names(EPIC_CpGstoGene_min3CpGs)],
	EPIC_CpGstoGene_min3CpGs
)
# It is, only in a different gene order. We will keep it alphabetised.
saveRDS(geneRegionsMin3_ls, "data/EPIC_CpGstoGene_min3CpGs.rds")

# NOTE: this parallel computing works in UNIX environments. For Windows, please
#   use parLapply or similar variants: 
# http://dept.stat.lsa.umich.edu/~jerrick/courses/stat701/notes/parallel.html

system.time(
	regions_min3_200bp_ls <- mclapply(
		geneRegionsMin3_ls,
		safely(coMethDMR::CloseBySingleRegion),
		arrayType = "EPIC",
		maxGap = 200,
		minCpGs = 3,
		mc.cores = 4L
	)
)
# Serial: 1.936317 min for the first 1000; 49.04342 min for full 25,578 genes
# Parallel: 1.369833 min for first 1000 over 4 cores; 19.24535 for full list

# The structure of this output is quite complicated. Within each gene, there is
#   a "result" and an "error". The "error" will often be NULL. We check how many
#   genes we lose:
regions_min3_200bp_ls %>% 
	map("error") %>% 
	map(is.null) %>% 
	as.logical() %>% 
	`!` %>% 
	which()
regions_min3_200bp_ls[13]

# # Now, we can look inside the "results" section: it will be a list for each
# #   cluster within the gene, so we want to unlist and unname *within* each gene
# geneClusters_min3_200bp_ls <- 
# 	regions_min3_200bp_ls %>% 
# 	map("result") %>% 
# 	map(unlist) %>% 
# 	map(unname) %>% 
# 	compact()

# All clusters within genes (considered independently)
allClusters_min3_200bp_ls <- 
	regions_min3_200bp_ls %>% 
	map("result") %>% 
	unlist(recursive = FALSE) %>% 
	unique()



######  Name Regions  #########################################################
# system.time(
# 	geneClustersDF_ls <- mclapply(
# 		geneClusters_min3_200bp_ls,
# 		coMethDMR::OrderCpGsByLocation,
# 		arrayType = "EPIC",
# 		output = "dataframe",
# 		mc.cores = 4L
# 	)
# )
# # Serial: 38.8541 min for all 25578
# # Parallel: 15.09952 min
# 
# names(geneClusters_min3_200bp_ls) <- map(geneClustersDF_ls, NameRegion)
# # This doesn't match any of the data sets that are included in coMethDMR at this
# #   time. I'm not sure what this code was supposed to do.


system.time(
	allClustersDF_ls <- mclapply(
		allClusters_min3_200bp_ls,
		coMethDMR::OrderCpGsByLocation,
		arrayType = "EPIC",
		output = "dataframe",
		mc.cores = 4L
	)
)
# Parallel: 34.4306 min for 45,351 clusters

names(allClusters_min3_200bp_ls) <- map(allClustersDF_ls, NameRegion)
# We believe this is the updated version of EPIC_Gene_3_200.rds; it has 2152
#   fewer clusters (1 less MB), but this could be due to changes in the Illumina
#   data set since 2018 (we are re-creating this data set in 2021, the original
#   was created three years ago).
# Because I can't re-create the old version exactly, I will include them both.
saveRDS(allClusters_min3_200bp_ls, "data/EPIC_Gene_3_200_new.rds")
