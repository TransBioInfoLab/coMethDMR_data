# Create EPIC_CpGstoGene_min3CpGs.rds & EPIC_Gene_3_200.rds Data Sets
# Gabriel Odom and Lissette Gomez
# 2021-05-02
# UPDATED: 2022-02-22


library(parallel)
library(coMethDMR)
library(tidyverse)



######  Map Genes to CpGs  ####################################################
cpgGene_df <- 
  IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other %>% 
	# IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Other %>% 
	subset(select = "UCSC_RefGene_Name") %>% 
	as.data.frame() %>% 
	rownames_to_column(var = "CpG") %>% 
	as_tibble() %>% # 865,589
	filter(UCSC_RefGene_Name != "") %>% # 616,598
	filter(str_detect(CpG, "cg")) %>% # 615,443
	separate_rows(UCSC_RefGene_Name, sep = ";") %>% # 1,357,524
	unique() # 683,797

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
saveRDS(
  geneRegionsMin3_ls,
  "data/EPIC_10b4_CpGstoGene_min3CpGs.rds"
)
geneRegionsMin3_ls <- readRDS(
  "data/EPIC_10b4_CpGstoGene_min3CpGs.rds"
)


# NOTE: this parallel computing works in UNIX environments. For Windows, please
#   use parLapply or similar variants: 
# http://dept.stat.lsa.umich.edu/~jerrick/courses/stat701/notes/parallel.html

# Mac
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

# Windows

system.time(
  regions_min3_200bp_ls <- parLapply(
    cl = makeCluster(3L),
    X = geneRegionsMin3_ls,
    fun = safely(coMethDMR::CloseBySingleRegion),
    arrayType = "EPIC",
    maxGap = 200,
    minCpGs = 3
  )
)
# Serial: 43.71 sec for the first 100; 49.04342 min for full 25,566 genes
# Parallel: 40.31 sec for first 100 over 4 cores
# Parallel: 60.11 sec for first 100 over 8 cores
# Parallel: 37.94 sec for first 100 over 2 cores
# Parallel: 37.93 sec for first 100 over 3 cores; 27.08433 for full list


# The structure of this output is quite complicated. Within each gene, there is
#   a "result" and an "error". The "error" will often be NULL. We check how many
#   genes we lose:
regions_min3_200bp_ls %>% 
	map("error") %>% 
	map(is.null) %>% 
	as.logical() %>% 
	`!` %>% 
	which()
# UPDATE 2022-02-22: we lose none.
# regions_min3_200bp_ls[13]

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
	map("result") %>% # 25,566
	unlist(recursive = FALSE) %>% # 53,224
	unique() #48,564


saveRDS(
  allClusters_min3_200bp_ls,
  "data/EPIC_10b4_Gene_3_200_unnamed.rds"
)
allClusters_min3_200bp_ls <- readRDS(
  "data/EPIC_10b4_Gene_3_200_unnamed.rds"
)



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


# Mac
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

# Windows
system.time(
  allClustersDF_ls <- parLapply(
    cl = makeCluster(2L),
    X = allClusters_min3_200bp_ls,
    fun = coMethDMR::OrderCpGsByLocation,
    arrayType = "EPIC",
    output = "dataframe"
  )
)
# Parallel: 36.97 sec for first 100 over 2 cores; 62.92 min for 48,564 clusters
# Parallel: 38.19 sec for first 100 over 3 cores

saveRDS(allClustersDF_ls, "data/EPIC_10b4_3_200_ordered_clusters.rds")

names(allClusters_min3_200bp_ls) <- map(allClustersDF_ls, NameRegion)
# We believe this is the updated version of EPIC_Gene_3_200.rds; it has 2152
#   fewer clusters (1 less MB), but this could be due to changes in the Illumina
#   data set since 2018 (we are re-creating this data set in 2021, the original
#   was created three years ago).
# Because I can't re-create the old version exactly, I will include them both.
saveRDS(allClusters_min3_200bp_ls, "data/EPIC_10b4_Gene_3_200.rds")
