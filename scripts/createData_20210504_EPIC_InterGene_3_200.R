# Create EPIC_InterGene_3_200.rds Data Sets
# Gabriel Odom and Lissette Gomez
# 2021-05-02
# UPDATED: 2022-04-20

library(parallel)
library(coMethDMR)
library(tidyverse)



######  Map Genes to CpGs  ####################################################
cpgInterGene_df <- 
  IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other %>% 
	# IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Other %>% 
	subset(select = "UCSC_RefGene_Name") %>% 
	as.data.frame() %>% 
	rownames_to_column(var = "CpG") %>% 
	as_tibble() %>% # 865,859
	filter(UCSC_RefGene_Name == "") %>% # 249,261
	filter(str_detect(CpG, "cg")) # 247,484

locations_df <- 
  IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations %>% 
	# IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations %>% 
	as.data.frame() %>% 
	rownames_to_column(var = "CpG") %>% 
	as_tibble() 

cpgLocations_df <- 
  left_join(
    cpgInterGene_df,
    locations_df
  ) %>% 
	arrange(CpG)

allChrRegions_ls <- 
	cpgLocations_df %>% 
	split(f = cpgLocations_df$chr) %>% 
	map(pull, "CpG")



######  Find Regions with Clusters of CpGs  ###################################

system.time(
  regions_min3_200bp_ls <- map(
    allChrRegions_ls,
    coMethDMR::CloseBySingleRegion,
    arrayType = "EPIC",
    maxGap = 200,
    minCpGs = 3
  )
)
# 14.5 sec


# There may be errors here. For example, "cg01420942" is included in the source
#   data (IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Other), but not
#   in the sesameDataGet("EPIC.hg19.manifest") database. We saw similar behavior
#   for the Gene Body regions (we lost about 2k genes because of this same
#   "subscript contains invalid names" error). The OrderCpGsByLocation()
#   function calls this, so that's where the error handling should be.
# UPDATE: see issue #4: https://github.com/TransBioInfoLab/coMethDMR/issues/4


# All clusters within chromosomes (considered independently)
allClusters_min3_200bp_ls <- 
	regions_min3_200bp_ls %>% 
	unlist(recursive = FALSE) %>% 
	unique()



######  Name Regions  #########################################################
# NOTE: this parallel computing works in UNIX environments. For Windows, please
#   use parLapply or similar variants: 
# http://dept.stat.lsa.umich.edu/~jerrick/courses/stat701/notes/parallel.html


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
# Parallel: 15.1388 min for 10,908 clusters


# Windows
system.time(
  allClustersDF_ls <- parLapply(
    cl = makeCluster(4L),
    X = allClusters_min3_200bp_ls,
    fun = coMethDMR::OrderCpGsByLocation,
    arrayType = "EPIC",
    output = "dataframe"
  )
)
# Parallel: 8.226667 min for 10,908 clusters (this is after we updated the
#   OrderCpGsByLocation() function, so it's not apples to oranges to compare Mac
#   to Windows)
# 42.21 for 100/2 cores; 39.28 for 100/3; 38.32 for 100/4; 41.86 for 100/5; 


names(allClusters_min3_200bp_ls) <- map(allClustersDF_ls, NameRegion)
# OLD COMMENT:
# We believe this is the updated version of EPIC_InterGene_3_200.rds; it has 189
#   additional clusters (0.1 additional MB), but this could be due to changes in
#   the Illumina data set since 2018 (we are re-creating this data set in 2021,
#   the original was created three years ago).
# Because I can't re-create the old version exactly, I will include them both.
# END OLD COMMENT
# saveRDS(allClusters_min3_200bp_ls, "data/EPIC_InterGene_3_200_new.rds")

saveRDS(allClusters_min3_200bp_ls, "data/EPIC_10b4_InterGene_3_200.rds")
