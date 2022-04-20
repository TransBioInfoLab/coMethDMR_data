# coMethDMR_data
This repository contains utility data sets for additional examples in the vignettes of <https://github.com/TransBioInfoLab/coMethDMR>.

- `scripts/` contains the R scripts necessary to create these data sets
- `data/` contains the compressed data sets (in .RDS form)
    + We recently recalculated the EPIC regions using the 10b4 dataset. Please use `EPIC_10b4_Gene_3_200.rds`. 

## Necessary Packages and Data Sets
You will need the following packages:

- the Tidyverse
- our package coMethDMR
- Illumina EPIC annotation `IlluminaHumanMethylationEPICanno.ilm10b2.hg19`
- the parallel package
- IMA's 450k region-level annotation library (specifically the `ISLANDSInd` data set) from <https://rforge.net/IMA/#sec-4>
