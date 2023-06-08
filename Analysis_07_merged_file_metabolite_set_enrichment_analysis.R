# Install/update MetaboAanalystR (if needed) ------------------------------

metanr_packages <- function(){
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR", "fgsea", "devtools", "crmn")
  list_installed <- installed.packages()
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  if(length(new_pkgs)!=0){if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}

metanr_packages()

#' Install MetaboAnalystR.
require(devtools)

devtools::install_github('xia-lab/MetaboAnalystR', build = T, build_vignettes = T, build_manual = T)

# Load TINMAN data --------------------------------------------------------

require(tidyverse)

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

#' Concentration table
tinman <- readRDS('TINMAN_merged_20230606.rds') %>% 
  filter(str_detect(GROUP_ID, 'POST'))

#' Metabolite names
metabs <- readRDS('TINMAN_merged_feces_metadata.rds') %>% 
  filter(!is.na(HMDB), SUPER_PATHWAY != 'Xenobiotics' | (SUPER_PATHWAY == 'Xenobiotics' & SUB_PATHWAY == 'Bacterial/Fungal'))

tinman.msea <- tinman %>% select(PARENT_SAMPLE_NAME, GROUP_NAME, all_of(metabs$CHEM_ID))

names(tinman.msea) <- c(names(tinman.msea)[1:2], metabs$HMDB)

#' Some metabolites are linked to more than one HMDB ID (e.g., glyceric acid and l-glyceric acid). Strip off everything after the first.
names(tinman.msea) <- names(tinman.msea) %>% 
  str_remove_all('(?<=,).{1,}') %>% # Remove everything after commas...
  str_remove_all(',') # ...then remove the commas

#' Some of these HMDB IDs don't seem to work with the MetaboAnalyst web interface. 
#' Metabolites may have secondary accession numbers that MetaboAnalyst recognizes.
#' Update those where this is the case. They were discovered through trial and error.
names(tinman.msea) <- case_when(names(tinman.msea) == 'HMDB0029155' ~ 'HMDB34367',
                                .default = names(tinman.msea))

# Create MSEA data table.
write_csv(tinman.msea, 'TINMAN_merged_MSEA_input_data_20230608.csv')

# Run MSEA ----------------------------------------------------------------

require(MetaboAnalystR)

# Create mSetObj
mset <- InitDataObjects('conc', 'msetqea', F)

# Read in data table
mset <- Read.TextData(mSetObj = mset, filePath = 'TINMAN_merged_MSEA_post_treatment_input_data_20230608.csv')

# Perform cross-referencing of compound HMDB IDs. Looks like 565 have matches.
mset <- CrossReferencing(mSetObj = mset, q.type = 'hmdb')

# Create mapping results table
mset <- CreateMappingResultTable(mset)

# Mandatory check of data
mset.check <- SanityCheckData(mset)

# Set the metabolome filter
#mset.filter <- SetMetabolomeFilter(mset, F)

# Specify pathway database and minimum number of compounds per set.
mset <- SetCurrentMsetLib(mset, 'smpdb_pathway', 2)

# Calculate global test and score.
#' TODO: troubleshoot this step. 
mset <- CalculateGlobalTestScore(mset)



tmp <- mset$dataSet
