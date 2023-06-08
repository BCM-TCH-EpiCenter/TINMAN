require(tidyverse); require(readxl)

# Load merged phase one and two data --------------------------------------

setwd('//smb-main.ad.bcm.edu/genepi2/TINMAN/Datasets/BAYL-03-22MD_merged/')

#' Load Metabolon data.
feces <- read_xlsx('BAYL-13-22MD CO DATA TABLES.XLSX', sheet = 'Batch-norm Data Common') 
names(feces) <- c(names(feces)[1], paste0('C',names(feces)[2:ncol(feces)]))

metab.metadata <- read_xlsx('BAYL-13-22MD CO DATA TABLES.XLSX', sheet = 'Chemical Annotation Common') %>% 
  mutate(CHEM_ID = paste0('C', CHEM_ID))

#' Load and combine clinical metadata from each phase.
clin.one.meta <- readRDS('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_sample_metadata.rds') %>% 
  select(PARENT_SAMPLE_NAME, CLIENT_IDENTIFIER, starts_with('GROUP'))

clin.two.meta <- readRDS('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_phase_2_sample_metadata.rds') %>% 
  select(PARENT_SAMPLE_NAME, CLIENT_IDENTIFIER, starts_with('GROUP'))

clinical.metadata <- bind_rows(clin.one.meta, clin.two.meta)

#' Per Josie's email on 9/26/2022, for patients 53 and 61, the second sample (D, E, or F) is the baseline sample and
#' the first sample (A, B, or C) is the neutropenic sample. Metabolon group assignments do not reflect this.
clinical.metadata <- clinical.metadata %>% 
  mutate(GROUP_ID = ifelse(CLIENT_IDENTIFIER %in% c('053F1','061F2'), 'NEUT_PRE',
                    ifelse(CLIENT_IDENTIFIER %in% c('053B2','061C2'), 'NEUT_POST', GROUP_ID)))

#' Standardize group names across the two files.
clinical.metadata <- clinical.metadata %>% 
  mutate(GROUP_ID = toupper(str_replace(GROUP_ID, 'Con', 'CTRL')),
         GROUP_NAME = case_when(GROUP_NAME == 'neutropenic' ~ 'Neutropenic', 
                                .default = GROUP_NAME),
         TIMEPOINT = ifelse(str_detect(GROUP_ID, 'POST'), 'POST', 'PRE'))

saveRDS(clinical.metadata, '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_merged_clinical_data.rds')

saveRDS(feces, '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_merged_feces_raw_data.rds')

saveRDS(metab.metadata, '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_merged_feces_metadata.rds')

rm(list = ls()); gc()
