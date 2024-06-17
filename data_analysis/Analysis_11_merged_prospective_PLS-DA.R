# Prepare data ------------------------------------------------------------

require(tidyverse)

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

#' Metabolite metadata
metab.metadata <- readRDS('TINMAN_merged_feces_metadata.rds') 

metabs <- metab.metadata %>% 
  filter(SUPER_PATHWAY != 'Xenobiotics' | (SUPER_PATHWAY == 'Xenobiotics' & SUB_PATHWAY == 'Bacterial/Fungal')) %>% 
  pull(CHEM_ID)

#' Concentration table.
#' Choice of scaled or unscaled data does not impact the algorithm because it performs its own centering and scaling.
tinman <- readRDS('TINMAN_prospective_merged_metab_imputed_unscaled_20240611.rds') %>% 
  dplyr::select(PARENT_SAMPLE_NAME, all_of(metabs))

sort.order <- tinman$PARENT_SAMPLE_NAME

#' Load and append clinical metadata.
clinical.metadata <- readRDS('TINMAN_merged_clinical_data_prospective_samples.rds') %>% 
  arrange(sort.order)

tinman <- tinman %>% 
  left_join(dplyr::select(clinical.metadata, PARENT_SAMPLE_NAME:GROUP_NAME), by = 'PARENT_SAMPLE_NAME') %>% 
  dplyr::select(PARENT_SAMPLE_NAME, GROUP_ID, GROUP_NAME, SUBJECT, all_of(metabs))

# Matrix of metabolite abundances (x) and class memberships (y) needs for sparse PLS-DA.
x <- as.matrix(tinman[ , 5:ncol(tinman)])
y <- as.factor(tinman$GROUP_ID)

# 'Stock' PLS-DA ---------------------------------------------------------

#' Run SPLS-DA analysis.
#' Note that if you are running multilevel PLS-DA, you cannot use R version 4.3.1. 
#' Version 4.2.1 seems to work.
MyResult.splsda.final <- splsda(x, y, ncomp = 2, keepX = c(50,50))

#' Save PLS-DA plot
svg('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/PLS-DA/prospective.samples.merged.plsa.da.plot.20240612.svg',
    width = 10, height = 10)

plotIndiv(MyResult.splsda.final, 
          ind.names = FALSE, 
          legend=TRUE,
          ellipse = TRUE, 
          title="sPLS-DA - final result")

dev.off()

# Loading plots -----------------------------------------------------------

#' Recover and plot loading vectors for components 1 and 2
svg(paste0('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/PLS-DA/prospective.samples.merged.plsa.da.comp1.loadings', format(Sys.Date(), '%Y%m%d'), '.svg'),
    width = 10, height = 10)

plotLoadings(MyResult.splsda.final, 
             contrib = 'max', 
             method = 'mean', 
             comp = 1,
             name.var = metabs)

dev.off()

svg(paste0('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/PLS-DA/prospective.samples.merged.plsa.da.comp2.loadings', format(Sys.Date(), '%Y%m%d'), '.svg'),
    width = 10, height = 10)

plotLoadings(MyResult.splsda.final, 
             contrib = 'max', 
             method = 'mean', 
             comp = 2,
             name.var = metabs)

dev.off()

# Extract and save sPLDS-DA loadings --------------------------------------

comp1.loadings <- selectVar(MyResult.splsda.final, comp=1)$value %>% 
  mutate(CHEM_ID = rownames(.),
         component = 1) %>% 
  arrange(desc(abs(value.var))) %>% 
  mutate(rank = 1:nrow(.))

comp2.loadings <- selectVar(MyResult.splsda.final, comp=2)$value %>% 
  mutate(CHEM_ID = rownames(.),
         component = 2) %>% 
  arrange(desc(abs(value.var))) %>% 
  mutate(rank = 1:nrow(.))

loadings <- bind_rows(comp1.loadings, comp2.loadings) %>% 
  left_join(dplyr::select(metab.metadata, CHEM_ID, CHEMICAL_NAME)) %>% 
  rename(loading = value.var) %>% 
  dplyr::select(CHEMICAL_NAME, CHEM_ID, component, loading, rank) %>% 
  as_tibble()

write_csv(loadings, 
          paste0('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/prospective_samples_plsa_da_component_loadings_', format(Sys.Date(), '%Y%m%d'), '.csv'))

# Compare PLS-DA metabolites from total and prospective data --------------

full.sample.loadings <- read_csv('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/plsa_da_component_loadings_20240111.csv')

#' 39 of 100 metabolites in both sets.
comp1pro <- loadings %>% filter(component == 1) %>% pull(CHEM_ID)
comp1all <- full.sample.loadings %>% filter(component == 1) %>% pull(CHEM_ID)

#' 38/50 component 1 metabolites are the same comparing the full and prospective analyses, whereas only 1 of 50 component 2 metabolites are the same.
intersect(comp1all, comp1pro) %>% length()

# 'Tuned' PLS-DA; generally the same findings -----------------------------

#' Install mixOmics, if neded.
#' BiocManager::install('mixOmics')

require(mixOmics)

tinman.splsda <- splsda(x, y, ncomp = 10)

# Run tuning procedure in order to tune the number of components to use
perf.splsda <- perf(tinman.splsda, validation = "Mfold", 
                    folds = 5, nrepeat = 50, # use repeated cross-validation
                    progressBar = FALSE, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10,  seq(20, 200, 10))

# undergo the tuning process to determine the optimal number of variables
tune.splsda <- tune.splsda(x, y, ncomp = 4, # calculate for first 4 components; most error measurements are lower/lowest at this value of ncomp.
                           validation = 'Mfold',
                           folds = 5, nrepeat = 50, # use repeated cross-validation
                           dist = 'max.dist', # use max.dist measure
                           measure = "BER", # use balanced error rate of dist measure
                           test.keepX = list.keepX,
                           cpus = 2) # allow for parallelization to decrease runtime

plot(tune.splsda, col = color.jet(4)) # plot output of variable number tuning

tune.splsda$choice.ncomp$ncomp # what is the optimal value of components according to tune.splsda()
tune.splsda$choice.keepX # what are the optimal values of variables according to tune.splsda()

optimal.ncomp <- tune.splsda$choice.ncomp$ncomp
optimal.keepX <- tune.splsda$choice.keepX[1:optimal.ncomp]

final.splsda <- splsda(x, y, 
                       ncomp = optimal.ncomp, 
                       keepX = optimal.keepX)

plotIndiv(final.splsda, 
          ind.names = FALSE, 
          legend=TRUE,
          ellipse = TRUE, 
          title="sPLS-DA - final result")


#' Export PubChem IDs for metabolites in the background set, as well as metabolites included in components 1 and 2.
#' These will be used to generate name mapping files, which will subsequently be used for overrepresentation analysis (ORA).
background <- metab.metadata %>% 
  filter((CHEM_ID %in% metabs), !is.na(PUBCHEM)) %>% 
  pull(PUBCHEM)

write.table(background,
            paste0('TINMAN_ORA_background_set_', format(Sys.Date(), '%Y%m%d'), '.txt'),
            row.names = F, col.names = F, quote = F)

comp1 <- loadings %>% 
  filter(component == 1) %>% 
  arrange(abs(loading)) %>% 
  left_join(dplyr::select(metab.metadata, CHEM_ID, PUBCHEM)) %>% 
  pull(PUBCHEM)

write.table(comp1,
            paste0('TINMAN_ORA_sPLDS_DA_comp1_', format(Sys.Date(), '%Y%m%d'), '.txt'),
            row.names = F, col.names = F, quote = F)

comp2 <- loadings %>% 
  filter(component == 2) %>% 
  arrange(abs(loading)) %>% 
  left_join(dplyr::select(metab.metadata, CHEM_ID, PUBCHEM)) %>% 
  pull(PUBCHEM)

write.table(comp2,
            paste0('TINMAN_ORA_sPLDS_DA_comp2_', format(Sys.Date(), '%Y%m%d'), '.txt'),
            row.names = F, col.names = F, quote = F)

rm(list = ls()); gc()

#' Load in the name mapping files produced by the MetaboAnalyst compound conversion tool.
setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

background <- read.csv('TINMAN_ORA_PLS_DA_background_name_map.csv') %>% 
  filter(!is.na(Match)) %>% 
  pull(Match)

comp1 <- read.csv('TINMAN_ORA_PLS_DA_comp1_name_map.csv') %>% 
  filter(!is.na(Match)) %>% 
  pull(Match)

comp2 <- read.csv('TINMAN_ORA_PLS_DA_comp2_name_map.csv') %>% 
  filter(!is.na(Match)) %>% 
  pull(Match)

write.table(background,
            paste0('TINMAN_ORA_background_set_', format(Sys.Date(), '%Y%m%d'), '.txt'),
            row.names = F, col.names = F, quote = F)

write.table(comp1,
            paste0('TINMAN_ORA_sPLDS_DA_comp1_', format(Sys.Date(), '%Y%m%d'), '.txt'),
            row.names = F, col.names = F, quote = F)

write.table(comp2,
            paste0('TINMAN_ORA_sPLDS_DA_comp2_', format(Sys.Date(), '%Y%m%d'), '.txt'),
            row.names = F, col.names = F, quote = F)

# Prepare data for ORA ----------------------------------------------------

#' Export PubChem IDs for metabolites in the background set, as well as metabolites included in components 1 and 2.
#' These will be used to generate name mapping files, which will subsequently be used for overrepresentation analysis (ORA).
background <- metab.metadata %>% 
  filter((CHEM_ID %in% metabs), !is.na(PUBCHEM)) %>% 
  pull(PUBCHEM)

write.table(background,
            paste0('TINMAN_prospective_ORA_background_set_', format(Sys.Date(), '%Y%m%d'), '.txt'),
            row.names = F, col.names = F, quote = F)

comp1 <- loadings %>% 
  filter(component == 1) %>% 
  arrange(abs(loading)) %>% 
  left_join(dplyr::select(metab.metadata, CHEM_ID, PUBCHEM)) %>% 
  pull(PUBCHEM) %>% 
  format(scientific = F)

write.table(comp1,
            paste0('TINMAN_prospective_ORA_sPLDS_DA_comp1_', format(Sys.Date(), '%Y%m%d'), '.txt'),
            row.names = F, col.names = F, quote = F)

comp2 <- loadings %>% 
  filter(component == 2) %>% 
  arrange(abs(loading)) %>% 
  left_join(dplyr::select(metab.metadata, CHEM_ID, PUBCHEM)) %>% 
  pull(PUBCHEM) %>% 
  format(scientific = F)

write.table(comp2,
            paste0('TINMAN_prospective_ORA_sPLDS_DA_comp2_', format(Sys.Date(), '%Y%m%d'), '.txt'),
            row.names = F, col.names = F, quote = F)

rm(list = ls()); gc()