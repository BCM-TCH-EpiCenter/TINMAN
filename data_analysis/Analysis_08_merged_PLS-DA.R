# Prepare data ------------------------------------------------------------

require(tidyverse); require(xlsx)

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

#' Metabolite metadata
metab.metadata <- readRDS('TINMAN_merged_feces_metadata.rds') 

metabs <- metab.metadata %>% 
  filter(SUPER_PATHWAY != 'Xenobiotics' | (SUPER_PATHWAY == 'Xenobiotics' & SUB_PATHWAY == 'Bacterial/Fungal')) %>% 
  pull(CHEM_ID)

#' Concentration table.
#' Choice of scaled or unscaled data does not impact the algorithm because it performs its own centering and scaling.
tinman <- readRDS('TINMAN_merged_metab_imputed_unscaled_20230606.rds') %>% 
  dplyr::select(PARENT_SAMPLE_NAME, all_of(metabs))

sort.order <- tinman$PARENT_SAMPLE_NAME

#' Load and append clinical metadata.
clinical.metadata <- readRDS('TINMAN_merged_clinical_data.rds') %>% 
  arrange(sort.order)

tinman <- tinman %>% 
  left_join(dplyr::select(clinical.metadata, PARENT_SAMPLE_NAME:GROUP_NAME), by = 'PARENT_SAMPLE_NAME') %>% 
  dplyr::select(PARENT_SAMPLE_NAME, GROUP_ID, GROUP_NAME, SUBJECT, all_of(metabs))

# Matrix of metabolite abundances (x) and class memberships (y) needs for sparse PLS-DA.
x <- as.matrix(tinman[ , 5:ncol(tinman)])
y <- as.factor(tinman$GROUP_ID)

# Tuning keepX, the number of features that will be used for each component.
#list.keepX <- c(5:10,  seq(20, 50, 10))

#set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
#tune.splsda.srbct <- tune.splsda(x, y, ncomp = 2, # we suggest to push ncomp a bit more, e.g. 4
#                                 validation = 'Mfold', folds = 3, 
#                                 dist = 'max.dist', progressBar = T,
#                                 measure = "BER", test.keepX = list.keepX,
#                                 nrepeat = 50)   # we suggest nrepeat = 50

#error <- tune.splsda.srbct$error.rate
#ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
#ncomp

#select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
#select.keepX

# Perform PLS-DA ----------------------------------------------------------

#' Install mixOmics, if neded.
#' BiocManager::install('mixOmics')

require(mixOmics)

#' Run SPLS-DA analysis.
#' Note that if you are running multilevel PLS-DA, you cannot use R version 4.3.1. 
#' Version 4.2.1 seems to work.
MyResult.splsda.final <- splsda(x, y, ncomp = 2, keepX = c(50,50))

#' Save PLS-DA plot
svg('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/PLS-DA/merged.plsa.da.plot.20240111.svg',
    width = 10, height = 10)

plotIndiv(MyResult.splsda.final, 
          ind.names = FALSE, 
          legend=TRUE,
          ellipse = TRUE, 
          title="sPLS-DA - final result")

dev.off()

y <- as.factor(clinical.metadata$TIMEPOINT)

design <- data.frame(sample = clinical.metadata$SUBJECT)

splsda.ml <- spls(x, y, ncomp = 2, keepX = c(50,50))

plotIndiv(splsda.ml, 
          ind.names = FALSE, 
          legend=TRUE,
          ellipse = TRUE, 
          title="sPLS-DA - final result")

# Loading plots -----------------------------------------------------------

#' Recover and plot loading vectors for components 1 and 2
svg(paste0('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/PLS-DA/merged.plsa.da.comp1.loadings', format(Sys.Date(), '%Y%m%d'), '.svg'),
    width = 10, height = 10)

plotLoadings(MyResult.splsda.final, 
             contrib = 'max', 
             method = 'mean', 
             comp = 1,
             name.var = metabs)

dev.off()

svg(paste0('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/PLS-DA/merged.plsa.da.comp2.loadings', format(Sys.Date(), '%Y%m%d'), '.svg'),
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
          paste0('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/plsa_da_component_loadings_', format(Sys.Date(), '%Y%m%d'), '.csv'))

# Correlation circle plot -------------------------------------------------

#' Metabolite metadata
metabs <- readRDS('TINMAN_merged_feces_metadata.rds')

var.names <- metabs$CHEMICAL_NAME

svg('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/PLS-DA/merged.plsa.da.correlation.circle.plot.svg',
    width = 10, height = 10)

plotVar(MyResult.splsda.final, 
        comp = c(1,2),
        cex = 3,
        cutoff = 0.6,
        var.names = list(var.names),
        X.label = 'Correlation with Component 1',
        Y.label = 'Correlation with Component 2') 

dev.off()

# CIM plot ----------------------------------------------------------------

svg('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/PLS-DA/merged.plsa.da.hierarchical.heatmap.one.component.svg',
    width = 10, height = 10)

cim(MyResult.splsda.final,
    row.sideColors = color.mixo(MyResult.splsda.final$Y),
    row.names = F, 
    transpose = T,
    col.names = var.names)

dev.off()








#' Append metadata to loadings
comp1.loadings <- left_join(comp1.loadings, metabs, by = 'CHEM_ID')

comp2.loadings <- left_join(comp2.loadings, metabs, by = 'CHEM_ID')

write.xlsx(comp1.loadings,
           '//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/TINMAN_merged_PLS-DA_loadings_20230614.xlsx',
           sheetName = 'component1',
           row.names = F)

write.xlsx(comp2.loadings,
           '//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/TINMAN_merged_PLS-DA_loadings_20230614.xlsx',
           sheetName = 'component2',
           row.names = F,
           append = T)

auc.plsda <- auroc(MyResult.splsda.final)

#' Evaluate overlap between loadings and univariate results
uni <- read.xlsx('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/TINMAN_merged_metabolomics_analysis_20230614.xlsx',
                 sheetIndex = 1,
                 colIndex = 2:13)

uni <- uni %>% filter(t.test.p < 0.05 | kruskal.test.p < 0.05)

loadings <- c(comp1.loadings$CHEM_ID, comp2.loadings$CHEM_ID)

intersect(uni$CHEM_ID, loadings)


# Calculate delta values --------------------------------------------------

#' For the 100 metabolites that have the highest loadings on components 1 and 2, collate and arrange pre-/post- abundances.
tinman.wide <- tinman %>% 
  left_join(dplyr::select(clinical.metadata, PARENT_SAMPLE_NAME, CLIENT_IDENTIFIER, TIMEPOINT), by = 'PARENT_SAMPLE_NAME') %>% 
  dplyr::select(CLIENT_IDENTIFIER, GROUP_NAME, TIMEPOINT, all_of(loadings$CHEM_ID)) %>% 
  mutate(subject = str_trim(str_remove(CLIENT_IDENTIFIER, '[:alpha:]{1}[:digit:]{1}$'), 'right')) %>% # Create a variable that identifies the subjects.
  dplyr::select(-CLIENT_IDENTIFIER) %>% 
  dplyr::select(subject, GROUP_NAME, TIMEPOINT, starts_with('C')) %>% 
  pivot_wider(id_cols = c(subject, GROUP_NAME),
              names_from = TIMEPOINT,
              values_from = starts_with('C'))

#' Compute difference in each metab between timepoints for each subject.
#' Define pre- and post- column pairs.
i <- 3
j <- 4
max.col <- ncol(tinman.wide)

#' Until we reach the end of these sets, calculate post - pre fold change.
while (j <= max.col){
  
  my.var <- names(tinman.wide)[i] %>% str_remove('_[:alpha:]{3,}')
  
  #' Subtract pre-induction values from post-induction values to get the delta value. 
  my.values <- (tinman.wide[,i] - tinman.wide[,j]) %>% as.vector() %>% unlist()
  
  #' Assign these to a new column.
  tinman.wide[[my.var]] <- my.values
  
  #tinman.wide <- cbind(tinman.wide, my.values)
  #names(tinman.wide) <- c(names(tinman.wide)[1:(ncol(tinman.wide)-1)], my.var)
  
  i <- i + 2
  j <- j + 2
  
}

#' Remove columns holding post-treatment abundances. Retain pre-treatment abundances since these will be covariates in models.
#' Convert GROUP_NAME to factor.
tinman.wide <- tinman.wide %>% 
  dplyr::select(!contains('POST')) %>% 
  mutate(GROUP_NAME = factor(GROUP_NAME))

# Associate metabolites with neutropenia ----------------------------------

#' Empty data frame to hold results.
results <- data.frame()

for (i in seq_along(loadings$CHEM_ID)){

  #' Select the variables for neutropenia status, index metabolite baseline abundance, and index metabolite delta value.
  tmp <- dplyr::select(tinman.wide, GROUP_NAME, all_of(loadings$CHEM_ID)[i], all_of(paste0(loadings$CHEM_ID[i], '_PRE')))
  names(tmp) <- c('group', 'change', 'baseline')
  
  #' For those metabolites with non-zero variance, estimate odds of neutropenia per SD change across the study.
  if (var(tmp$baseline, na.rm = T) > 0) {
    
    baseline.model <- glm(group ~ baseline, data = tmp, family = 'binomial') %>% summary()
    
    delta.model <- glm(group ~ baseline + change, data = tmp, family = 'binomial') %>% summary()
    
    new.result <- data.frame(metab = loadings$CHEM_ID[i],
                             
                             coef.delta = delta.model$coefficients[3,1],
                             se.delta = delta.model$coefficients[3,2],
                             p.value.delta = delta.model$coefficients[3,4],
                             
                             coef.baseline = baseline.model$coefficients[2,1],
                             se.baseline = baseline.model$coefficients[2,2],
                             p.value.baseline = baseline.model$coefficients[2,4])
    
    results <- rbind(results, new.result)
    
  }
  
  else {
    
    next
    
  }
  
}

#' Compute ORs and CIs.
#' Note that the OR for the delta value is calculated from the -1 * the coefficient for the delta value, because we are interested in reporting the OR per one SD decrease (rather than increase) in metabolite abundance.
results <- results %>% 
  mutate(delta.or = exp(coef.delta),
         delta.ci.lower = exp(coef.delta - (1.96*se.delta)),
         delta.ci.upper = exp(coef.delta + (1.96*se.delta)),
         
         baseline.or = exp(coef.baseline),
         baseline.ci.lower = exp(coef.baseline - (1.96*se.baseline)),
         baseline.ci.upper = exp(coef.baseline + (1.96*se.baseline)),
         across(delta.or:baseline.ci.upper, ~round(.x, 2)))

#' Append metabolite metadata.
#' Append information on sPLS-DA loadings.
results <- results %>% 
  as_tibble() %>% 
  left_join(dplyr::select(metab.metadata, CHEM_ID, SUPER_PATHWAY, SUB_PATHWAY, CHEMICAL_NAME, HMDB, KEGG, PUBCHEM), by = c('metab' = 'CHEM_ID')) %>% 
  dplyr::select(metab, CHEMICAL_NAME, SUPER_PATHWAY, SUB_PATHWAY, HMDB, PUBCHEM, KEGG, starts_with('delta'), p.value.delta, starts_with('baseline'), p.value.baseline) %>% 
  left_join(dplyr::select(loadings, CHEM_ID, component, loading, rank), by = c('metab' = 'CHEM_ID')) %>% 
  mutate(delta.estimate = paste0(delta.or, ' (', delta.ci.lower, '-', delta.ci.upper, ')'),
         baseline.estimate = paste0(baseline.or, ' (', baseline.ci.lower, '-', baseline.ci.upper, ')'))

results %>% 
  filter(p.value.delta < 0.1 | p.value.baseline < 0.1) %>% 
  dplyr::select(CHEMICAL_NAME, component, rank, delta.or, baseline.or, p.value.delta, p.value.baseline) %>% 
  arrange(p.value.delta)

#' Re-format for an export to Excel.
final.results <- results %>% 
  dplyr::select(CHEMICAL_NAME, HMDB, PUBCHEM, KEGG, SUPER_PATHWAY, SUB_PATHWAY, component, loading, rank, delta.estimate, p.value.delta, baseline.estimate, p.value.baseline) %>% 
  rename_all(toupper) %>% 
  rename_all(~str_replace_all(.x, '\\Q.\\E', '_')) %>% 
  rename(`BASELINE_OR_(95% CI)` = BASELINE_ESTIMATE,
         `DELTA_OR_(95% CI)` = DELTA_ESTIMATE,
         COMPONENT_LOADING = LOADING,
         COMPONENT_RANK = RANK)

write_csv(final.results,
          paste0('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/odds_of_neutropenia_', format(Sys.Date(), '%Y%m%d'), '.csv'))

# Random forests ----------------------------------------------------------

require(caret)

# Cross-validation parameters for random forests models.
train_control <- trainControl(method = "repeatedcv", 
                              number = 5, 
                              repeats = 3,
                              predictionBounds = c(TRUE, TRUE))

#' Initliaze an empty data frame to hold performance statistics.
perform <- data.frame()

#' Initialize an empty list to hold variable importance statistics.
important.vars <- list()

rf.data <- tinman.wide %>% 
  dplyr::select(-subject)

#' Train models.
rf<- train(GROUP_NAME ~ ., 
           data = rf.data, 
           trControl = train_control,
           method = 'rf',
           tuneLength = 10,
           preProcess = c('center', 'scale'))

rf$finalModel

probs <- cbind(rf$finalModel$votes, rf.data$GROUP_NAME)

# Prepare data for ORA ----------------------------------------------------

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
