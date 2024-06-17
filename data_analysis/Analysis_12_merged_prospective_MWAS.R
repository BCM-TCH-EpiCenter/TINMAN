require(tidyverse); require(ggrepel)

# Load and prepare data ---------------------------------------------------

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

#' Metabolite metadata
metab.metadata <- readRDS('TINMAN_merged_feces_metadata.rds') 

metabs <- metab.metadata %>% 
  filter(SUPER_PATHWAY != 'Xenobiotics' | (SUPER_PATHWAY == 'Xenobiotics' & SUB_PATHWAY == 'Bacterial/Fungal')) %>% 
  pull(CHEM_ID)

#' Concentration table.
#' Choice of scaled or unscaled data does not impact the algorithm because it performs its own centering and scaling.
tinman <- readRDS('TINMAN_prospective_merged_metab_imputed_unscaled_20240611.rds') 

sort.order <- tinman$PARENT_SAMPLE_NAME

#' Load and append clinical metadata.
clinical.metadata <- readRDS('TINMAN_merged_clinical_data_prospective_samples.rds') %>% 
  arrange(sort.order)

tinman <- tinman %>% 
  left_join(select(clinical.metadata, PARENT_SAMPLE_NAME:GROUP_NAME, TIMEPOINT), by = 'PARENT_SAMPLE_NAME') %>% 
  select(PARENT_SAMPLE_NAME, CLIENT_IDENTIFIER, GROUP_ID, GROUP_NAME, TIMEPOINT, SUBJECT, all_of(metabs))

#' For the 100 metabolites that have the highest loadings on components 1 and 2, collate and arrange pre-/post- abundances.
tinman.wide <- tinman %>% 
  select(CLIENT_IDENTIFIER, GROUP_NAME, TIMEPOINT, all_of(metabs)) %>% 
  mutate(subject = str_trim(str_remove(CLIENT_IDENTIFIER, '[:alpha:]{1}[:digit:]{1}$'), 'right')) %>% # Create a variable that identifies the subjects but not the timepoint.
  select(-CLIENT_IDENTIFIER) %>% 
  select(subject, GROUP_NAME, TIMEPOINT, starts_with('C')) %>% 
  pivot_wider(id_cols = c(subject, GROUP_NAME),
              names_from = TIMEPOINT,
              values_from = starts_with('C'))

#' Compute difference in each metab between timepoints for each subject.
#' Define pre- and post- column pairs.
i <- 3 # First column with "post" values
j <- 4 # First column with "pre" values
max.col <- ncol(tinman.wide)

#' Until we reach the end of these sets, calculate post - pre fold change.
while (j <= max.col){
  
  #' CHEM_ID for the index metabolite.
  my.var <- names(tinman.wide)[i] %>% str_remove('_[:alpha:]{3,}')
  
  #' Subtract pre-exposure values from post-exposure values to get the delta value. 
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

# Perform MWAS ------------------------------------------------------------

#' A vector of metabolites to iterate over.
metabs <- tinman.wide %>% 
  select(starts_with('C')) %>% 
  names() %>% 
  str_remove('_PRE') %>% 
  subset(!duplicated(.))

#' Empty data frame to hold results.
results <- data.frame()

for (i in seq_along(metabs)){
  
  #' Select the variables for neutropenia status, index metabolite baseline abundance, and index metabolite delta value.
  tmp <- select(tinman.wide, GROUP_NAME, metabs[i], paste0(metabs[i], '_PRE'))
  names(tmp) <- c('group', 'delta', 'baseline')
  
  #' For those metabolites with non-zero variance, estimate odds of neutropenia per SD change across the study.
  if (var(tmp$baseline, na.rm = T) > 0) {
    
    baseline.model <- glm(group ~ baseline, data = tmp, family = 'binomial') %>% summary()
    
    delta.model <- glm(group ~ baseline + delta, data = tmp, family = 'binomial') %>% summary()
    
    new.result <- data.frame(metab = metabs[i],
                             
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
  left_join(select(metab.metadata, CHEM_ID, SUPER_PATHWAY, SUB_PATHWAY, CHEMICAL_NAME, HMDB, KEGG, PUBCHEM), by = c('metab' = 'CHEM_ID')) %>% 
  select(metab, CHEMICAL_NAME, SUPER_PATHWAY, SUB_PATHWAY, HMDB, PUBCHEM, KEGG, starts_with('delta'), p.value.delta, starts_with('baseline'), p.value.baseline) %>% 
  mutate(delta.estimate = paste0(delta.or, ' (', delta.ci.lower, '-', delta.ci.upper, ')'),
         baseline.estimate = paste0(baseline.or, ' (', baseline.ci.lower, '-', baseline.ci.upper, ')'))

# Annotate results --------------------------------------------------------

#' Loading scores for metabolites in PLS-DA components 1 and 2.
loadings <- read_csv('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/prospective_samples_plsa_da_component_loadings_20240612.csv')

results %>% 
  filter(p.value.delta < 0.05 | p.value.baseline < 0.05) %>% 
  select(CHEMICAL_NAME, delta.or, baseline.or, p.value.delta, p.value.baseline) %>% 
  arrange(p.value.delta)

#' Re-format for an export to Excel.
results <- results %>% 
  left_join(select(loadings, CHEMICAL_NAME, component, rank), by = 'CHEMICAL_NAME') %>% 
  #dplyr::select(CHEMICAL_NAME, HMDB, PUBCHEM, KEGG, SUPER_PATHWAY, SUB_PATHWAY, component, loading, rank, delta.estimate, p.value.delta, baseline.estimate, p.value.baseline) %>% 
  rename_all(toupper) %>% 
  rename_all(~str_replace_all(.x, '\\Q.\\E', '_')) %>% 
  rename(`BASELINE_OR_(95% CI)` = BASELINE_ESTIMATE,
         `DELTA_OR_(95% CI)` = DELTA_ESTIMATE)

results %>% 
  filter(P_VALUE_BASELINE < 0.05 | P_VALUE_DELTA < 0.05) %>% 
  filter(METAB %in% loadings$CHEM_ID)

write_csv(results,
          paste0('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/prospective_mwas_results_unscaled_', format(Sys.Date(), '%Y%m%d'), '.csv'))

# Compare to prior results from whole cohort ------------------------------

results %>% filter(P_VALUE_DELTA < 0.05 | P_VALUE_BASELINE < 0.05)

#' Results from whole cohort.
full.results <- read_csv('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/mwas_results_unscaled_20240228.csv')

full.results %>% filter(P_VALUE_DELTA < 0.05 | P_VALUE_BASELINE < 0.05)

comparison <- full.results %>% 
  filter(P_VALUE_DELTA < 0.05 | P_VALUE_BASELINE < 0.05) %>% 
  select(METAB, CHEMICAL_NAME, `BASELINE_OR_(95% CI)`, `DELTA_OR_(95% CI)`) %>%
  rename(Baseline_value_OR_and_CI_prospective = `BASELINE_OR_(95% CI)`,
         Delta_value_OR_and_CI_prospective = `DELTA_OR_(95% CI)`) %>% 
  left_join(select(results, METAB, `BASELINE_OR_(95% CI)`,`DELTA_OR_(95% CI)`), by = 'METAB') %>%
  rename(Baseline_value_OR_and_CI_all_samples = `BASELINE_OR_(95% CI)`,
         Delta_value_OR_and_CI_all_samples = `DELTA_OR_(95% CI)`) %>% 
  select(METAB, CHEMICAL_NAME, starts_with('Baseline'), starts_with('Delta'))

write_csv(comparison,
          '//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/TINMAN_comparison_of_log_reg_models_full_vs_prospective_data_20240612.csv')
