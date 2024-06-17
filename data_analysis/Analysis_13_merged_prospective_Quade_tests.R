require(tidyverse); require(ggrepel)

# Load and prepare data ---------------------------------------------------

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/')

#' Metabolite metadata
metab.metadata <- readRDS('Datasets/TINMAN_merged_feces_metadata.rds') 

metabs <- metab.metadata %>% 
  filter(SUPER_PATHWAY != 'Xenobiotics' | (SUPER_PATHWAY == 'Xenobiotics' & SUB_PATHWAY == 'Bacterial/Fungal')) %>% 
  pull(CHEM_ID)

#' Concentration table.
tinman <- readRDS('Datasets/TINMAN_prospective_merged_metab_imputed_unscaled_20240611.rds') 

sort.order <- tinman$PARENT_SAMPLE_NAME

#' Load and append clinical metadata.
clinical.metadata <- readRDS('Datasets/TINMAN_merged_clinical_data_prospective_samples.rds') %>% 
  arrange(sort.order)

tinman <- tinman %>% 
  left_join(select(clinical.metadata, PARENT_SAMPLE_NAME:GROUP_NAME, TIMEPOINT), by = 'PARENT_SAMPLE_NAME') %>% 
  select(PARENT_SAMPLE_NAME, CLIENT_IDENTIFIER, GROUP_ID, GROUP_NAME, TIMEPOINT, SUBJECT, all_of(metabs))

# Perform Quade tests -----------------------------------------------------

#' Convert to wide. 
tinman.block <- tinman %>% 
  select(CLIENT_IDENTIFIER, GROUP_NAME, TIMEPOINT, all_of(metabs)) %>% 
  mutate(subject = str_trim(str_remove(CLIENT_IDENTIFIER, '[:alpha:]{1}[:digit:]{1}$'), 'right')) %>% # Create a variable that identifies the subjects but not the timepoint.
  select(-CLIENT_IDENTIFIER) %>% 
  select(subject, GROUP_NAME, TIMEPOINT, starts_with('C')) 

groups <- list(Neutropenic = tinman.block %>% 
                 filter(GROUP_NAME == 'Neutropenic'),
               Control = tinman.block %>% 
                 filter(GROUP_NAME == 'Control'))

#' Empty data frame to hold results from Quade tests.
quade.p.values <- data.frame()

for (j in seq_along(groups)){
  
  #' Iterate over metabolites, conducting Quade tests for pre-/post- differences with subjects as blocks.
  for (i in metabs){
    
    my.data <- groups[[j]] %>% 
      select(subject, TIMEPOINT, all_of(i))
    
    names(my.data) <- c('subject', 'time', 'value')
    
    my.test <- quade.test(value ~ time | subject, data = my.data)
    
    new.result <- data.frame(group = names(groups)[j],
                             metab = i,
                             comparison = 'post/pre',
                             test = 'Quade test',
                             p.value = my.test$p.value)
    
    quade.p.values <- rbind(quade.p.values, new.result)
    
  }
  
}

# Generate fold change data -----------------------------------------------

#' Empty data frame to hold results.
fold.changes <- data.frame()

for (i in metabs){
  
  #my.data <- select(neutropenic, TIMEPOINT, all_of(i))
  my.data <- select(tinman, GROUP_NAME, TIMEPOINT, all_of(i))
  names(my.data) <- c('group', 'time', 'value')
  
  means <- aggregate(value ~ group + time, mean, data = my.data, na.rm = T) %>% 
    pivot_wider(id_cols = group, names_from = time, values_from = value) %>% 
    mutate(log2fc = log2(POST/PRE),
           metab = i)
  
  fold.changes <- rbind(fold.changes, means)
  
}

fold.changes <- fold.changes %>% 
  select(group, metab, POST:log2fc) %>% 
  left_join(quade.p.values, by = c('group', 'metab'))

#' "Widen" the data so it is easier to compare p-values and fold changes between the two groups.
fold.changes.wide <- fold.changes %>% 
  select(group, metab, log2fc, p.value) %>% 
  pivot_wider(id_cols = metab,
              names_from = group,
              values_from = c(log2fc, p.value))

#' Create vectors of hits in control and neutropenic subjects and note which are unique to one group.
neutropenic.hits <- fold.changes.wide %>% filter(p.value_Neutropenic < 0.05)

neut.only <- fold.changes.wide %>% 
  mutate(concordant = case_when(sign(log2fc_Control) == sign(log2fc_Neutropenic) ~ 1, .default = 0)) %>% 
  filter(p.value_Neutropenic < 0.05 & p.value_Control >= 0.05 | 
           p.value_Neutropenic < 0.05 & p.value_Control < 0.05 & concordant == 0) %>% 
  pull(metab)

control.hits <- fold.changes.wide %>% filter(p.value_Control < 0.05)

control.only <- fold.changes.wide %>% 
  mutate(concordant = case_when(sign(log2fc_Control) == sign(log2fc_Neutropenic) ~ 1, .default = 0)) %>% 
  filter(p.value_Neutropenic >= 0.05 & p.value_Control < 0.05 | 
           p.value_Neutropenic < 0.05 & p.value_Control < 0.05 & concordant == 0) %>% 
  pull(metab)

metab.metadata %>% filter(CHEM_ID %in% neut.only) %>% pull(CHEMICAL_NAME)

# Annotate neutropenic results --------------------------------------------

#' The metabolites Josie would like annotated.
neut.hits <- fold.changes %>% 
  filter(group == 'Neutropenic', p.value < 0.05, log2fc < 0) %>% 
  pull(metab)

#' Loading scores for metabolites in PLS-DA components 1 and 2.
loadings <- read_csv('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/prospective_samples_plsa_da_component_loadings_20240612.csv')

#' Filter to those in our set of interest.
loadings <- loadings %>% 
  select(-CHEMICAL_NAME) %>% 
  filter(CHEM_ID %in% neut.hits)

# Generate volcano plot ---------------------------------------------------

#' These are Metabolon metabolite accession numbers for a handful of species that Josie indicated she wanted highlighted in an email on 6/21/23.
#' (This was later updated to exclude creatine).
#my.metabs <- c('C391','C100001104', 'C100006435', 'C100000098', 'C925', 'C926', 'C100003239', 'C100001416', 'C100021504')

#' Annotated fold changes and p-values for plotting.
plot.data <- fold.changes %>% 
  filter(group == 'Neutropenic') %>% 
  left_join(select(metab.metadata, CHEM_ID, SUPER_PATHWAY, SUB_PATHWAY, CHEMICAL_NAME), by = c('metab' = 'CHEM_ID')) %>% 
  left_join(loadings, by = c('metab' = 'CHEM_ID')) %>% 
  mutate(color = case_when(component == 1 ~ 'a', 
                           component == 2 ~ 'b',
                           is.na(component) & log2fc < 0 & p.value < 0.05 ~ 'c',
                           .default = 'd'),
         CHEMICAL_NAME = str_remove_all(CHEMICAL_NAME, '\\s[:punct:].{1,}$'))

plot <- ggplot(data = plot.data, aes(x = log2fc, y = -log10(p.value))) +
  
  geom_point(aes(color = color)) +
  
  geom_hline(yintercept = -log10(0.05), linetype = 'longdash') +
  geom_vline(xintercept = 0, linetype = 'longdash') +
  
  geom_text_repel(aes(label = ifelse(p.value < 0.05 & log2fc < 0, CHEMICAL_NAME, ''),
                      color = color),
                  size = 4,
                  force = 2,
                  max.iter = 20000,
                  max.overlaps = 20) +
  
  #annotate('text', x = 5, y = 3.8, label = 'More abundant at follow-up', fontface = 'bold') +
  #annotate('text', x = -5, y = 3.8, label = 'More abundant at baseline', fontface = 'bold') +
  
  labs(x = 'Log2 fold change',
       y = "-Log10 p-value (Quade test)",
       title = 'Neutropenic Subjects') +
  
  scale_x_continuous(limits = c(-8,8),
                     breaks = c(-8,-6,-4,-2,0,2,4,6,8)) +
  
  scale_y_continuous(limits = c(0,4)) +
  
  scale_color_manual(values = c('forestgreen', 'darkorange', 'firebrick', 'cornflowerblue')) +
  
  theme_classic() +
  
  theme(axis.title = element_text(size = 14, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        plot.title = element_text(size = 16, face = 'bold', hjust = 0.5),
        
        legend.position = 'none') 

svg('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/Volcano_plots/volcano_plot_prospective_neutropenic_unscaled_top_metabolites_annotated_20240612.svg', 
    height = 8, width = 8)

print(plot)

dev.off()

sig.results <- plot.data %>% filter(metab %in% neut.only) %>% arrange(p.value)

write_csv(sig.results,
          '//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/TINMAN_prospective_neutropenic_subjects_fold_changes_20240612.csv')

# Volcano plot: controls --------------------------------------------------

#' Annotated fold changes and p-values for plotting.
control.plot.data <- fold.changes %>% 
  filter(group == 'Control') %>% 
  left_join(select(metab.metadata, CHEM_ID, SUPER_PATHWAY, SUB_PATHWAY, CHEMICAL_NAME), by = c('metab' = 'CHEM_ID')) %>% 
  left_join(loadings, by = c('metab' = 'CHEM_ID')) %>% 
  mutate(color = case_when(component == 1 ~ 'a', 
                           component == 2 ~ 'b',
                           is.na(component) & log2fc < 0 & p.value < 0.05 ~ 'c',
                           .default = 'd'),
         CHEMICAL_NAME = str_remove_all(CHEMICAL_NAME, '\\s[:punct:].{1,}$'))


control.plot <- ggplot(data = control.plot.data, aes(x = log2fc, y = -log10(p.value))) +
  
  geom_point(aes(color = color)) +
  
  geom_hline(yintercept = -log10(0.05), linetype = 'longdash') +
  geom_vline(xintercept = 0, linetype = 'longdash') +
  
  geom_text_repel(aes(label = ifelse(p.value < 0.05 & log2fc < 0, CHEMICAL_NAME, 
                              ifelse(p.value < 1E-4 & log2fc > 0, CHEMICAL_NAME, '')),  
                      color = color),
                  size = 4,
                  force = 2,
                  max.iter = 20000,
                  max.overlaps = 20) +
  
  #annotate('text', x = 5, y = 3.8, label = 'More abundant at follow-up', fontface = 'bold') +
  #annotate('text', x = -5, y = 3.8, label = 'More abundant at baseline', fontface = 'bold') +
  
  labs(x = 'Log2 fold change',
       y = "-Log10 p-value (Quade test)",
       title = 'Non-neutropenic Subjects') +
  
  scale_x_continuous(limits = c(-8,8),
                     breaks = c(-8,-6,-4,-2,0,2,4,6,8)) +
  
  scale_y_continuous(limits = c(0,8)) +
  
  scale_color_manual(values = c('forestgreen', 'darkorange', 'firebrick', 'cornflowerblue')) +
  
  theme_classic() +
  
  theme(axis.title = element_text(size = 14, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'),
        plot.title = element_text(size = 16, face = 'bold', hjust = 0.5),
        
        legend.position = 'none') 

control.plot

# Compare to logistic regression results ----------------------------------

log.reg.results <- read_csv('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/prospective_mwas_results_unscaled_20240612.csv')

log.reg.results <- log.reg.results %>% 
  filter(P_VALUE_DELTA < 0.05 | P_VALUE_BASELINE < 0.05)

log.reg.results %>% filter(METAB %in% neut.only)
