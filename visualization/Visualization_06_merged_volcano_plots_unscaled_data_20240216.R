require(tidyverse); require(ggrepel)

# Load and prepare data ---------------------------------------------------

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/')

#' Metabolite metadata
metab.metadata <- readRDS('Datasets/TINMAN_merged_feces_metadata.rds') 

metabs <- metab.metadata %>% 
  filter(SUPER_PATHWAY != 'Xenobiotics' | (SUPER_PATHWAY == 'Xenobiotics' & SUB_PATHWAY == 'Bacterial/Fungal')) %>% 
  pull(CHEM_ID)

#' Concentration table.
tinman <- readRDS('Datasets/TINMAN_merged_metab_imputed_unscaled_20230606.rds') 

sort.order <- tinman$PARENT_SAMPLE_NAME

#' Load and append clinical metadata.
clinical.metadata <- readRDS('Datasets/TINMAN_merged_clinical_data.rds') %>% 
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

# Annotate neutropenic results --------------------------------------------

#' The metabolites Josie would like annotated.
neut.hits <- fold.changes %>% 
  filter(group == 'Neutropenic', p.value < 0.05, log2fc < 0) %>% 
  pull(metab)

#con.hits <- fold.changes %>% 
#  filter(group == 'Control', p.value < 0.05) %>% 
#  pull(metab)

#' Compounds that are only altered in neutropenic patients.
#difference <- setdiff(neut.hits, con.hits)

#' Compounds that are altered in both groups, but in different directions.
#discordant <- fold.changes %>% 
#  filter(metab %in% intersect(neut.hits, con.hits)) %>% 
#  select(metab, group, log2fc) %>% 
#  pivot_wider(id_cols = metab, names_from = group, values_from = log2fc) %>% 
#  mutate(sign.Control = sign(Control),
#         sign.Neutropenic = sign(Neutropenic),
#         discordant = case_when(sign.Control == sign.Neutropenic ~ 0,
#                                sign.Control != sign.Neutropenic ~ 1)) %>% 
#  filter(discordant == 1) %>% 
#  pull(metab)

#' The union of those two sets.
#my.hits <- c(difference, discordant)

#' Append metabolite metadata for the metabolites above.
#neut.hits <- fold.changes %>%
#  left_join(select(metab.metadata, CHEM_ID, SUPER_PATHWAY, SUB_PATHWAY, CHEMICAL_NAME, HMDB), 
#            by = c('metab' = 'CHEM_ID')) %>% 
#  filter(metab %in% neut.hits, group == 'Neutropenic') %>% 
#  arrange(p.value) %>% 
#  select(group, metab, CHEMICAL_NAME, HMDB, contains('PATHWAY'), POST:p.value)

#' Loading scores for metabolites in PLS-DA components 1 and 2.
loadings <- read_csv('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/plsa_da_component_loadings_20240111.csv')

#' Filter to those in our set of interest.
loadings <- loadings %>% 
  select(-CHEMICAL_NAME) %>% 
  filter(CHEM_ID %in% neut.hits)

#neut.hits <- neut.hits %>% 
#  left_join(loadings, by = c('metab' = 'CHEM_ID')) %>% 
#  mutate(pls.da.hit = ifelse(is.na(component), 'No', 'Yes')) %>% 
#  select(group:p.value, pls.da.hit, component:rank)

#' Files containing the names and HMDB IDs of compounds annotated to pathways that came up in MetaboAnalyst ORA.
#pathways <- list(
#    sphingolipid = read_csv('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/SMPDB_pathways/SMP0000034_metabolites.csv'),
#    pyrimidine = read_csv('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/SMPDB_pathways/SMP0000046_metabolites.csv'),
#    butyrate = read_csv('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/SMPDB_pathways/SMP0000073_metabolites.csv'),
#    thiamine = read_csv('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/SMPDB_pathways/SMP0000076_metabolites.csv'),
#    pantothenate = read_csv('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/SMPDB_pathways/SMP0000027_metabolites.csv')
#  )

#pathway.metabolites <- lapply(pathways, '[', , 'HMDB ID') %>% unlist()
#pathway.names <- lapply(pathways, '[', , 'Pathway Name') %>% unlist()

#pathways <- data.frame(metab = pathway.metabolites,
#                       ORA.pathway = pathway.names)

#neut.hits <- neut.hits %>% 
#  mutate(HMDB = str_remove(HMDB, ',[:alnum:]{1,}')) %>% 
#  left_join(pathways, by = c('HMDB' = 'metab')) %>% 
#mutate(ora.hit = ifelse(is.na(ORA.pathway), 0, 1)) %>% 
#  select(group:rank, ora.hit, ORA.pathway)

#write_csv(neut.hits, file = '//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/TINMAN_neutropenic_hits_20240220.csv')

# Generate volcano plot ---------------------------------------------------

#' These are Metabolon metabolite accession numbers for a handful of species that Josie indicated she wanted highlighted in an email on 6/21/23.
#' (This was later updated to exclude creatine).
#my.metabs <- c('C391','C100001104', 'C100006435', 'C100000098', 'C925', 'C926', 'C100003239', 'C100001416', 'C100021504')

#' Annotated fold changes and p-values for plotting.
plot.data <- fold.changes %>% 
  filter(group == 'Neutropenic') %>% 
  left_join(select(metab.metadata, CHEM_ID, SUPER_PATHWAY, CHEMICAL_NAME), by = c('metab' = 'CHEM_ID')) %>% 
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

svg('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/Volcano_plots/volcano_plot_neutropenic_unscaled_top_metabolites_annotated_20240227.svg', 
    height = 8, width = 8)

print(plot)

dev.off()

sig.results <- plot.data %>% filter(p.value < 0.05) %>% arrange(p.value)
write_csv(sig.results,
          '//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/TINMAN_neutropenic_subjects_fold_changes_20240216.csv')
