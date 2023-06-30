# Prep environment --------------------------------------------------------
require(tidyverse); require(readxl)

setwd('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/')

#' Imputed, autoscaled metabolite data from the merged file.
metabs <- readRDS('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_merged_metab_imputed_autoscaled_20230606.rds')

#' Metadata concerning the metabolites above, notably names and accession numbers.
metab.metadata <- readRDS('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_merged_feces_metadata.rds') %>% 
  mutate(CHEMICAL_NAME = str_trim(tolower(CHEMICAL_NAME), 'both'))

#' P-values for 

#' Clinical metadata like group and timepoint.
clinical.metadata <- readRDS('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_merged_clinical_data.rds')

#' A list of pathways and metabolites that turned up in the most recent qMSEA.
msea.hits <- read_xlsx(path = 'MSEA_results/merged/MSEA_metabolite_hits_20230620.xlsx', sheet = 1) %>% 
  mutate(metabolites = tolower(metabolites),
         pathway = str_replace_all(pathway, ' ', '_'))

# Pull and plot metabolites -----------------------------------------------

#' A vector of the timepoints at which each pathway was found to differ between control & neutropenic subjects.
my.timepoints <- toupper(msea.hits$timepoint)

#' For every pathway, extract relevant metabolites and generate box plots of control versus neutropenic abundances.
#' Save these as .svg and .png files.
for (i in 1:nrow(msea.hits)){
  
  index.metabs <- msea.hits$metabolites[i] %>% str_split_1(',') %>% str_trim('both')
  
  ids <- metab.metadata %>% 
    filter(CHEMICAL_NAME %in% index.metabs) %>% 
    pull(CHEM_ID)
  
  #' Run t-tests on each metabolite so that we can annotate with the p-value.
  test.data <- metabs %>% 
    select(PARENT_SAMPLE_NAME, all_of(ids)) %>% 
    left_join(select(clinical.metadata, PARENT_SAMPLE_NAME, GROUP_NAME, TIMEPOINT)) %>% 
    filter(TIMEPOINT == my.timepoints[i])
  
  test.results <- tibble()
  
  for (j in seq_along(ids)){
    
    index.test.data <- select(test.data, GROUP_NAME, ids[j])
    
    #' Only metabolite abundance columns. Used to find y-axis value at which to annotate plots that are produced below.
    metabs.only <- test.data %>% 
      select(starts_with('C'))
    
    names(index.test.data) <- c('group','metab')
    
    test <- t.test(metab ~ group, data = index.test.data)
    
    test.result <- tibble(id = ids[j],
                          timepoint = my.timepoints[i],
                          t.test.p.value = paste0('p=', round(test$p.value, 2)),
                          value = max(metabs.only)*0.8) # 'value' will become the point on the y-axis at which the p-value will be plotted.
    
    test.results <- rbind(test.results, test.result)
    
  }
  
  test.results <- test.results %>% 
    mutate(name = index.metabs,
           GROUP_NAME = 1.5) # point along the x-axis at which the p-value will be plotted. This puts it between the two groups.
  
  plot.data <- metabs %>% 
    select(PARENT_SAMPLE_NAME, all_of(ids)) 
  
  names(plot.data) <- c('PARENT_SAMPLE_NAME', index.metabs)
  
  plot.data <- plot.data %>% 
    pivot_longer(cols = !PARENT_SAMPLE_NAME) %>% 
    left_join(select(clinical.metadata, PARENT_SAMPLE_NAME, GROUP_NAME, TIMEPOINT), by = 'PARENT_SAMPLE_NAME')
  
  index.plot.data <- plot.data %>% 
    filter(TIMEPOINT == my.timepoints[i])
  
  plot <- ggplot(data = index.plot.data, aes(x = GROUP_NAME, y = value)) + 
    
    geom_boxplot(aes(fill = GROUP_NAME)) +
    
    geom_text(aes(label = t.test.p.value), data = test.results, size = 5, fontface = 'bold') +
    
    theme_classic() + 
    
    labs(title = paste0(str_replace_all(toupper(msea.hits$pathway[i]), '_', ' '), ':\n', 'CONTROL VERSUS NEUTROPENIC', ' (', my.timepoints[i], '-TREATMENT', ')'), 
         x = '', 
         y = 'Autoscaled Metabolite abundance') +
    
    facet_wrap(~name) + 
    
    theme(legend.position = 'none',
          
          plot.title = element_text(face = 'bold', size = 12),
          
          strip.text = element_text(face = 'bold', size = 12),
          
          axis.text = element_text(face = 'bold', size = 12),
          axis.title = element_text(face = 'bold', size = 12),
          axis.title.y = element_text(margin = margin(0,10,0,0)))

  svg(paste0('Figures/MSEA/boxplot_',msea.hits$pathway[i],'_',format(Sys.Date(), '%Y%m%d'),'.svg'),
      8, 8)
  
  print(plot)
  
  dev.off()
  
  png(paste0('Figures/MSEA/boxplot_',msea.hits$pathway[i],'_',format(Sys.Date(), '%Y%m%d'),'.png'),
             width = 2400, height = 1600, res = 300)
  
  print(plot)
  
  dev.off()
  
}

# Scratch paper -----------------------------------------------------------
