# Prepare environment -----------------------------------------------------

require(tidyverse); require(readxl)
require(ggrepel); require(plotly); require(htmlwidgets)

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

#' A function to perform t tests and Kruskal-Wallis tests comparing metabolite abundances pre- vs. post- separately in control and neutropenic subjects.
univariate.tests <- function(x){
  
  #' Initialize an empty data frame to hold results for pre-/post- comparisons in the control and neutropenic groups.
  results <- data.frame()
  
  #' For each metabolite...
  for (i in seq_along(metabolites)){
    
    print(i)
    
    data <- select(x, GROUP_NAME, TIMEPOINT, all_of(metabolites[i]))
    names(data) <- c('GROUP_NAME','TIMEPOINT','METABOLITE')
    data <- data %>% filter(!is.na(METABOLITE))
    
    counts <- data %>% count(TIMEPOINT, .drop = F)
    
    #' If it was detected in more than one subject from the index group at each timepoint, calculate the means at baseline and follow up and compare by univariate tests.
    if (nrow(counts) == 2 & all(counts$n > 1) & var(data$METABOLITE) != 0){  
      
      t.test <- data %>% 
        t.test(METABOLITE ~ TIMEPOINT, data = .)
      
      kruskal.test <- data %>% 
        kruskal.test(METABOLITE ~ TIMEPOINT, data = .)
      
      result <- tibble(metabolite = metabolites[i],
                       mean.pre = t.test$estimate[2],
                       mean.post = t.test$estimate[1],
                       t.test.p = t.test$p.value,
                       kruskal.test.p = kruskal.test$p.value)
      
      results <- rbind(results, result)
      
    }
    
    #' If that metabolte was not detected in the index group, mark as missing.
    else if (nrow(data) == 0){
      
      result <- tibble(metabolite = metabolites[i],
                       mean.pre = NA,
                       mean.post = NA,
                       t.test.p = NA,
                       kruskal.test.p = NA)
      
      results <- rbind(results, result)
      
    }
    
    #' Otherwise, report means at baseline and follow up.
    else {
      
      means <- aggregate(METABOLITE ~ TIMEPOINT, data = data, mean, drop = F)
      
      result <- tibble(metabolite = metabolites[i],
                       mean.pre = means[2,2],
                       mean.post = means[1,2],
                       t.test.p = NA,
                       kruskal.test.p = NA)
      
      results <- rbind(results, result)
      
    }
    
  }
  
  return(results)
  
}

#' A function that modifies a data frame prior to plotting.
prep.data <- function(data, string){
  
  new <- data %>% 
    select(metabolite, CHEMICAL_NAME, starts_with(string)) %>% 
    rename_with(~str_remove(.x, paste0(string,'.')), starts_with(paste0(string,'.'))) %>% 
    filter(n.pre > 1, n.post > 1) %>% 
    mutate(log2.fold.change = log(mean.post/mean.pre, base = 2),
           color = factor(ifelse(t.test.p >= 0.05, 0,
                          ifelse(t.test.p < 0.05 & abs(log2.fold.change) > 2, 2, 1)),
                               labels = c('blue', 'red', 'green')))
  
  return(new)
  
}

#' A function to generate a customized volcano plot.
volcano <- function(data, title, cutoff){
  
  plot <- ggplot(data = data, aes(x = log2.fold.change, y = -log10(t.test.p))) +
    
    geom_hline(yintercept = -log10(0.05), linetype = 'longdash') +
    
    geom_vline(xintercept = -2, linetype = 'longdash') +
    geom_vline(xintercept = 2, linetype = 'longdash') +
    
    geom_point(aes(color = color)) +
    
    
    geom_text_repel(aes(label = ifelse(t.test.p < cutoff, CHEMICAL_NAME, '')),
                    max.iter = 20000) +
    
    annotate('text', x = 5, y = 3.8, label = 'More abundant at follow-up', fontface = 'bold') +
    annotate('text', x = -5, y = 3.8, label = 'More abundant at baseline', fontface = 'bold') +
    
    labs(x = 'Log2 fold change (follow-up/baseline, unscaled values)',
         y = "-Log10 p-value (Student's t-test)",
         title = title) +
    
    scale_x_continuous(limits = c(-8,8),
                       breaks = c(-8,-6,-4,-2,0,2,4,6,8)) +
    
    scale_y_continuous(limits = c(0,4)) +
    
    scale_color_manual(values = c('cornflowerblue','forestgreen','firebrick1')) +
    
    theme_classic() +
    
    theme(axis.title = element_text(size = 14, face = 'bold'),
          axis.text = element_text(size = 12, face = 'bold'),
          plot.title = element_text(size = 16, face = 'bold'),
          
          legend.position = 'none')
  
  return(plot)
  
}

#' A variant of the function above that produces an interactive plot at the cost of the color aesthetic.
volcano.plotly <- function(data, title){
  
  plot <- ggplot(data = data, aes(x = log2.fold.change, y = -log10(t.test.p), color = CHEMICAL_NAME)) +
    
    geom_hline(yintercept = -log10(0.05), linetype = 'longdash') +
    
    geom_vline(xintercept = -2, linetype = 'longdash') +
    geom_vline(xintercept = 2, linetype = 'longdash') +
    
    geom_point() +
    
    annotate('text', x = 5, y = 3.8, label = 'More abundant at follow-up', fontface = 'bold') +
    annotate('text', x = -5, y = 3.8, label = 'More abundant at baseline', fontface = 'bold') +
    
    labs(x = 'Log2 fold change (follow-up/baseline, unscaled values)',
         y = "-Log10 p-value (Student's t-test)",
         title = title) +
    
    scale_x_continuous(limits = c(-8,8),
                       breaks = c(-8,-6,-4,-2,0,2,4,6,8)) +
    
    scale_y_continuous(limits = c(0,4)) +
    
    scale_color_manual(values = rep('cornflowerblue', times = nrow(data))) +
    
    theme_classic() +
    
    theme(axis.title = element_text(size = 14, face = 'bold'),
          axis.text = element_text(size = 12, face = 'bold'),
          plot.title = element_text(size = 16, face = 'bold'),
          
          legend.position = 'none')
  
  return(plot)
  
}

#' Unscaled, imputed data will be used for purposes of volcano plots.
tinman <- readRDS('TINMAN_merged_metab_imputed_unscaled_20230606.rds')

#' File with clinical metadata for merged dataset.
clinical.metadata <- readRDS('TINMAN_merged_clinical_data.rds')

#' File with metabolite names.
metadata <- readRDS('TINMAN_merged_feces_metadata.rds') 

#' Metabolites to test.
metabolites <- metadata %>% 
  #filter(SUPER_PATHWAY != 'Xenobiotics' | (SUPER_PATHWAY == 'Xenobiotics' & SUB_PATHWAY == 'Bacterial/Fungal') ) %>% 
  pull(CHEM_ID)

#' Append group information.
tinman <- tinman %>% 
  left_join(select(clinical.metadata, PARENT_SAMPLE_NAME, GROUP_ID, GROUP_NAME, TIMEPOINT), by = 'PARENT_SAMPLE_NAME')

# Generate fold change data -----------------------------------------------

#' Split into control and neutropenic datasets.
control <- filter(tinman, GROUP_NAME == 'Control')

neutropenic <- filter(tinman, GROUP_NAME == 'Neutropenic')

#' Perform pre- vs. post- tests in each set, then bind them together.
control.results <- univariate.tests(control) %>% 
  rename_with(~ paste0('control.', .x), !metabolite)

neutropenic.results <- univariate.tests(neutropenic) %>% 
  rename_with(~ paste0('neutropenic.', .x), !metabolite)

#' Combine and append counts.
results <- left_join(control.results, neutropenic.results, by = 'metabolite')

counts <- readRDS('TINMAN_merged_metabolite_detection_counts_20230614.rds')

results <- results %>% left_join(counts, by = 'metabolite')

#' Add compound names and split by status.
results <- results %>% 
  left_join(select(metadata, CHEM_ID, CHEMICAL_NAME), by = c('metabolite' = 'CHEM_ID'))

control <- prep.data(results, 'control')

neutropenic <- prep.data(results, 'neutropenic')

# Static .SVG plots -------------------------------------------------------

setwd('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/Volcano_plots/')

#' .svg versions of the static plots.
svg('volcano_plot_controls_unscaled_20230614.svg', height = 8, width = 8)

volcano(control, 'CONTROLS', 0.005)

dev.off()

svg('volcano_plot_neutropenic_unscaled_20230614.svg', height = 8, width = 8)

volcano(neutropenic, 'NEUTROPENIC SUBJECTS', 0.025)

dev.off()rm

# Interactive .HTML plots -------------------------------------------------

#' .html versions of the interactive plots.
pc <- volcano.plotly(control, 'CONTROLS') %>% 
  ggplotly() %>% 
  saveWidget(file = 'interactive_volcano_plot_controls_unscaled_20230614.html',
             selfcontained = T)

pn <- volcano.plotly(neutropenic, 'NEUTROPENIC SUBJECTS') %>% 
  ggplotly() %>% 
  saveWidget(file = 'interative_volcano_plot_neutropenic_unscaled_20230614.html',
             selfcontained = T)
