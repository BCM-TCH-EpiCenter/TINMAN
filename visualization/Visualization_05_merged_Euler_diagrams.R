require(tidyverse); require(eulerr)

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

results <- readRDS('TINMAN_merged_metabolite_detection_counts_20230608.rds')

# Identify metabolites detected only pre-treatment ------------------------

pre.transplant.controls <- results %>% 
  filter(control.n.pre > 0 & control.n.post == 0) %>% 
  pull(metabolite)

pre.transplant.neutropenia <- results %>% 
  filter(neutropenic.n.pre > 0 & neutropenic.n.post == 0) %>% 
  pull(metabolite)

pre.transplant.common <- intersect(pre.transplant.controls, pre.transplant.neutropenia)

# Identify metabolites detected only post-transplant ----------------------

post.transplant.controls <- results %>% 
  filter(control.n.post > 0 & control.n.pre == 0) %>% 
  pull(metabolite)

post.transplant.neutropenia <- results %>% 
  filter(neutropenic.n.post > 0 & neutropenic.n.pre == 0) %>% 
  pull(metabolite)

post.transplant.common <- intersect(post.transplant.controls, post.transplant.neutropenia)

# Euler diagrams ----------------------------------------------------------

#' Input data.
#' Numbers come from the lengths of the vectors above.
pre.euler <- c('Control' = length(pre.transplant.controls), 
               'Neutropenic' = length(pre.transplant.neutropenia), 
               'Control&Neutropenic' = length(pre.transplant.common))

post.euler <- c('Control' = length(post.transplant.controls), 
                'Neutropenic' = length(post.transplant.neutropenia), 
                'Control&Neutropenic' = length(post.transplant.common))

setwd('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/Euler_diagrams/')

#' Save plots.
svg('Euler_diagram_species_detected_pre_TX_only_20221012.svg', height = 8, width = 8)

plot(euler(pre.euler), quantities = T, main = 'Metabolites Detected Only at Baseline')

dev.off()

svg('Euler_diagram_species_detected_post_TX_only_20221012.svg', height = 8, width = 8)

plot(euler(post.euler), quantities = T, main = 'Metabolites Detected Only at Follow-up')

dev.off()

# Count species detected in each group ------------------------------------

count.metabolites <- function(x){
  
  out <- results %>% filter(.data[[x]] > 0) %>% pull(metabolite)
  
  return(out)
  
}

control.post <- count.metabolites('control.n.post')
control.pre <- count.metabolites('control.n.pre')
control.intersect <- intersect(control.post, control.pre)
control.euler <- c('Baseline' = length(control.pre) - length(control.intersect), 
                   'Follow-up' = length(control.post) - length(control.intersect), 
                   'Baseline&Follow-up' = length(control.intersect))

neutropenic.post <- count.metabolites('neutropenic.n.post')
neutropenic.pre <- count.metabolites('neutropenic.n.pre')
neutropenic.intersect <- intersect(neutropenic.post, neutropenic.pre)
neutropenic.euler <- c('Baseline' = length(neutropenic.pre) - length(neutropenic.intersect), 
                       'Follow-up' = length(neutropenic.post) - length(neutropenic.intersect), 
                       'Baseline&Follow-up' = length(neutropenic.intersect))

baseline.intersect <- intersect(control.pre, neutropenic.pre)
baseline.euler <- c('Control' = length(control.pre)-length(baseline.intersect), 
                    'Neutropenic' = length(neutropenic.pre)-length(baseline.intersect), 
                    'Control&Neutropenic' = length(baseline.intersect))

follow.up.intersect <- intersect(control.post, neutropenic.post)
follow.up.euler <- c('Control' = length(control.post)-length(follow.up.intersect), 
                     'Neutropenic' = length(neutropenic.post)-length(follow.up.intersect), 
                     'Control&Neutropenic' = length(follow.up.intersect))

setwd('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/Euler_diagrams/')

#' Save plots.
svg('Euler_diagram_control_pre_post_20230608.svg', height = 8, width = 8)

plot(venn(control.euler), quantities = T, main = 'Metabolites Detected at Baseline and Follow-up (Controls)')

dev.off()

svg('Euler_diagram_neutropenic_pre_post_20230608.svg', height = 8, width = 8)

plot(venn(neutropenic.euler), quantities = T, main = 'Metabolites Detected at Baseline and Follow-up (Neutropenic Subjects)')

dev.off()

svg('Euler_diagram_baseline_20230608.svg', height = 8, width = 8)

plot(venn(baseline.euler), quantities = T, main = 'Richness of the Control & Neutropenic Metabolomes (Baseline)')

dev.off()

svg('Euler_diagram_follow_up_20230608.svg', height = 8, width = 8)

plot(venn(follow.up.euler), quantities = T, main = 'Richness of the Control & Neutropenic Metabolomes (Follow-up)')

dev.off()
