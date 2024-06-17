require(tidyverse); require(mixOmics); require(xlsx)

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

#' Concentration table
tinman <- readRDS('TINMAN_feces_20221129_2.rds')

#' Metabolite metadata
metabs <- readRDS('TINMAN_feces_metadata.rds') %>% 
  filter(SUPER_PATHWAY != 'Xenobiotics' | (SUPER_PATHWAY == 'Xenobiotics' & SUB_PATHWAY == 'Bacterial/Fungal')) %>% 
  pull(CHEM_ID)

#' Remove xenobiotics, except for bacterial/fungal metabolites..
tinman <- tinman %>% dplyr::select(PARENT_SAMPLE_NAME:TIMEPOINT, all_of(metabs))

# Matrix of metabolite abundances (x) and class memberships (y) needs for sparse PLS-DA.
x <- as.matrix(tinman[,6:ncol(tinman)])
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

MyResult.splsda.final <- splsda(x, y, ncomp = 2, keepX = c(50,50))

my.plot <- plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
                     ellipse = TRUE, title="sPLS-DA - final result")

#' Recover and plot loading vectors for components 1 and 2
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean', comp = 1)
comp1.loadings <- selectVar(MyResult.splsda.final, comp=1)$value %>% 
  mutate(CHEM_ID = rownames(.))

plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean', comp = 2)
comp2.loadings <- selectVar(MyResult.splsda.final, comp=2)$value %>% 
  mutate(CHEM_ID = rownames(.))

#' Metabolite metadata
metabs <- readRDS('TINMAN_feces_metadata.rds')

#' Append metadata to loadings
comp1.loadings <- left_join(comp1.loadings, metabs, by = 'CHEM_ID')

comp2.loadings <- left_join(comp2.loadings, metabs, by = 'CHEM_ID')

#' Evaluate overlap between loadings and univariate results
uni <- read.xlsx('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/TINMAN_metabolomics_analysis_20221017.xlsx',
                 sheetIndex = 2)

uni <- uni %>% filter(t.test.p < 0.05 | kruskal.test.p < 0.05)

loadings <- c(comp1.loadings$CHEM_ID, comp2.loadings$CHEM_ID)

intersect(uni$CHEM_ID, loadings)

# Save the loadings -------------------------------------------------------

write.xlsx(comp1.loadings,
           '//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/TINMAN_PLS-DA_loadings_20221212.xlsx',
           sheetName = 'component1',
           row.names = F)

write.xlsx(comp2.loadings,
           '//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/TINMAN_PLS-DA_loadings_20221212.xlsx',
           sheetName = 'component2',
           row.names = F,
           append = T)

#auc.plsda <- auroc(MyResult.splsda.final)

# Save plots --------------------------------------------------------------

#' Save plots
svg('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/plsa.da.plot.20221212.svg',
    width = 10, height = 10)

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - final result")

dev.off()

svg('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/plsa.da.comp1.loadings.svg',
    width = 10, height = 10)

plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean', comp = 1)

dev.off()

#' Save plots
svg('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/plsa.da.comp2.loadings.svg',
    width = 10, height = 10)

plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean', comp = 2)

dev.off()

# Calculate pre-post change in top metabs ---------------------------------

#' For the 100 metabolites that have the highest loadings, calculate pre- to post- change.
tinman.wide <- tinman %>% 
  dplyr::select(CLIENT_IDENTIFIER, GROUP_NAME, TIMEPOINT, all_of(loadings)) %>% 
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

#' Until we reach the end of these sets, calculate post - pre fold change.
while (j <= ncol(tinman.wide)){
  
  my.var <- names(tinman.wide)[i] %>% str_remove('_[:alpha:]{3,}')
  
  #' Subtract pre-induction values from post-induction values to get the delta value. 
  my.values <- (tinman.wide[,i] - tinman.wide[,j]) %>% as.vector() %>% unlist()
  
  #' Assign these to a new column.
  tinman.wide[[my.var]] <- my.values
  
  i <- i + 2
  j <- j + 2
  
}

#' Remove columns holding post-treatment abundances. Retain pre-treatment abundances since these will be covariates in models.
#' Convert GROUP_NAME to factor.
tinman.wide <- tinman.wide %>% 
  dplyr::select(!contains('POST')) %>% 
  mutate(GROUP_NAME = factor(GROUP_NAME))

tmp <- tinman.wide %>% as.data.frame()

# Test changes for association with neutropenia ---------------------------

#' Empty data frame to hold results.
results <- data.frame()

for (i in seq_along(loadings)){
  
  tmp <- tinman.wide %>% dplyr::select(GROUP_NAME, contains(loadings[i]))
  names(tmp) <- c('group', 'baseline', 'change')
  
  if (var(tmp$baseline) > 0) {
    
    
    my.model <- glm(group ~ baseline + change, data = tmp, family = 'binomial') %>% summary()
    
    new.result <- data.frame(metab = loadings[i],
                             coef.delta = my.model$coefficients[3,1],
                             se.delta = my.model$coefficients[3,2],
                             p.value.delta = my.model$coefficients[3,4])
    
    results <- rbind(results, new.result)
    
  }
  
  else {
    
    next
    
  }
  
}

# Multilevel PLS-DA -------------------------------------------------------

#' A tutorial on multilevel PLS-DA here:
#' http://mixomics.org/methods/multilevel/
#' This is appropriate when you have a repeated measures design and the variation by subject may be greater than the variation by treatment/exposure.
