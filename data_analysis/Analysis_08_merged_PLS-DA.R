# Prepare data ------------------------------------------------------------

require(tidyverse); require(xlsx)

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

#' Metabolite metadata
metabs <- readRDS('TINMAN_merged_feces_metadata.rds') %>% 
  #filter(SUPER_PATHWAY != 'Xenobiotics' | (SUPER_PATHWAY == 'Xenobiotics' & SUB_PATHWAY == 'Bacterial/Fungal')) %>% 
  pull(CHEM_ID)

#' Concentration table.
#' Choice of scaled or unscaled data does not impact the algorithm because it performs its own centering and scaling.
tinman <- readRDS('TINMAN_merged_metab_imputed_unscaled_20230606.rds') %>% 
  select(PARENT_SAMPLE_NAME, all_of(metabs))

#' Load and append clinical metadata.
clinical.metadata <- readRDS('TINMAN_merged_clinical_data.rds')

tinman <- tinman %>% 
  left_join(select(clinical.metadata, PARENT_SAMPLE_NAME:GROUP_NAME), by = 'PARENT_SAMPLE_NAME') %>% 
  select(PARENT_SAMPLE_NAME, GROUP_ID, GROUP_NAME, all_of(metabs))

# Matrix of metabolite abundances (x) and class memberships (y) needs for sparse PLS-DA.
x <- as.matrix(tinman[,4:ncol(tinman)])
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

MyResult.splsda.final <- splsda(x, y, ncomp = 2, keepX = c(50,50))

my.plot <- plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
                     ellipse = TRUE, title="sPLS-DA - final result",
                     group = y, 
                     col.per.group = c('red','white','blue','ghostwhite'),
                     pch = c(16,1,15,2),
                     cex = c(2,-2,2,-2))

my.plot

#' Recover and plot loading vectors for components 1 and 2
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean', comp = 1)
comp1.loadings <- selectVar(MyResult.splsda.final, comp=1)$value %>% 
  mutate(CHEM_ID = rownames(.))

plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean', comp = 2)
comp2.loadings <- selectVar(MyResult.splsda.final, comp=2)$value %>% 
  mutate(CHEM_ID = rownames(.))

#' Save plots
svg('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/PLS-DA/merged.plsa.da.plot.20230614.svg',
    width = 10, height = 10)

plotIndiv(MyResult.splsda.final, 
          ind.names = FALSE, 
          legend=TRUE,
          ellipse = TRUE, 
          title="sPLS-DA - final result")

dev.off()

svg('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/PLS-DA/merged.plsa.da.comp1.loadings.svg',
    width = 10, height = 10)

plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean', comp = 1)

dev.off()

svg('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/PLS-DA/merged.plsa.da.comp2.loadings.svg',
    width = 10, height = 10)

plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean', comp = 2)

dev.off()

#' Metabolite metadata
metabs <- readRDS('TINMAN_merged_feces_metadata.rds')

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
