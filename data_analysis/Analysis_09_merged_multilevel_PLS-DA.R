require(mixOmics)

#' This is a worked example of a multilevel SPLS procedure from the mixOmics authors.
#' See: http://mixomics.org/methods/multilevel/
#' and: http://mixomics.org/case-studies/multilevel-liver-toxicity-case-study/
data("liver.toxicity")

X <- liver.toxicity$gene # use the gene expression data as the X matrix
Y <- liver.toxicity$clinic # use the clinical data as the Y matrix

repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                  6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                  10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                  13, 14, 15, 16, 15, 16, 15, 16, 15, 16)

design <- data.frame(sample = repeat.indiv) # load this into a dataframe

summary(as.factor(repeat.indiv)) # 16 rats, 4 measurements each

spls.liver.multilevel <- spls(X, Y, 
                              ncomp = 2, 
                              keepX = c(20,50),
                              keepY = c(5, 10),
                              multilevel = design,
                              mode = 'canonical')

plotIndiv(spls.liver.multilevel, rep.space = "XY", 
          group = liver.toxicity$treatment$Time.Group, 
          pch = as.factor(liver.toxicity$treatment$Dose.Group),
          col.per.group = color.mixo(1:4), 
          legend = TRUE, legend.title = 'Time', legend.title.pch = 'Dose')

# My data -----------------------------------------------------------------

#' Concentration table.
tinman <- readRDS('Datasets/TINMAN_merged_metab_imputed_autoscaled_20230606.rds') 

sort.order <- tinman$PARENT_SAMPLE_NAME

#' Load and append clinical metadata.
clinical.metadata <- readRDS('Datasets/TINMAN_merged_clinical_data.rds') %>% 
  arrange(sort.order)

#' This is me working through the same procedure with my data (and a categorical outcome/classification task).
my.names <- tinman$PARENT_SAMPLE_NAME
my.x <- tinman[, 2:ncol(tinman)] %>% as.data.frame()
rownames(my.x) <- my.names

my.y <- as.factor(clinical.metadata$TIMEPOINT)
#my.y <- as.factor(clinical.metadata$GROUP_ID)

my.design <- data.frame(samples = rep(1:39, each = 2))

splsda <- splsda(X= my.x, Y = my.y, ncomp = 2, keepX = c(50,50))

svg(paste0('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/PLS-DA/PLS-DA_pre_versus_post_', format(Sys.Date(), '%Y%m%d'), '.svg'),
    width = 10, height = 10)

plotIndiv(splsda, 
          ind.names = FALSE, 
          legend=TRUE,
          ellipse = TRUE, 
          title="PLS-DA of timepoint on TINMAN data")

dev.off()

splsda.multilevel <- splsda(X= my.x, Y = my.y, ncomp = 2, keepX = c(50,50), multilevel = my.design)

svg(paste0('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/PLS-DA/multilevel_PLS-DA_pre_versus_post_', format(Sys.Date(), '%Y%m%d'), '.svg'),
    width = 10, height = 10)

my.pch <- clinical.metadata$GROUP_NAME
my.pch <- ifelse(my.pch == 'Control', 16, 3)

plotIndiv(splsda.multilevel, 
          ind.names = FALSE, 
          legend=TRUE,
          ellipse = TRUE, 
          pch = my.pch,
          title="Multilevel PLS-DA of timepoint on TINMAN data")

dev.off()