########################################################################
#
# process_rnaseq_cellines
# process RNAseq data from Cell Lines. From raw data to voom and ComBat
# batch corrected.
#
# Copyright (C) 2019 saezlab, Uniklinikum Heidelberg
# Author: Luz Garcia-Alonso
# minor changes and exploratory analysis added by Rosa Hernansaiz-Ballesteros
# (https://github.com/saezlab/PathActivities/blob/master/src/transcriptomics/rnaseq_cell_lines/process_rnaseq_cellines.r)
#
# process_rnaseq_cellines.r
# Pipeline to process RNAseq data from a batch of 3 different datasets:
# GDSC from EGAS00001000828, CCLE from PRJNA169425 and Genentech from EGAS00001000610.
# Data available at https://rwth-aachen.sciebo.de/s/KhEshwnoFarWhKM?path=%2FLuz_data
# From the raw data, the pipeline:
# - Imputes 0 to NA values
# - Removes samples exeding genes with 0
# - Removes non-expressed genes with average CPM <=0
# - Normalises through TMM method
# - Runs Voom transformation
# - Uses ComBat for Batch correction of voom transformed data 
# - Merges duplicates by the mean
# 
# process_rnaseq_cellines is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# process_rnaseq_cellines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CARNIVAL1000.  If not, see <https://www.gnu.org/licenses/>.
#
########################################################################

rm(list=ls());cat("\014");if(length(dev.list())>0){dev.off()};options(warn=-1)

home = '~/Projects/process_rnaseq_cellines/GDSC_CLLE_GENENTECH'
setwd(home)



library(limma)
library(sva)
library(edgeR)
#library(GMD)
library(diptest)

library(ggplot2)
library(ggforce)
library(ggfortify)


cat('Load expression data \n')
load('~/Projects/00-Luz/merged_rawcounts.RData')
dim(EXP)
# [1] 21823  2001
# length(grep('CCLE',colnames(EXP)))
# GDSC samples: 447
# GENENTECH samples: 622
# CCLE samples: 932

cat('Load clinics data \n')
load('~/Projects/00-Luz/merged_cellcovariates.RData')
dim(CLI)
# [1] 2090   16



cat('Remove samples with unknown tissue\n')
tissue = names(which(table(CLI$gdsc_desc_1) > 3))
samples = intersect(colnames(EXP),  CLI$expression_matrix_name[ which(CLI$gdsc_desc_1 %in% tissue) ])
EXP_rawcounts = EXP[ , colnames(EXP) %in% samples ]
cat('Total genes and cell lines:', dim(EXP_rawcounts), '\n')
# Total genes and cell lines: 21823 1986 
# length(grep('GENENTECH',colnames(EXP_rawcounts)))
# GDSC samples: 447
# GENENTECH samples: 620
# CCLE samples: 919



cat('Imputing NA (ZEROmethod) \n') 
# There are no NA in these datasets
EXP_rawcounts[ which(is.na(EXP_rawcounts)) ] = 0
cat('Total genes and cell lines:', dim(EXP_rawcounts), '\n')

png('xplr/rawcounts_distribution.png', height = 800, width = 900, res = 150)
hist(log2(as.numeric(EXP_rawcounts)), breaks = 100, main = 'Cell lines raw counts', xlab = 'log2(raw counts)')
dev.off()

png('xplr/sequncing_depth.png', height = 800, width = 1000, res = 150)
par(mfrow = c(1,3))
barplot(colSums(EXP_rawcounts[,grep('GDSC',colnames(EXP_rawcounts))]), main = 'raw counts GDSC', xlab = 'Sum raw counts per Sample')
barplot(colSums(EXP_rawcounts[,grep('CCLE',colnames(EXP_rawcounts))]), main = 'raw counts CCLE', xlab = 'Sum raw counts per Sample')
barplot(colSums(EXP_rawcounts[,grep('GENENTECH',colnames(EXP_rawcounts))]), main = 'raw counts GENENTECH', xlab = 'Sum raw counts per Sample')
dev.off()



cat('Remove samples exeding genes with 0\n')
zero.count = apply(EXP_rawcounts, 2, function(x) sum(x==0) ) / nrow(EXP_rawcounts)
summary(zero.count)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1251  0.2148  0.2330  0.2376  0.2527  0.6088 

EXP_rawcounts = EXP_rawcounts[ , ! colnames(EXP_rawcounts) %in% names(zero.count)[ zero.count > .4 ] ]
cat('Total genes and cell lines:', dim(EXP_rawcounts), '\n')
# Total genes and cell lines: 21823 1985


cat('Revome excluded samples from CLI\n')
all(colnames(EXP) %in% CLI$expression_matrix_name)
CLI = unique(CLI[ CLI$expression_matrix_name %in% colnames(EXP) &  (CLI$Analysis_status != 'Exclude' | is.na(CLI$Analysis_status) ) , grep('^sample', names(CLI), invert = T) ])



message('Remove non-expressed genes with average CPM <= 0 ') # "The limma-voom method assumes that rows with zero or very low counts have been removed"
keep = aveLogCPM(EXP_rawcounts) > 0
EXP_rawcounts = EXP_rawcounts[ keep, ]
cat('Total genes and cell lines:', dim(EXP_rawcounts), '\n')
# Total genes and cell lines: 15388 1985 




cat('Normalization\n')
# Although it is also possible to give a matrix of counts directly to voom without TMM normalization, limma package reommends it
EXP_rawcounts = DGEList(counts=EXP_rawcounts,genes=rownames(EXP_rawcounts))
EXP_normcounts = calcNormFactors(EXP_rawcounts, method = 'TMM')
save(EXP_normcounts, file = 'tmm_norm.RData')
#     "method="TMM" is the weighted trimmed mean of M-values (to the reference) proposed by Robinson
#           and Oshlack (2010), where the weights are from the delta method on Binomial data. If refColumn
#           is unspecified, the library whose upper quartile is closest to the mean upper quartile is used.
#     method="RLE" is the scaling factor method proposed by Anders and Huber (2010). We call it
#           "relative log expression", as median library is calculated from the geometric mean of all columns
#           and the median ratio of each sample to the median library is taken as the scale factor.
#     method="upperquartile" is the upper-quartile normalization method of Bullard et al (2010), in
#           which the scale factors are calculated from the 75% quantile of the counts for each library, after
#           removing genes which are zero"





cat('Voom transformation\n')
EXP_voom = voom(EXP_normcounts, plot = T)$E # $E extract the "data corrected for variance"
save(EXP_voom, file = 'voom.RData')
# load('voom.RData')

#PCA
exp.pca <- prcomp(t(EXP_voom))
str(exp.pca)

exp.pca.variances <- ((exp.pca$sdev^2) / (sum(exp.pca$sdev^2)))*100

png('xplr/voom_PCA_barplot.png', height = 800, width = 900, res = 150)
barplot(exp.pca.variances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(exp.pca$sdev)), 
        ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,20))
dev.off()

png('xplr/voom_PCA.png', height = 800, width = 900, res = 150)
qplot(exp.pca$x[,1], exp.pca$x[,2], data=CLI[ match(colnames(EXP_voom), CLI$expression_matrix_name), ],
      color= project_code, geom = c('point'), label = project_code,#,'text'
      xlab=paste("PC1, ", round(exp.pca.variances[1], 2), "%"), 
      ylab=paste("PC2, ", round(exp.pca.variances[2], 2), "%"))
dev.off()




cat('Batch correction of voom transfromed data using ComBat\n')
cat('Using "gdsc_desc_1" and triplicates as covariates\n')
batch = CLI$project_code[ match(colnames(EXP_voom), CLI$expression_matrix_name) ]
cov_tissue = CLI$gdsc_desc_1[ match(colnames(EXP_voom), CLI$expression_matrix_name) ]
triplicates = names(which(sort(table(unique(CLI[, c('COSMIC_ID', 'project_code') ])$COSMIC_ID)) == 3))
cov_triplicates = sapply(triplicates, function(id) (CLI$COSMIC_ID[ match(colnames(EXP_voom), CLI$expression_matrix_name) ] == id) + 0)
cov_triplicates[ is.na(cov_triplicates) ] = 0
covariates = model.matrix(~cov_tissue+cov_triplicates)
EXP_corvoom = ComBat(EXP_voom, batch = batch, mod = covariates, par.prior = T, prior.plots = T)
save(EXP_corvoom, file = 'voom_batchcor.RData')
# load(file = 'voom_batchcor.RData')


## PCA

exp.pca <- prcomp(t(EXP_corvoom))
exp.pca.variances <- ((exp.pca$sdev^2) / (sum(exp.pca$sdev^2)))*100

png('xplr/PCAcomBat_barplot.png', height = 800, width = 900, res = 150)
barplot(exp.pca.variances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(exp.pca$sdev)), 
        ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,20))
dev.off()

png('xplr/PCAcomBat.png', height = 800, width = 900, res = 150)
qplot(exp.pca$x[,1], exp.pca$x[,2], data=CLI[ match(colnames(EXP_corvoom), CLI$expression_matrix_name), ],
      color= project_code, geom = c('point'), label = project_code,#,'text'
      xlab=paste("PC1, ", round(exp.pca.variances[1], 2), "%"), 
      ylab=paste("PC2, ", round(exp.pca.variances[2], 2), "%"))
dev.off()

###


cat('Merging duplicates with the mean')
samples = unique(unlist(lapply(strsplit(colnames(EXP_corvoom), '\\.'), function(x) paste(x[2], x[3], sep = '.') )))
genes = rownames(EXP_corvoom)
X = matrix(NA, nrow = length(genes), ncol = length(samples), dimnames = list(genes, samples))
for( s in samples){
  x = EXP_corvoom[ , grep(paste(s, '$', sep=''), colnames(EXP_corvoom)) ]
  if ( is.matrix(x) )
    x = apply(x, 1, mean)
  X[, s] = x
}
EXPmerged = X
save(EXPmerged, file = 'voom_batchcor_merged.RData')
