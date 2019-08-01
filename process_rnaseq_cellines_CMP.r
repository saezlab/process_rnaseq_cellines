########################################################################
#
# process_rnaseq_cellines
# process RNAseq data from Cell Lines. From raw data to voom and ComBat
# batch corrected.
#
# Copyright (C) 2019 saezlab, Uniklinikum Heidelberg
# Author: Rosa Hernansaiz-Ballesteros
# addapted verion of process_rnaseq_cellines.r for CMP releases
# (https://github.com/saezlab/PathActivities/blob/master/src/transcriptomics/rnaseq_cell_lines/process_rnaseq_cellines.r)
#
# process_rnaseq_cellines_CMP.r
# Pipeline to process RNAseq data from Cell Model Pasports releases
# (https://cellmodelpassports.sanger.ac.uk/)
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

library(data.table)
library(dplyr)
library(stringr)
library(limma)
library(sva)
library(edgeR)

library(ggplot2)
library(ggforce)
library(ggfortify)



###### set paths and variables ######

wd = '~/Projects/process_rnaseq_cellines/CMP' #working directory
setwd(wd)

metasample = '../../00-CMP/metadata/model_list_2019-04-04_0916.csv' #vector containing the path to the metainformation about the cell lines
# metagene = '../00-CMP/metadata/gene_identifiers_2019-02-19_1024.csv' #vector containing the path to the metainformation about the genes used in the cell model passports, with HUGO, Ensembl, Entrez and Refseq annotation
rnadata = '../../00-CMP/expression/rnaseq_2019-04-15_1133.csv'#vector containing the path to the RNAseq file of cell lines




###### begin of the preprocessing ######
  
cat('Load expression data \n')
cmpRAW <- fread(rnadata, sep = ',')

cat('Load clinics data \n')
metaSMP <- read.delim(metasample, sep = ',')
# genes <- read.delim(metagene, sep = ',')




message('Create matrix Gene x Sample \n')

models <- unique(cmpRAW$model_id) #unique samples

## as the original dataset from Cell Model Pasports contains samples from 2 different sources (Sanger and Broad),
## the samples are separarated based on the source to run the pipeline, as they have different number of genes:

pb <- txtProgressBar(min = 0, max = length(models), style = 3) # create progress bar
for (i in 1:length(models)) { #for each of the unique samples
  setTxtProgressBar(pb, i) #update progress bar
  model <- models[i]
  df <-  cmpRAW[which(cmpRAW$model_id==model),c('gene_id','read_count', 'dataset_name')] #get the relevant information
  dset <- unique(df$dataset_name)
  for(institute in dset){
    dfsub <- df[which(df$dataset_name==institute),c('gene_id', 'read_count')]
    colnames(dfsub) <- c('gene_id', model)
    if(is.null(EXP_rawcounts[[institute]])){
      EXP_rawcounts[[institute]] <- dfsub
    }else{
      EXP_rawcounts[[institute]] <- merge(EXP_rawcounts[[institute]],dfsub, by = 'gene_id')
    }
  }
}
close(pb)

rm(pb, df, dfsub, i, dset, institute, model)

cat('Total genes and cell lines: \n')
print(lapply(EXP_rawcounts, dim))
# $`Sanger RNASeq`
# [1] 35004   448
# 
# $`Broad (CCLE) RNASeq`
# [1] 37260   707

# transforming data.table to data.frame 
names(EXP_rawcounts) <- gsub("Sanger RNASeq", 'sanger', names(EXP_rawcounts))
names(EXP_rawcounts) <- gsub("Broad \\(CCLE\\) RNASeq", 'broad', names(EXP_rawcounts))

EXP_rawcounts[['sanger']] <- setDF(EXP_rawcounts[['sanger']], rownames=EXP_rawcounts[['sanger']]$gene_id)
EXP_rawcounts[['sanger']] <- EXP_rawcounts[['sanger']][,-1]

EXP_rawcounts[['broad']] <- setDF(EXP_rawcounts[['broad']], rownames=EXP_rawcounts[['broad']]$gene_id)
EXP_rawcounts[['broad']] <- EXP_rawcounts[['broad']][,-1]

saveRDS(EXP_rawcounts, file = '1907_release1904/rnaseq_2019-06-26_geneVSmodel_list.rds')
# EXP_rawcounts <- readRDS('1907_release1904/rnaseq_2019-06-26_geneVSmodel_list.rds')




message('Remove samples with unknown tissue\n')
#There are none for this dataset
tissue = unique(metaSMP$tissue)
samples = intersect(unique(c(colnames(EXP_rawcounts[['sanger']]), colnames(EXP_rawcounts[['broad']]))),  
                    metaSMP$model_id[ which(metaSMP$tissue %in% tissue) ])
EXP_rawcounts <- lapply(EXP_rawcounts, function(EXP){EXP[ , colnames(EXP) %in% samples ]})




message('Imputing NA (ZEROmethod) \n')
#Therea are none for this dataset
exp_na <- lapply(EXP_rawcounts, function(EXP){# save gene:sample pairs that are NAs
                index <- which(is.na(EXP), arr.ind=T)
                paste(rownames(EXP)[index[,'row']],  colnames(EXP)[index[,'col']], sep=':')
          })
if ( !is.null(do.call(c,lapply(exp_na,dim))) ){
  EXP_rawcounts <- lapply(EXP_rawcounts, function(EXP){EXP[ which(is.na(EXP)) ] = 0})
}

cat('Total genes and cell lines: \n')
print(lapply(EXP_rawcounts, dim))



#### Exploratory analysis of the raw counts distribution and the sequencing depth ####

# Raw counts distribution
png('xplr/rnaseq_2019-06-26_rawcounts_distribution.png', height = 800, width = 900, res = 150)
par(mfrow = c(1,2))
hist(log2(as.numeric(unlist(EXP_rawcounts[['sanger']]))), breaks = 100, main = 'Cell lines raw counts Sanger', xlab = 'log2(raw counts)')
hist(log2(as.numeric(unlist(EXP_rawcounts[['broad']]))), breaks = 100, main = 'Cell lines raw counts Broad', xlab = 'log2(raw counts)')
dev.off()

# Sequencing Depth
png('xplr/rnaseq_2019-06-26_sequncing_depth.png', height = 800, width = 900, res = 150)
par(mfrow = c(1,2))
barplot(colSums(EXP_rawcounts$sanger), main = 'Cell lines raw counts Sanger', xlab = 'Sum raw counts per Sample')
barplot(colSums(EXP_rawcounts$broad), main = 'Cell lines raw counts Broad', xlab = 'Sum raw counts per Sample')
dev.off()

#### IMPORTANT 
### Because the sequencing depth is different, the datasets won't be merge before the normalitation 
####


rmSamples = F # Turn it T to remove samples that exeed genes with 0

if (rmSamples){
  message('Remove samples exeding genes with 0\n')
  
  zero.count = lapply(EXP_rawcounts,function(EXP){
    apply(EXP, 2, function(x) sum(x==0) ) / nrow(EXP)
  })
  
  print(lapply(zero.count, summary))
  
  zeroThreshold = NULL#Threshold to remove samples. Original script used .4
  
  EXP_rawcounts[['sanger']] = EXP_rawcounts[['sanger']][ , ! colnames(EXP_rawcounts[['sanger']]) %in% names(zero.count[['sanger']])[ zero.count[['sanger']] > zeroThreshold ] ]
  EXP_rawcounts[['broad']] = EXP_rawcounts[['broad']][ , ! colnames(EXP_rawcounts[['broad']]) %in% names(zero.count[['broad']])[ zero.count[['broad']] > zeroThreshold ] ]
}



message('Remove non-expressed genes with average CPM <= 0 ') 
# "The limma-voom method assumes that rows with zero or very low counts have been removed"
exp.zero <- lapply(EXP_rawcounts, function(EXP){
  index = aveLogCPM(EXP) <= 0
  rownames(EXP)[index]
})

saveRDS(exp.zero, file = '1907_release1904/rnaseq_2019-06-26_expZero_list.rds')

EXP_rawcounts <- lapply(EXP_rawcounts, function(EXP){
  keep = aveLogCPM(EXP) > 0
  EXP [ keep, ]
})

cat('Total genes and cell lines: \n')
print(lapply(EXP_rawcounts, dim))

saveRDS(EXP_rawcounts, file = '1907_release1904/rnaseq_2019-06-26_rmZeroGenes_list.rds')




message('Normalization\n')
# Although it is also possible to give a matrix of counts directly to voom without TMM normalization, 
# limma package reommends it
EXP_rawcounts = lapply(EXP_rawcounts, function(EXP){DGEList(counts=EXP,genes=rownames(EXP))})
EXP_normcounts = lapply(EXP_rawcounts, function(EXP){calcNormFactors(EXP, method = 'TMM')})
saveRDS(EXP_normcounts, file = '1907_release1904/rnaseq_2019-06-26_tmm_norm_list.rds')
# EXP_normcounts <- readRDS(file = 'rnaseq_2019-06-26_tmm_norm_list.rds')

#     "method="TMM" is the weighted trimmed mean of M-values (to the reference) proposed by Robinson
#           and Oshlack (2010), where the weights are from the delta method on Binomial data. If refColumn
#           is unspecified, the library whose upper quartile is closest to the mean upper quartile is used.
#     method="RLE" is the scaling factor method proposed by Anders and Huber (2010). We call it
#           "relative log expression", as median library is calculated from the geometric mean of all columns
#           and the median ratio of each sample to the median library is taken as the scale factor.
#     method="upperquartile" is the upper-quartile normalization method of Bullard et al (2010), in
#           which the scale factors are calculated from the 75% quantile of the counts for each library, after
#           removing genes which are zero"




message('Voom transformation\n')
EXP_voom = lapply(EXP_normcounts, function(EXP){voom(EXP, plot = T)$E}) # $E extract the "data corrected for variance"
saveRDS(EXP_voom, file = '1907_release1904/rnaseq_2019-06-26_voom_list.rds')
# EXP_voom <- readRDS(file = 'rnaseq_2019-06-26_voom_list.rds')





message('Building merged matrix \n')
genes = Reduce(intersect, list(rownames(EXP_voom$sanger), rownames(EXP_voom$broad)))

tmp <- lapply(EXP_voom, function(EXP, genes){
  EXP <- as.data.frame(cbind(genes=rownames(EXP), EXP),stringsAsFactors=F)
  EXP[ genes, ]},genes)

cat('Total genes and cell lines: \n')
print(lapply(tmp, dim))

colnames(tmp$sanger) = paste0('S.',colnames(tmp$sanger))
colnames(tmp$broad) = paste0('B.',colnames(tmp$broad))

EXP <- merge(tmp$sanger, tmp$broad, by.x = 'S.genes', by.y = 'B.genes', row)
rownames(EXP) <- EXP$S.genes
EXP <- EXP[,2:ncol(EXP)]
EXP <- data.matrix(EXP, rownames.force = T)

cat('Total genes and cell lines:', dim(EXP), '\n')



#### Exploratory analysis through PCA ####

## variables
catgroup <- as.data.frame(cbind(id=colnames(EXP),sample=gsub('[S/B]\\.','',colnames(EXP)),dataset=substr(colnames(EXP),1,1)),stringsAsFactors = F)
catgroup <- merge(catgroup, metaSMP[,c("model_id","tissue","tissue_status", "cancer_type")], by.x='sample', by.y='model_id')
catgroup <- catgroup[,-1]
# catgroup$suppliers <-  unlist(lapply(strsplit(as.character(catgroup$suppliers), split=':'),function(X){X[1]}))

### sequencing Depth
# sequencingDepth <- do.call(c,lapply(EXP_rawcounts, colSums))
# names(sequencingDepth) <- gsub('sanger','S',names(sequencingDepth))
# names(sequencingDepth) <- gsub('broad','B',names(sequencingDepth))
# 
# catgroup <- merge(catgroup, as.data.frame(cbind(id=names(sequencingDepth),sDepth=sequencingDepth),stringsAsFactors=F), by='id')

## PCA

exp.pca <- prcomp(t(EXP))
str(exp.pca)

exp.pca.variances <- ((exp.pca$sdev^2) / (sum(exp.pca$sdev^2)))*100

png('xplr/rnaseq_2019-06-26_PCA_barplot.png', height = 800, width = 900, res = 150)
barplot(exp.pca.variances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(exp.pca$sdev)), 
        ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,20))
dev.off()

png('xplr/rnaseq_2019-06-26_PCA.png', height = 800, width = 900, res = 150)
qplot(exp.pca$x[,1], exp.pca$x[,2], data=catgroup,
      color= tissue, geom = c('point'), label = dataset,#,'text'
      xlab=paste("PC1, ", round(exp.pca.variances[1], 2), "%"), 
      ylab=paste("PC2, ", round(exp.pca.variances[2], 2), "%"))
dev.off()

##



message('Batch correction of voom transfromed data using ComBat\n')
cat('Using "tisue" and duplicates as covariates\n')
batch = substr(colnames(EXP),1,1)
cov_tissue = as.character(metaSMP$tissue[ match( gsub('[S/B]\\.','',colnames(EXP)), metaSMP$model_id) ])
duplicates = gsub('[S/B]\\.','',colnames(EXP)[duplicated( gsub('[S/B]\\.','',colnames(EXP)) )])
cov_duplicates = sapply(duplicates, function(id) (metaSMP$model_id[ match(gsub('[S/B]\\.','',colnames(EXP)), metaSMP$model_id) ] == id) + 0)
cov_duplicates[ is.na(cov_duplicates) ] = 0
covariates = model.matrix(~cov_tissue+cov_duplicates)
EXP_corvoom = ComBat(EXP, batch = batch, mod = covariates, par.prior = T, prior.plots = T)
saveRDS(EXP_corvoom, file = '1907_release1904/rnaseq_2019-06-26_voom_batchcor.rds')
# EXP_corvoom <- readRDS(file = 'rnaseq_2019-06-26_voom_batchcor.rds')



## PCA

exp.pca <- prcomp(t(EXP_corvoom))
str(exp.pca)

exp.pca.variances <- ((exp.pca$sdev^2) / (sum(exp.pca$sdev^2)))*100

png('xplr/rnaseq_2019-06-26_PCAcomBat_barplot.png', height = 800, width = 900, res = 150)
barplot(exp.pca.variances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(exp.pca$sdev)), 
        ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,20))
dev.off()

png('xplr/rnaseq_2019-06-26_PCAcomBat.png', height = 800, width = 900, res = 150)
qplot(exp.pca$x[,1], exp.pca$x[,2], data=catgroup,
      color= dataset, geom = c('point'), label = dataset,#,'text'
      xlab=paste("PC1, ", round(exp.pca.variances[1], 2), "%"), 
      ylab=paste("PC2, ", round(exp.pca.variances[2], 2), "%"))
dev.off()

##

message('Merging duplicates with the mean')
samples = unique( gsub('[S/B]\\.','',colnames(EXP)) )
genes = rownames(EXP_corvoom)
X = matrix(NA, nrow = length(genes), ncol = length(samples), dimnames = list(genes, samples))
for( s in samples){
  x = EXP_corvoom[ , grep(s, colnames(EXP_corvoom)) ]
  if ( is.matrix(x) )
    x = apply(x, 1, mean)
  X[, s] = x
}
EXPmerged = X
saveRDS(EXPmerged, file = '1907_release1904/rnaseq_2019-06-26_voom_batchcor_merged.rds')
