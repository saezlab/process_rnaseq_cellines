########################################################################
#
# process_rnaseq_cellines
# process RNAseq data from Cell Lines. From raw data to voom and ComBat
# batch corrected.
#
# Copyright (C) 2019 saezlab, Uniklinikum Heidelberg
# Author: Rosa Hernansaiz-Ballesteros
#
# gdscVScmp.r
# Comparision of old dataset with the new release on CMP
# (https://cellmodelpassports.sanger.ac.uk/)
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


library('data.table')
library('ggplot2')
library('ggforce')
library('openxlsx')

rm(list=ls());cat("\014");if(length(dev.list())>0){dev.off()};options(warn=-1)


###### set paths and variables ######

wd = '~/Projects/process_rnaseq_cellines' #working directory
setwd(wd)

# CMPs
fcmp.rna.processed <- 'CMP/1907_release1904/rnaseq_2019-06-26_voom_batchcor_merged.rds'
fcmp.meta.sample <- '../00-CMP/metadata/model_list_2019-04-04_0916.csv'
fcmp.meta.gene <- '../00-CMP/metadata/gene_identifiers_2019-02-19_1024.csv'

# GDSC_CLLE_GENENTECH
lgcg.rna.processed <- 'GDSC_CLLE_GENENTECH/voom_batchcor_merged.RData'
lgcg.meta.sample <-  '../00-Luz/Cell_Lines_Details.xlsx'


# Loading datasets

cmp <- readRDS(fcmp.rna.processed)
load(lgcg.rna.processed)
gcg <- EXPmerged
rm(EXPmerged)
colnames(gcg) <- gsub('cosmic\\.','',colnames(gcg))

cmp.mSample <- read.delim(fcmp.meta.sample, sep = ',')
cmp.mGenes <- read.delim(fcmp.meta.gene, sep = ',')


#### shared Samples ####

cmp_gcg.samples <- intersect(colnames(gcg), cmp.mSample$COSMIC_ID[ match(colnames(cmp), cmp.mSample$model_id) ])
cmpS <- cmp.mSample$COSMIC_ID[ match(colnames(cmp), cmp.mSample$model_id) ]
cmpS <- as.character(cmpS [ which( !cmpS%in%cmp_gcg.samples ) ])
gcgS <- colnames(gcg)[ which( !colnames(gcg)%in%cmp_gcg.samples ) ]

gg.venn <- data.frame(x = c(-0.866, 0.866),
                      y = c(0, 0),
                      labels = c('CMP', 'GDSC+CLLE+GENENTECH'))



gg.sample <- as.data.frame(matrix(data=NA, nrow=3, ncol = 4, dimnames = list(NULL, c('label', 'x', 'y', 'count'))), stringsAsFactors = F)
gg.sample$label <- c('shared', 'CMP', 'GDSC+CLLE+GENENTECH')
gg.sample$x <- c(0, -1.5, 1.5)
gg.sample$y <- c(0, 0, 0)
gg.sample$count <- c(length(cmp_gcg.samples), length(cmpS), length(gcgS))


png('comparison/samples_voomComBat.png', height = 600, width = 600, res = 150)
ggplot(gg.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
    geom_circle(alpha = .3, size = 1, colour = 'grey') +
    coord_fixed() +
    labs(fill = NULL) +
    theme_void() + theme(legend.position = 'bottom') +
    annotate("text", x = gg.sample$x, y = gg.sample$y, label = gg.sample$count , size = 10) +
    ggtitle('Samples CMPs (2019-04-15_1133) and \n GDSC+CLLE+GENENTECH RNAseq datasets') + 
    scale_fill_manual(values = c('#E69F00', '#009E73',  '#0072B2'))
dev.off()

#### shared Genes ####

cmp_gcg.genes <- intersect(rownames(gcg), cmp.mGenes$ensembl_gene_id[ match(rownames(cmp), cmp.mGenes$gene_id) ])
cmpG <- cmp.mGenes$ensembl_gene_id[ match(rownames(cmp), cmp.mGenes$gene_id) ]
cmpG <- as.character(cmpG [ which( !cmpG%in%cmp_gcg.genes ) ])
gcgG <- colnames(gcg)[ which( !rownames(gcg)%in%cmp_gcg.genes ) ]

gg.gene <- as.data.frame(matrix(data=NA, nrow=3, ncol = 4, dimnames = list(NULL, c('label', 'x', 'y', 'count'))), stringsAsFactors = F)
gg.gene$label <- c('shared', 'CMP', 'GDSC+CLLE+GENENTECH')
gg.gene$x <- c(0, -1.5, 1.5)
gg.gene$y <- c(0, 0, 0)
gg.gene$count <- c(length(cmp_gcg.genes), length(cmpG), length(gcgG))


png('comparison/genes_voomComBat.png', height = 600, width = 600, res = 150)
ggplot(gg.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .3, size = 1, colour = 'grey') +
  coord_fixed() +
  labs(fill = NULL) +
  theme_void() + theme(legend.position = 'bottom') +
  annotate("text", x = gg.gene$x, y = gg.gene$y, label = gg.gene$count , size = 10) +
  ggtitle('Genes CMPs (2019-04-15_1133) and \n GDSC+CLLE+GENENTECH RNAseq datasets') + 
  scale_fill_manual(values = c('#E69F00', '#009E73',  '#0072B2'))
dev.off()


# distribution of gene expression per dataset on the shared genes (14933) on the shared samples (1047)

#sum of expressions by gene
cmp.tmp <- rowSums( cmp[ as.character(cmp.mGenes$gene_id[ match(cmp_gcg.genes, cmp.mGenes$ensembl_gene_id) ]), 
                as.character(cmp.mSample$model_id[ match(cmp_gcg.samples, cmp.mSample$COSMIC_ID) ]) ] )
names(cmp.tmp) <- cmp.mGenes$ensembl_gene_id [ match(names(cmp.tmp), cmp.mGenes$gene_id) ]

cmp.tmp <- cmp.tmp [order(factor(names(cmp.tmp)))]

gcg.tmp <- rowSums( gcg[ cmp_gcg.genes, cmp_gcg.samples ] ) 
gcg.tmp <- gcg.tmp [order(factor(names(gcg.tmp)))]

# melt and merge by gene
gg.dat <- merge( cbind(gene = names(gcg.tmp), melt(gcg.tmp,value.name='gcg')), 
                 cbind(gene = names(cmp.tmp), melt(cmp.tmp,value.name='cmp')), by = 'gene')
gg.dat <- melt(gg.dat, variable.name = 'dataset', value.name = 'sumExpr', id.var = 'gene')

png('comparison/density_expression_voomComBat.png', height = , width = 800, res = 150)
ggplot(data = gg.dat, aes(x = sumExpr, color = dataset)) +
  geom_density() + 
  scale_color_manual(values = c('#009E73','#E69F00'),
                     labels = c("GDSC+CCLE+GENENTECH", "CMP"))
dev.off()






















###### begin of the comparison ######

cat('Load expression data \n')
cmpRAW <- fread(rnadata, sep = ',')

cmpSamples <- unique(cmpRAW[,c('model_id','dataset_name')])

load(gdscDATA)
gdscVOOM <- EXPmerged # HGNC vs Cosmic
rm(EXPmerged)

gdscR7 <- fread(gdscR7file)

cat('Load metainformation \n')

samples <- read.delim(metasample, sep = ',')
genes <- read.delim(metagene, sep = ',')

gdscmeta <- read.xlsx(gdscMETA)

cat('Exploration of Samples and Genes across CMP \n')

cmpSamples <- merge(cmpSamples, samples[,c('model_id','model_name','COSMIC_ID')], by='model_id')

cat('Sample-based venn diagrams \n')
samplesCMP <- setDF(unique.data.frame(cmpRAW[,c('model_id','dataset_name')]))

## CMP different datasets
sharedSamples <- intersect(samplesCMP$model_id[which(samplesCMP$dataset_name=='Sanger RNASeq')], samplesCMP$model_id[which(samplesCMP$dataset_name=='Broad (CCLE) RNASeq')])

gg.venn <- data.frame(x = c(-0.866, 0.866),
                      y = c(0, 0),
                      labels = c('Sanger', 'Broad'))

gg.venndat <- as.data.frame(matrix(data=NA, nrow=3, ncol = 4, dimnames = list(NULL, c('label', 'x', 'y', 'count'))), stringsAsFactors = F)
gg.venndat$label <- c('shared', 'Sanger', 'Broad')
gg.venndat$x <- c(0, -1.5, 1.5)
gg.venndat$y <- c(0, 0, 0)
gg.venndat$count <- c(length(sharedSamples), 
                      sum(samplesCMP$dataset_name=='Sanger RNASeq')-length(sharedSamples),
                      sum(samplesCMP$dataset_name=='Broad (CCLE) RNASeq')-length(sharedSamples))

(venn.Samples <- ggplot(gg.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .3, size = 1, colour = 'grey') +
  coord_fixed() +
  labs(fill = NULL) +
  theme_void() + theme(legend.position = 'bottom') +
  annotate("text", x = gg.venndat$x, y = gg.venndat$y, label = gg.venndat$count , size = 12) +
  ggtitle('Samples in CMPs RNAseq dataset 2019-04-15_1133'))


ggsave(filename = '1907_process_rnaseq_CMP/venn_samples_BS.svg', plot = venn.Samples) #w=450, h=400

# Comparison CMP and GDSC (processed batch)

samplesCMP <-  merge(samplesCMP, samples[,c('model_id','model_type','COSMIC_ID')], by='model_id')

G <- gsub('cosmic\\.','', colnames(gdscVOOM))
S <- as.character(samplesCMP$COSMIC_ID[which(samplesCMP$dataset_name=='Sanger RNASeq')])
# there is a double cosmic Id in S, "1331031;1331030", which map to the same model_id , SIDM00400. 
# it is replace by 1331030, as it the onle mapping into the models of the old GDSC dataset
S[which(S=="1331031;1331030")] <- 1331030
B <- as.character(samplesCMP$COSMIC_ID[which(samplesCMP$dataset_name=='Broad (CCLE) RNASeq')])

SGB <- intersect(intersect(S,B),G)
BS <- intersect(S,B)
BS <- BS[!BS%in%SGB]
BG <- intersect(G,B)
BG <- BG[!BG%in%SGB]
SG <- intersect(S,G)
SG <- SG[!SG%in%SGB]
G <- G[!(G%in%SGB | G%in%BG | G%in%SG)]
S <- S[! (S%in%SGB | S%in%SG | S%in%BS )]
B <- B[! (B%in%SGB | B%in%BG | B%in%BS )]

gg.venndat <- as.data.frame(matrix(data=NA, nrow=7, ncol = 4, dimnames = list(NULL, c('label', 'x', 'y', 'count'))), stringsAsFactors = F)
gg.venndat$label <- c('GDSC', 'Sanger', 'Broad', 'Broad-Sanger', 'Broad-GDSC', 'Sanger-GDSC', 'All')
gg.venndat$x <- c(0, -1.2, 1.2, 0, 0.8, -0.8, 0)
gg.venndat$y <- c(1.2, -0.6, -0.6, -1, 0.5, 0.5, 0)
gg.venndat$count <- c(length(G), length(S), length(B),
                      length(BS), length(BG), length(SG), length(SGB))

gg.venn3G <- data.frame(x = c(0, 0.866, -0.866),
                        y = c(1, -0.5, -0.5),
                        labels = c('GDSC + CCLE + Genentech (ours)', 'CMP-Broad', 'CMP-Sanger'))


(venn.SBG <- ggplot(gg.venn3G, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .3, size = 1, colour = 'grey') +
  coord_fixed() +
  labs(fill = NULL) +
  theme_void() + theme(legend.position = 'bottom') +
  annotate("text", x = gg.venndat$x, y = gg.venndat$y, label = gg.venndat$count , size = 8) +
  scale_fill_manual(values = c('#E69F00', '#009E73',  '#0072B2')) +
  ggtitle('Samples in CMP and GDSC RNAseq datasets'))

ggsave(filename = '1907_process_rnaseq_CMP/venn_samples_GBS.svg', plot = venn.Samples) #w=500, h=565

#check samples not found in GDSC
S
# [1] "907284"  "1303910"

samples[samples$COSMIC_ID=='905967',] #K052 #Sarc9371

gdscmeta[gdscmeta$Sample.Name=='HUH6',]

B
# "905967" 'NCI-H322M'
head(samples[which(samples$COSMIC_ID=='' & samples$model_id%in%samplesCMP$model_id[which(samplesCMP$dataset_name=='Broad (CCLE) RNASeq')]),c('CCLE_ID','model_name','synonyms')])

missingCosmic <- rbind(cmpSamples[which(cmpSamples$COSMIC_ID%in%c("907284", "1303910", "905967")),],
                       cmpSamples[which(cmpSamples$COSMIC_ID==''),])
write.table(missingCosmic,file = '~/Projects/00-CMP/expression/mismatching_samples.tsv',quote = F, sep = '\t', row.names = F)

#comparison GDSC R7 with processed GDSC

samplesR7 <- gsub('DATA\\.','',colnames(gdscR7)[3:ncol(gdscR7)])
samplesVOOM <- gsub('cosmic\\.','',colnames(gdscVOOM))

gg.venndat <- as.data.frame(matrix(data=NA, nrow=3, ncol = 4, dimnames = list(NULL, c('label', 'x', 'y', 'count'))), stringsAsFactors = F)
gg.venndat$label <- c('shared', 'release7', 'merged')
gg.venndat$x <- c(0, -1.5, 1.5)
gg.venndat$y <- c(0, 0, 0)
gg.venndat$count <- c(length(intersect(samplesVOOM,samplesR7)), 
                      length(samplesR7[which(!samplesR7%in%intersect(samplesVOOM,samplesR7))]),
                      length(samplesVOOM[which(!samplesVOOM%in%intersect(samplesVOOM,samplesR7))]))


(venn.Samples <- ggplot(gg.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
    geom_circle(alpha = .3, size = 1, colour = 'grey') +
    coord_fixed() +
    labs(fill = NULL) +
    theme_void() + theme(legend.position = 'bottom') +
    annotate("text", x = gg.venndat$x, y = gg.venndat$y, label = gg.venndat$count , size = 12) +
    ggtitle('Samples GDSC release 7 and merged Luz dataset'))


rm(samplesCMP, sharedSamples, venn.Samples, venn.SBG, G,S,B,SGB,BS,BG,SB)

cat('Gene-based venn diagrams \n')
##venn diagram of the genes distributed across the Sanger and the Broad institute.
geneCMP <- setDF(unique.data.frame(cmpRAW[,c('gene_id','dataset_name')]))

sharedGenes <- intersect(geneCMP$gene_id[which(geneCMP$dataset_name=='Sanger RNASeq')], geneCMP$gene_id[which(geneCMP$dataset_name=='Broad (CCLE) RNASeq')])

gg.venndat <- as.data.frame(matrix(data=NA, nrow=3, ncol = 4, dimnames = list(NULL, c('label', 'x', 'y', 'count'))), stringsAsFactors = F)
gg.venndat$label <- c('shared', 'Sanger', 'Broad')
gg.venndat$x <- c(0, -1.5, 1.5)
gg.venndat$y <- c(0, 0, 0)
gg.venndat$count <- c(length(sharedGenes), 
                      sum(geneCMP$dataset_name=='Sanger RNASeq')-length(sharedGenes),
                      sum(geneCMP$dataset_name=='Broad (CCLE) RNASeq')-length(sharedGenes))

venn.Genes <- ggplot(gg.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .3, size = 1, colour = 'grey') +
  coord_fixed() +
  labs(fill = NULL) +
  theme_void() + theme(legend.position = 'bottom') +
  annotate("text", x = gg.venndat$x, y = gg.venndat$y, label = gg.venndat$count , size = 10) +
  ggtitle('Genes in CMPs RNAseq dataset 2019-04-15_1133')


ggsave(filename = '1907_process_rnaseq_CMP/venn_genes_BS.svg', plot = venn.Samples) #w=450, h=400


#Exploration of Samples and Genes across CMP and GDSC

geneCMP <- merge(geneCMP, genes, by='gene_id')
S <- as.character(geneCMP$hgnc_symbol[which(geneCMP$dataset_name=='Sanger RNASeq')])
B <- as.character(geneCMP$hgnc_symbol[which(geneCMP$dataset_name=='Broad (CCLE) RNASeq')])
SGB <- intersect(intersect(S,B),rownames(gdscVOOM))
BS <- intersect(S,B)
BS <- BS[!BS%in%SGB]
BG <- intersect(as.character(geneCMP$hgnc_symbol[which(geneCMP$dataset_name=='Broad (CCLE) RNASeq')]),rownames(gdscVOOM))
BG <- BG[!BG%in%SGB]
SG <- intersect(as.character(geneCMP$hgnc_symbol[which(geneCMP$dataset_name=='Sanger RNASeq')]),rownames(gdscVOOM))
SG <- SG[!SG%in%SGB]
G <- rownames(gdscVOOM)[!(rownames(gdscVOOM)%in%SGB | rownames(gdscVOOM)%in%BG | rownames(gdscVOOM)%in%SG)]
S <- S[! (S%in%SGB | S%in%SG | S%in%BS )]
B <- B[! (B%in%SGB | B%in%BG | B%in%BS )]

ggplot(gg.venn3G, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .3, size = 1, colour = 'grey') +
  coord_fixed() +
  theme_void()  + theme(legend.position = 'bottom')


gg.venndat <- as.data.frame(matrix(data=NA, nrow=7, ncol = 4, dimnames = list(NULL, c('label', 'x', 'y', 'count'))), stringsAsFactors = F)
gg.venndat$label <- c('GDSC', 'Sanger', 'Broad', 'Broad-Sanger', 'Broad-GDSC', 'Sanger-GDSC', 'All')
gg.venndat$x <- c(0, -1.2, 1.2, 0, 0.8, -0.8, 0)
gg.venndat$y <- c(1.2, -0.6, -0.6, -1, 0.5, 0.5, 0)
gg.venndat$count <- c(length(G), length(S), length(B),
                      length(BS), length(BG), length(SG), length(SGB))

venn.SBG <- ggplot(gg.venn3G, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .3, size = 1, colour = 'grey') +
  coord_fixed() +
  labs(fill = NULL) +
  theme_void() + theme(legend.position = 'bottom') +
  annotate("text", x = gg.venndat$x, y = gg.venndat$y, label = gg.venndat$count , size = 8) +
  scale_fill_manual(values = c('#E69F00', '#009E73',  '#0072B2')) +
  ggtitle('Genes in CMP and GDSC RNAseq datasets') 

ggsave(filename = '1907_process_rnaseq_CMP/venn_genes_GBS.svg', plot = venn.Samples) #w=500, h=565

rm(geneCMP, sharedGenes, gg.venn, gg.venndat, gg.venn3G, venn.Genes, venn.SBG, G,S,B,SGB,BS,BG,SB)

