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
