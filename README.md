# process_rnaseq_cellines
Processing RNAseq data from Cell Lines. The process_rnaseq_cellines_* files
are pipelines that follow the next steps:

 - Imputes 0 to NA values
 - Removes samples exeding genes with 0
 - Removes non-expressed genes with average CPM <=0
 - Normalises through TMM method
 - Runs Voom transformation
 - Uses ComBat for Batch correction of voom transformed data 
 - Merges duplicates by the mean

In the wiki of this repo the exploratory analysis of the data can be found.

# Technical details

R version 3.6.0 (2019-04-26)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.5

Packages:

- data.table 1.12.2
- dplyr 0.8.1
- stringr 1.4.0
- limma 3.40.2
- sva 3.32.1
- edgeR 3.26.4
- ggplot2 3.1.1
- ggforce 0.2.2
- ggfortify 0.4.7

