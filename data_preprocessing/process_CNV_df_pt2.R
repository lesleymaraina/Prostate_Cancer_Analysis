#--------------------------------------------#
#   Title: Add dataframe case control labels
#   Author: Lesley Chapman Hannah
#   Notes:
#   Example Input:
#     arg1: 'pCA_del_common_AFR_no_gene.tsv'
#     arg2: 'pCA_del_common_EUR_no_gene.tsv'
#     arg3: 'annotSV_pca_gene_list_AFR.tsv'
#     arg4: 'annotSV_pca_gene_list_EUR.tsv'
#     arg5: 'df_AFR_ca_co_labels.tsv'
#     arg6: 'df_EUR_ca_co_labels.tsv'
#-------------------------------------------#

#Import Libraries
library(data.table)
library(tidyverse)
library(rstanarm)

dir <- getwd()
setwd(dir)

args <- commandArgs(trailingOnly = TRUE)
pCA_AFR_gt_df <- args[1]
pCA_EUR_gt_df <- args[2]
pCA_AFR_annotsv <- args[3]
pCA_EUR_annotsv <- args[4]
pCA_AFR_caco_labels <- args[5]
pCA_EUR_caco_labels <- args[6]

#Import dataframes
#genotype dataframes
pCA_AFR <- fread(pCA_AFR_gt_df)
pCA_EUR <- fread(pCA_EUR_gt_df)
#annotSV gene names dataframes
pCA_AFR_gene <- fread(pCA_AFR_annotsv)
pCA_EUR_gene <- fread(pCA_EUR_annotsv)
#case/control labels dataframes
pCA_AFR_labels <- fread(pCA_AFR_caco_labels)
pCA_EUR_labels <- fread(pCA_EUR_caco_labels)

# Process samples of African Ancestry
#------------------------------------#
pCA_AFR_gene_v2 <- pCA_AFR_gene %>% select(-c(SV_chrom, Annotation_mode, Gene_name, count)) %>% rename(START = SV_start, END = SV_end) %>% mutate(START = START - 1)

#Add gene names to GT dataframe
pCA_AFR_gene <- pCA_AFR %>% inner_join(pCA_AFR_gene_v2, by=c('CHROM', 'START', 'END')) %>% select(-c(CHROM,START,REF,ALT,END,POPMAX_AF))
pCA_AFR_gene_gt <- as.data.frame(t(pCA_AFR_gene))
colnames(pCA_AFR_gene_gt) <- pCA_AFR_gene$ID
pCA_AFR_gene_gt$ID <- rownames(pCA_AFR_gene_gt)
pCA_AFR_gene_gt <- pCA_AFR_gene_gt %>% mutate(across(c(ID), as.integer)) %>% inner_join(pCA_AFR_labels, by='ID')
pCA_AFR_gene_gt$Label_var <- ifelse(grepl("case",pCA_AFR_gene_gt$Label),'1','0')
pCA_AFR_gene_gt <- pCA_AFR_gene_gt %>% select(-c(ID, Label))

pCA_AFR_gene_gt2 <- pCA_AFR_gene_gt %>% mutate_if(is.character, as.integer)
#Drop columns with all zero values
pCA_AFR_gene_gt2 <- pCA_AFR_gene_gt2[colSums(abs(pCA_AFR_gene_gt2), na.rm = TRUE) > 0]
print(head(pCA_AFR_gene_gt2)

write.table(pCA_AFR_gene_gt2, row.names=FALSE, col.names = T, quote=FALSE, sep = '\t', file='pCA_del_common_AFR_input_df.tsv')

# Process samples of European Ancestry
#------------------------------------#
pCA_EUR_gene_v2 <- pCA_EUR_gene %>% select(-c(SV_chrom, Annotation_mode, Gene_name, count)) %>% rename(START = SV_start, END = SV_end) %>% mutate(START = START - 1)

#Add gene names to GT dataframe
pCA_EUR_gene <- pCA_EUR %>% inner_join(pCA_EUR_gene_v2, by=c('CHROM', 'START', 'END')) %>% select(-c(CHROM,START,REF,ALT,END,POPMAX_AF))
pCA_EUR_gene_gt <- as.data.frame(t(pCA_EUR_gene))
colnames(pCA_EUR_gene_gt) <- pCA_EUR_gene$ID
pCA_EUR_gene_gt$ID <- rownames(pCA_EUR_gene_gt)
pCA_EUR_gene_gt <- pCA_EUR_gene_gt %>% mutate(across(c(ID), as.integer)) %>% inner_join(pCA_EUR_labels, by='ID')
pCA_EUR_gene_gt$Label_var <- ifelse(grepl("case",pCA_EUR_gene_gt$Label),'1','0')
pCA_EUR_gene_gt <- pCA_EUR_gene_gt %>% select(-c(ID, Label))

pCA_EUR_gene_gt2 <- pCA_EUR_gene_gt %>% mutate_if(is.character, as.integer)
#Drop columns with all zero values
pCA_EUR_gene_gt2 <- pCA_EUR_gene_gt2[colSums(abs(pCA_EUR_gene_gt2), na.rm = TRUE) > 0]

write.table(pCA_EUR_gene_gt2, row.names=FALSE, col.names = T, quote=FALSE, sep = '\t', file='pCA_del_common_EUR_input_df.tsv')




