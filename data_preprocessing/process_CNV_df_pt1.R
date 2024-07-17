#--------------------------------------------#
#   Title: Add gene names to input dataframe
#   Author: Lesley Chapman Hannah
#   Notes:
#   Example Input:
#     arg1: 'AnnotSV dataframe'
#     arg2: 'file name for gene list'
#     arg3: 'Dataframe with select varaints for cases and controls'
#     arg4: 'EUR or AFR'
#-------------------------------------------#

library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

file1 <- args[1]
name1 <- args[2]
file2 <- args[3]
name2 <- args[4]


#Process AnnotSV dataframe
df <- fread(file1)
df3 <- df %>% select(c(SV_chrom,SV_start,SV_end,Annotation_mode,Gene_name)) %>% filter(Annotation_mode == 'full')
df3 <- df3 %>% filter(!grepl(';', Gene_name))
df3 <- df3 %>% mutate(CHROM = paste0('chr', SV_chrom))
df3 <- df3 %>% mutate(ID = paste0(CHROM,'_',SV_start,'_',SV_end,'_',Gene_name))
df4 <- as.data.frame(df3)
df4$count <- nchar(df4$Gene_name)
df4 <- df4 %>% filter(count != 0)
write.table(df4,row.names=FALSE, col.names = T, quote=FALSE, sep = '\t', file=paste0(name1,'.tsv'))

#Process CNV gt dataframe
df <- fread(file2)
names(df) <- sub('^\\[[0-9]+\\]', '', names(df))
names(df) <- sub(':GT$', '', names(df))
df <- df %>%
  mutate_all(funs(str_replace(., "0/0", "0"))) %>%
  mutate_all(funs(str_replace(., "0/1", "1"))) %>%
  mutate_all(funs(str_replace(., "1/0", "1"))) %>%
  mutate_all(funs(str_replace(., "1/1", "2"))) %>%
  mutate_all(funs(str_replace(., "0\\|0", "0"))) %>%
  mutate_all(funs(str_replace(., "0\\|1", "1"))) %>%
  mutate_all(funs(str_replace(., "1\\|0", "1"))) %>%
  mutate_all(funs(str_replace(., "1\\|1", "2"))) %>%
  mutate_all(funs(str_replace(., "\\.", "0"))) %>%
  mutate_all(funs(str_replace(., "\\./\\.", "0"))) %>%
  mutate_all(funs(str_replace(., "./.", "0"))) %>%
rename(START = POS)

#EUR: V1051
#AFR: V1403
df_del_common <- df %>% select(-c(V1403)) %>% rename(CHROM = `# [1]CHROM`) %>% filter(ALT == "<DEL>") %>% filter(POPMAX_AF >= 0.01)
df_ins_common <- df %>% select(-c(V1403)) %>% rename(CHROM = `# [1]CHROM`) %>% filter(ALT == "<INS>") %>% filter(POPMAX_AF >= 0.01)
df_del_rare <- df %>% select(-c(V1403)) %>% rename(CHROM = `# [1]CHROM`) %>% filter(ALT == "<DEL>") %>% filter(POPMAX_AF < 0.01)
df_ins_rare <- df %>% select(-c(V1403)) %>% rename(CHROM = `# [1]CHROM`) %>% filter(ALT == "<INS>") %>% filter(POPMAX_AF < 0.01)
write.table(df_del_common,row.names=FALSE, col.names = T, quote=FALSE, sep = '\t', file=paste0('pCA_del_common_',name2,'_no_gene.tsv'))
write.table(df_ins_common,row.names=FALSE, col.names = T, quote=FALSE, sep = '\t', file=paste0('pCA_ins_common_',name2,'_no_gene.tsv'))
write.table(df_del_rare,row.names=FALSE, col.names = T, quote=FALSE, sep = '\t', file=paste0('pCA_del_rare_',name2,'_no_gene.tsv'))
write.table(df_ins_rare,row.names=FALSE, col.names = T, quote=FALSE, sep = '\t', file=paste0('pCA_ins_rare_',name2,'_no_gene.tsv'))




