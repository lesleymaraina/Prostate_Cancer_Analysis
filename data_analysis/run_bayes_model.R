#----------------------------------------------------------#
#   Title: Find CNVs associated with case control status
#   Author: Lesley Chapman Hannah
#   Notes: Script to run bayesian model
#   Example Input:
#     arg1: 'Genotype dataframe with case/control labels'
#     arg2: 'chromosome number'
#     arg3: 'Ancestry: EUR or AFR'
#---------------------------------------------------------#

library(data.table)
library(tidyverse)
library(rstanarm)


args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
chr_num <- args[2]
anc <- args[3]

file <- fread(input_file)

file_t <- file %>% select(contains("chr"))
file_t <- as.data.frame(t(file_t))
file_t <- file_t %>% mutate(sumVar = rowSums(select(., contains("V")))) %>% mutate(ratio = sumVar/nrow(xx)) %>% filter(ratio <= 0.01) %>% select(-c(sumVar,ratio))


file2 <- as.data.frame(t(file_t))
file2$Label_var <- file$Label_var
file <- file2

model_1 <- stan_glm(Label_var~., 
                    data=file, 
                    prior = normal(),
                    family = binomial,
                    prior_intercept = normal(0, 1),
                    chains = 4, iter = 5000*2, seed = 84735,
                    prior_PD = FALSE)
summary(model_1)


model_output_df <- as.data.frame(summary(model_1))
model_output_df$varID <- rownames(model_output_df)

write.table(model_output_df,row.names=FALSE, col.names = T, quote=FALSE, sep = '\t', file=paste0('chr',chr_num,'_',anc,'_bayes_model_output.tsv'))


df_post_80 <- as.data.frame(round(posterior_interval(model_1, prob = 0.8), 2))
df_post_95 <- as.data.frame(round(posterior_interval(model_1, prob = 0.95), 2))


df_post_80$varID <- rownames(df_post_80)
df_post_95$varID <- rownames(df_post_95)

write.table(df_post_80,row.names=FALSE, col.names = T, quote=FALSE, sep = '\t', file=paste0('chr',chr_num,'_',anc,'_bayes_80_CI_model_output.tsv'))

write.table(df_post_95,row.names=FALSE, col.names = T, quote=FALSE, sep = '\t', file=paste0('chr',chr_num,'_',anc,'_bayes_95_CI_model_output.tsv'))






