# Prostate Cancer (PCa) Analysis in an Ancestrally Diverse Cohort

The following repository lists code to evaluate germline copy number variants (CNVs) associated with prostate cancer risk in men of African and European ancestry using the NIH *All of Us* Cohort. 


## Overview
The scripts include code that run the following:

### Data Preprocessing [./data_preprocessing/]
Data (vcf) download 
```
sh download_cnv_vcf.sh
```

Scripts to add gene names and case control labels.

Suggested preliminary step(s):
- Use [AnnotSV](https://pages.github.com/) to list gene names based on: chromosome, start and end positions
- Generate list of cases and controls using your pre-selected group

```
sh process_CNV_df_pt1.sh
sh process_CNV_df_pt2.sh
```
### Data Analysis [./data_analysis/]
Use Bayesian logistic regression model to estimate CNVs associated with case control output in men of African and European ancestry independently

```
sh run_bayes_model.sh
```

### Data Analysis [./data_transfer/]
Save output dataframe in Google Cloud bucket

```
sh transfer_files_to_bucket.sh
```

