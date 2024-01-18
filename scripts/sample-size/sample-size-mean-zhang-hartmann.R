# Program to calculate the sample size of Gene Expression Glutaminases.

# Load important packages
library(dplyr)
library(effectsize)
library(pwr)

# Load original dataset
cesc <- read.table('database/original/original-gdc-tcga-cesc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))
hnsc <- read.table('database/original/original-gdc-tcga-hnsc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))
lihc <- read.table('database/original/original-gdc-tcga-lihc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))

# Rename gene expression variable
names(cesc)[names(cesc) == "ENSG00000115419.11"] <- "gls1_gene_expression"
names(cesc)[names(cesc) == "ENSG00000135423.11"] <- "gls2_gene_expression"
names(hnsc)[names(hnsc) == "ENSG00000115419.11"] <- "gls1_gene_expression"
names(hnsc)[names(hnsc) == "ENSG00000135423.11"] <- "gls2_gene_expression"
names(lihc)[names(lihc) == "ENSG00000115419.11"] <- "gls1_gene_expression"
names(lihc)[names(lihc) == "ENSG00000135423.11"] <- "gls2_gene_expression"

# Filter data by tumor vs control
primary_tumor_cesc <- cesc %>% filter(cesc$sample_type.samples == 'Primary Tumor')
solid_tissue_normal_cesc <- cesc %>% filter(cesc$sample_type.samples == 'Solid Tissue Normal')
primary_tumor_hnsc <- hnsc %>% filter(hnsc$sample_type.samples == 'Primary Tumor')
solid_tissue_normal_hnsc <- hnsc %>% filter(hnsc$sample_type.samples == 'Solid Tissue Normal')
primary_tumor_lihc <- lihc %>% filter(lihc$sample_type.samples == 'Primary Tumor')
solid_tissue_normal_lihc <- lihc %>% filter(lihc$sample_type.samples == 'Solid Tissue Normal')

# Remove NAs in gene expression
primary_tumor_cesc <- primary_tumor_cesc[complete.cases(primary_tumor_cesc$gls1_gene_expression), ]
solid_tissue_normal_cesc <- solid_tissue_normal_cesc[complete.cases(solid_tissue_normal_cesc$gls1_gene_expression), ]
primary_tumor_hnsc <- primary_tumor_hnsc[complete.cases(primary_tumor_hnsc$gls1_gene_expression), ]
solid_tissue_normal_hnsc <- solid_tissue_normal_hnsc[complete.cases(solid_tissue_normal_hnsc$gls1_gene_expression), ]
primary_tumor_lihc <- primary_tumor_lihc[complete.cases(primary_tumor_lihc$gls1_gene_expression), ]
solid_tissue_normal_lihc <- solid_tissue_normal_lihc[complete.cases(solid_tissue_normal_lihc$gls1_gene_expression), ]

# Calculate mean
# CESC
# GLS1
mean_gls1_exp_tumor_cesc <- mean(primary_tumor_cesc$gls1_gene_expression)
mean_gls1_exp_control_cesc <- mean(solid_tissue_normal_cesc$gls1_gene_expression)
# GLS2
mean_gls2_exp_tumor_cesc <- mean(primary_tumor_cesc$gls2_gene_expression)
mean_gls2_exp_control_cesc <- mean(solid_tissue_normal_cesc$gls2_gene_expression)
# HNSC
# GLS1
mean_gls1_exp_tumor_hnsc <- mean(primary_tumor_hnsc$gls1_gene_expression)
mean_gls1_exp_control_hnsc <- mean(solid_tissue_normal_hnsc$gls1_gene_expression)
# GLS2
mean_gls2_exp_tumor_hnsc <- mean(primary_tumor_hnsc$gls2_gene_expression)
mean_gls2_exp_control_hnsc <- mean(solid_tissue_normal_hnsc$gls2_gene_expression)
# LIHC
# GLS1
mean_gls1_exp_tumor_lihc <- mean(primary_tumor_lihc$gls1_gene_expression)
mean_gls1_exp_control_lihc <- mean(solid_tissue_normal_lihc$gls1_gene_expression)
# GLS2
mean_gls2_exp_tumor_lihc <- mean(primary_tumor_lihc$gls2_gene_expression)
mean_gls2_exp_control_lihc <- mean(solid_tissue_normal_lihc$gls2_gene_expression)

# Calculate sd
# CESC
# GLS1
sd_pooled_gls1_exp_cesc <- sd_pooled(primary_tumor_cesc$gls1_gene_expression, solid_tissue_normal_cesc$gls1_gene_expression)
# GLS2
sd_pooled_gls2_exp_cesc <- sd_pooled(primary_tumor_cesc$gls2_gene_expression, solid_tissue_normal_cesc$gls2_gene_expression)
# HNSC
# GLS1
sd_pooled_gls1_exp_hnsc <- sd_pooled(primary_tumor_hnsc$gls1_gene_expression, solid_tissue_normal_hnsc$gls1_gene_expression)
# GLS2
sd_pooled_gls2_exp_hnsc <- sd_pooled(primary_tumor_hnsc$gls2_gene_expression, solid_tissue_normal_hnsc$gls2_gene_expression)
# LIHC
# GLS1
sd_pooled_gls1_exp_lihc <- sd_pooled(primary_tumor_lihc$gls1_gene_expression, solid_tissue_normal_lihc$gls1_gene_expression)
# GLS2
sd_pooled_gls2_exp_lihc <- sd_pooled(primary_tumor_lihc$gls2_gene_expression, solid_tissue_normal_lihc$gls2_gene_expression)

# Calculate effect size d, tumor = group 1, control = group 2
# CESC
# GLS1
d_pooled_gls1_exp_cesc <- abs((mean_gls1_exp_tumor_cesc-mean_gls1_exp_control_cesc)/sd_pooled_gls1_exp_cesc)
# GLS2
d_pooled_gls2_exp_cesc <- abs((mean_gls2_exp_tumor_cesc-mean_gls2_exp_control_cesc)/sd_pooled_gls2_exp_cesc)
# HNSC
# GLS1
d_pooled_gls1_exp_hnsc <- abs((mean_gls1_exp_tumor_hnsc-mean_gls1_exp_control_hnsc)/sd_pooled_gls1_exp_hnsc)
# GLS2
d_pooled_gls2_exp_hnsc <- abs((mean_gls2_exp_tumor_hnsc-mean_gls2_exp_control_hnsc)/sd_pooled_gls2_exp_hnsc)
# LIHC
# GLS1
d_pooled_gls1_exp_lihc <- abs((mean_gls1_exp_tumor_lihc-mean_gls1_exp_control_lihc)/sd_pooled_gls1_exp_lihc)
# GLS2
d_pooled_gls2_exp_lihc <- abs((mean_gls2_exp_tumor_lihc-mean_gls2_exp_control_lihc)/sd_pooled_gls2_exp_lihc)

# Calculate sample size
# CESC
# GLS1
pwr.t.test(d=d_pooled_gls1_exp_cesc,power=0.95,sig.level=0.05,type="two.sample",alternative="two.sided")
cat('Total sample size GLS1 in CESC = ',round(118.9294*2, digits = 0), '\n')
# GLS2
pwr.t.test(d=d_pooled_gls2_exp_cesc,power = 0.95,sig.level = 0.05,type="two.sample",alternative="two.sided")
cat('Total sample size GLS2 in CESC = ',round(37.23257*2, digits = 0), '\n')
# HNSC
# GLS1
pwr.t.test(d=d_pooled_gls1_exp_hnsc,power=0.95,sig.level=0.05,type="two.sample",alternative="two.sided")
cat('Total sample size GLS1 in HNSC = ',round(15.89826*2, digits = 0), '\n')
# GLS2
pwr.t.test(d=d_pooled_gls2_exp_hnsc,power = 0.95,sig.level = 0.05,type="two.sample",alternative="two.sided")
cat('Total sample size GLS2 in HNSC = ',round(20450.78*2, digits = 0), '\n')
# LIHC
# GLS1
pwr.t.test(d=d_pooled_gls1_exp_lihc,power=0.95,sig.level=0.05,type="two.sample",alternative="two.sided")
cat('Total sample size GLS1 in LIHC = ',round(26.12672*2, digits = 0), '\n')
# GLS2
pwr.t.test(d=d_pooled_gls2_exp_lihc,power = 0.95,sig.level = 0.05,type="two.sample",alternative="two.sided")
cat('Total sample size GLS2 in LIHC = ',round(15.26547*2, digits = 0), '\n')