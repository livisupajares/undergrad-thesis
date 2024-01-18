# Load important packages
library(dplyr)
library(pwr)

# Load original dataset
cesc <- read.table('database/original/original-gdc-tcga-cesc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))
hnsc <- read.table('database/original/original-gdc-tcga-hnsc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))
lihc <- read.table('database/original/original-gdc-tcga-lihc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))

# Filter data by tumor vs control
primary_tumor_cesc <- cesc %>% filter(cesc$sample_type.samples == 'Primary Tumor')
solid_tissue_normal_cesc <- cesc %>% filter(cesc$sample_type.samples == 'Solid Tissue Normal')
primary_tumor_hnsc <- hnsc %>% filter(hnsc$sample_type.samples == 'Primary Tumor')
solid_tissue_normal_hnsc <- hnsc %>% filter(hnsc$sample_type.samples == 'Solid Tissue Normal')
primary_tumor_lihc <- lihc %>% filter(lihc$sample_type.samples == 'Primary Tumor')
solid_tissue_normal_lihc <- lihc %>% filter(lihc$sample_type.samples == 'Solid Tissue Normal')

# Calculate amount of total observations
total_cesc <- nrow(cesc)
total_hnsc <- nrow(hnsc)
total_lihc <- nrow(lihc)

# Calculate amount of observations for tumor
obs_primary_tumor_cesc <- nrow(primary_tumor_cesc)
obs_primary_tumor_hnsc <- nrow(primary_tumor_hnsc)
obs_primary_tumor_lihc <- nrow(primary_tumor_lihc)

# Calculate amount of observations for control
obs_solid_tissue_normal_cesc <- nrow(solid_tissue_normal_cesc)
obs_solid_tissue_normal_hnsc <- nrow(solid_tissue_normal_hnsc)
obs_solid_tissue_normal_lihc <- nrow(solid_tissue_normal_lihc)

# Calculate proportions for tumor (p1)
p1_primary_tumor_cesc <- obs_primary_tumor_cesc/total_cesc
p1_primary_tumor_hnsc <- obs_primary_tumor_hnsc/total_hnsc
p1_primary_tumor_lihc <- obs_primary_tumor_lihc/total_lihc

# Calculate proportions for normal (p2)
p2_solid_tissue_normal_cesc <- obs_solid_tissue_normal_cesc/total_cesc
p2_solid_tissue_normal_hnsc <- obs_solid_tissue_normal_hnsc/total_hnsc
p2_solid_tissue_normal_lihc <- obs_solid_tissue_normal_lihc/total_lihc

# Calculate sample size
# CESC
power.prop.test(p1 = p1_primary_tumor_cesc, p2 = p2_solid_tissue_normal_cesc, sig.level=0.05, power=0.95)
cat('Total sample size in CESC = ',round(3.530042*2, digits = 0), '\n')
# HNSC
power.prop.test(p1 = p1_primary_tumor_hnsc, p2 = p2_solid_tissue_normal_hnsc, sig.level=0.05, power=0.95)
cat('Total sample size in HNSC = ',round(8.967603*2, digits = 0), '\n')
# LIHC
power.prop.test(p1 = p1_primary_tumor_lihc, p2 = p2_solid_tissue_normal_lihc, sig.level=0.05, power=0.95)
cat('Total sample size in LIHC = ',round(14.0752*2, digits = 0), '\n')