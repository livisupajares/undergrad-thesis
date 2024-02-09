# Load necessry libraries
library(dplyr)

# Load original dataset
cesc <- read.table('database/original/original-gdc-tcga-cesc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))
hnsc <- read.table('database/original/original-gdc-tcga-hnsc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))
lihc <- read.table('database/original/original-gdc-tcga-lihc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))

# Fix category names
cesc <- cesc %>%
  mutate(across(where(is.character), ~ ifelse(. == "Primary Tumor", "primary_tumor", .)))
cesc <- cesc %>%
  mutate(across(where(is.character), ~ ifelse(. == "Solid Tissue Normal", "solid_tissue_normal", .)))
hnsc <- hnsc %>%
  mutate(across(where(is.character), ~ ifelse(. == "Primary Tumor", "primary_tumor", .)))
hnsc <- hnsc %>%
  mutate(across(where(is.character), ~ ifelse(. == "Solid Tissue Normal", "solid_tissue_normal", .)))
lihc <- lihc %>%
  mutate(across(where(is.character), ~ ifelse(. == "Primary Tumor", "primary_tumor", .)))
lihc <- lihc %>%
  mutate(across(where(is.character), ~ ifelse(. == "Solid Tissue Normal", "solid_tissue_normal", .)))

# Remove other values that aren't primary tumor or solid tissue normal
cesc <- subset(cesc, cesc$sample_type.samples != "Metastatic")
hnsc <- subset(hnsc, hnsc$sample_type.samples != "Metastatic")
lihc <- subset(lihc, lihc$sample_type.samples != "Metastatic")
lihc <- subset(lihc, lihc$sample_type.samples != "Recurrent Tumor")

# Calculate proportions
proportion_cesc <- prop.table(table(cesc$sample_type.samples))
proportion_hnsc <- prop.table(table(hnsc$sample_type.samples))
proportion_lihc <- prop.table(table(lihc$sample_type.samples))

# Get proportion values
for (sample_type.samples in names(proportion_cesc)) {
  assign(paste0("prop_cesc_", sample_type.samples), proportion_cesc[sample_type.samples])
}

for (sample_type.samples in names(proportion_hnsc)) {
  assign(paste0("prop_hnsc_", sample_type.samples), proportion_hnsc[sample_type.samples])
}

for (sample_type.samples in names(proportion_lihc)) {
  assign(paste0("prop_lihc_", sample_type.samples), proportion_lihc[sample_type.samples])
}

# Constants
z_score <- 1.96 # Z score for C.I = 0.95 (95%)
margin_of_error <- 0.05 # Margin of error of 5%

# Calculate sample size
cesc_sample_size <- (317*(z_score^2)*prop_cesc_primary_tumor*prop_cesc_solid_tissue_normal)/((margin_of_error^2)*(317-1)+(z_score^2)*prop_cesc_primary_tumor*prop_cesc_solid_tissue_normal)
cat("The final sample size of CESC is: ",cesc_sample_size, "\n")

hnsc_sample_size <- (612*(z_score^2)*prop_hnsc_primary_tumor*prop_hnsc_solid_tissue_normal)/((margin_of_error^2)*(612-1)+(z_score^2)*prop_hnsc_primary_tumor*prop_hnsc_solid_tissue_normal)
cat("The final sample size of HNSC is: ",hnsc_sample_size, "\n")

lihc_sample_size <- (469*(z_score^2)*prop_lihc_primary_tumor*prop_lihc_solid_tissue_normal)/((margin_of_error^2)*(469-1)+(z_score^2)*prop_lihc_primary_tumor*prop_lihc_solid_tissue_normal)
cat("The final sample size of LIHC is: ",lihc_sample_size)