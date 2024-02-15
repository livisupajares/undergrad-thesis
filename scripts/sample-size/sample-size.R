# Load original dataset
cesc <- read.table('database/original/original-gdc-tcga-cesc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))
hnsc <- read.table('database/original/original-gdc-tcga-hnsc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))
lihc <- read.table('database/original/original-gdc-tcga-lihc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))

# Calculate proportions
proportion_cesc <- prop.table(table(cesc$sample_type.samples))
proportion_hnsc <- prop.table(table(hnsc$sample_type.samples))
proportion_lihc <- prop.table(table(lihc$sample_type.samples))

# Get proportion values
prop_cesc_solid_tissue_normal <- proportion_cesc["Solid Tissue Normal"]
prop_hnsc_solid_tissue_normal <- proportion_hnsc["Solid Tissue Normal"]
prop_lihc_solid_tissue_normal <- proportion_lihc["Solid Tissue Normal"]

# Calculate sample size
# Constants
z_score <- 1.96 # Z score for C.I = 0.95 (95%)
margin_of_error <- 0.05 # Margin of error of 5%

cesc_sample_size <- (317*(z_score^2)*prop_cesc_primary_tumor*prop_cesc_solid_tissue_normal)/((margin_of_error^2)*(317-1)+(z_score^2)*prop_cesc_primary_tumor*prop_cesc_solid_tissue_normal)
cat("The final sample size of CESC is: ",cesc_sample_size, "\n")
# p = Solid Tissue Normal 
# q = 1 - p = Primary Tumor

hnsc_sample_size <- (612*(z_score^2)*prop_hnsc_primary_tumor*prop_hnsc_solid_tissue_normal)/((margin_of_error^2)*(612-1)+(z_score^2)*prop_hnsc_primary_tumor*prop_hnsc_solid_tissue_normal)
cat("The final sample size of HNSC is: ",hnsc_sample_size, "\n")

lihc_sample_size <- (469*(z_score^2)*prop_lihc_primary_tumor*prop_lihc_solid_tissue_normal)/((margin_of_error^2)*(469-1)+(z_score^2)*prop_lihc_primary_tumor*prop_lihc_solid_tissue_normal)
cat("The final sample size of LIHC is: ",lihc_sample_size)