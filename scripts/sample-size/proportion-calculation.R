# Load original dataset
cesc <- read.table('database/original/original-gdc-tcga-cesc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))
hnsc <- read.table('database/original/original-gdc-tcga-hnsc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))
lihc <- read.table('database/original/original-gdc-tcga-lihc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))

# Calculate proportions
proportion_cesc <- prop.table(table(cesc$sample_type.samples))
proportion_hnsc <- prop.table(table(hnsc$sample_type.samples))
proportion_lihc <- prop.table(table(lihc$sample_type.samples))

# Print proportions
cat("CESC\n")
print(proportion_cesc)
cat("HNSC\n")
print(proportion_hnsc)
cat("LIHC\n")
print(proportion_lihc)