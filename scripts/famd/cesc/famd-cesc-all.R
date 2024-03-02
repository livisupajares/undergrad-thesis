# ~~~~~ FAMED Analysis for CESC Cohort ~~~~~

# ===== Load libraries =====
library("tidyverse")
library("naniar")
library("missMDA")
library("factoextra")
library("FactoMineR")
library("ggpubr")

# ===== FIXING DATA =====
# Load fixed data set with empty data as NA
cesc <- read.table("database/fixed/fixed-gdc-tcga-cesc.csv",
                   header = TRUE,
                   sep = ",",
                   na.strings = c(NA, ""))

