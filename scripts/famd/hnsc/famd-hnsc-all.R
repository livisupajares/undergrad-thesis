# ~~~~~ FAMED Analysis for HNSC Cohort ~~~~~

# ===== Load libraries =====
library("tidyverse")
library("naniar")
library("missMDA")
library("factoextra")
library("FactoMineR")
library("ggpubr")

# ===== FIXING DATA =====
# Load fixed data set with empty data as NA
hnsc <- read.table("database/fixed/fixed-gdc-tcga-hnsc.csv",
                   header = TRUE,
                   sep = ",",
                   na.strings = c(NA, ""))

