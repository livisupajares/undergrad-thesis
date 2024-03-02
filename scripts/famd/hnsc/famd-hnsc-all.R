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

# Change col names for easier visualization
col_names_hnsc <- c("M",
                    "pHPV",
                    "tHPV",
                    "tHPV2",
                    "TM2",
                    "TM.PV",
                    "sncGLS1",
                    "encGLS1",
                    "2encGLS1",
                    "eGLS1",
                    "sncGLS2",
                    "encGLS2",
                    "2encGLS2",
                    "eGLS2",
                    "edad",
                    "etnicidad",
                    "raza",
                    "TM",
                    "GH",
                    "EC",
                    "OS",
                    "genero",
                    "CA",
                    "CT")

colnames(hnsc) <- col_names_hnsc

# Subset genotipic and phenotipyc variables
hnsc_active <- subset(hnsc,
                      select = c("pHPV",
                                 "tHPV",
                                 "TM.PV",
                                 "sncGLS1",
                                 "encGLS1",
                                 "eGLS1",
                                 "sncGLS2",
                                 "encGLS2",
                                 "eGLS2",
                                 "TM",
                                 "GH",
                                 "EC",
                                 "OS",
                                 "CA",
                                 "CT"))

# Switch OS with eGLS2
hnsc_active <- hnsc_active %>%
  relocate(OS, .after = eGLS2)

# Covert from character to factor
hnsc_active <- hnsc_active |>
  mutate(across(where(is.character), as.factor))

# ===== FIX MISSING VALUES =====
# See a table of number of missing values and percentages by variable
miss_var_summary(hnsc_active)

# Estimate ncp to predict the number of components used to predict the missing
# entries
estim_ncpFAMD(hnsc_active, sup.var = 5:10)

# Use imputeFAMD to impute (fill in) missing
# values in the dataset before applying FAMD.
# sup.var are GLS1/GLS2 gene expression,
# GLS1/GLS2 CNV and OS
active_impute_famd <- imputeFAMD(hnsc_active, sup.var = 4:10, ncp = 3)

