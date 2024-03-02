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

# Change col names for easier visualization
col_names_cesc <- c("M",
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
                    "IMC",
                    "cIMC",
                    "TM",
                    "GH",
                    "EC",
                    "genero",
                    "OS",
                    "CT")

colnames(cesc) <- col_names_cesc

# Subset genotipic and phenotipyc variables
cesc_active <- subset(cesc,
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
                                 "CT"))

# Switch OS with eGLS2
cesc_active <- cesc_active %>%
  relocate(OS, .after = eGLS2)

# Covert from character to factor
cesc_active <- cesc_active |>
  mutate(across(where(is.character), as.factor))

# ===== FIX MISSING VALUES =====
# See a table of number of missing values and percentages by variable
miss_var_summary(cesc_active)

# Estimate ncp to predict the number of components used to predict the missing
# entries
estim_ncpFAMD(cesc_active, sup.var = 5:10)

# Use imputeFAMD to impute (fill in) missing
# values in the dataset before applying FAMD.
# sup.var are GLS1/GLS2 gene expression,
# GLS1/GLS2 CNV and OS
active_impute_famd <- imputeFAMD(cesc_active, sup.var = 4:10, ncp = 1)

