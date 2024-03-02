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

