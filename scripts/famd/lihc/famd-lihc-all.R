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
lihc <- read.table("database/fixed/fixed-gdc-tcga-lihc.csv",
                   header = TRUE,
                   sep = ",",
                   na.strings = c(NA, ""))

# Change col names for easier visualization
col_names_lihc <- c("M",
                    "pHep",
                    "tHep",
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
                    "genero",
                    "IMC",
                    "cIMC",
                    "TM",
                    "GH",
                    "EC",
                    "OS")

colnames(lihc) <- col_names_lihc

