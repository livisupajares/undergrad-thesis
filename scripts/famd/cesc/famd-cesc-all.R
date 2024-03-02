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

# ===== FAMD =====
# Set the maximun number of overlaps for plots using ggrepel globally
# Default is 10
options(ggrepel.max.overlaps = 80)

# Create the FAMD model
cesc_famd <- FAMD(cesc_active,
                  ncp = 10,
                  sup.var = 4:10,
                  tab.disj = active_impute_famd$tab.disj,
                  graph = FALSE)
print(cesc_famd)

# ===== Inspect Principal Components or Dimensions =====

# Extract the eigenvalues/variances retained by each PC (axis).
eig_val <- get_eigenvalue(cesc_famd)
print(eig_val)

# Visualize the eigenvalues as a line
line_eig <- fviz_eig(cesc_famd,
                     choice = "eigenvalue",
                     title = "Eigenvalores de los PCs (CESC)",
                     xlab = "Componentes Principales",
                     ylab = "Eigenvalores",
                     ylim = c(0, 3),
                     geom = c("bar", "line"),
                     barfill = "#95C4CB",
                     barcolor = "#258B9A",
                     linecolor = "#5A7081",
                     addlabels = TRUE,
                     ggtheme = theme_gray())
print(line_eig)

# Visualize PC as percentage of explained variances
screeplot <- fviz_screeplot(cesc_famd,
                            title = "Scree Plot (CESC)",
                            xlab = "Componentes principales",
                            ylab = "Porcentaje de varianzas explicadas",
                            barfill = "#95C4CB",
                            barcolor = "#258B9A",
                            ggtheme = theme_gray(),
                            addlabels = TRUE)
print(screeplot)

