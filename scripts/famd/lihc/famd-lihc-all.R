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

# Subset genotipic and phenotipyc variables
lihc_active <- subset(lihc,
                      select = c("pHep",
                                 "tHep",
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
                                 "OS"))

# Switch OS with eGLS2
lihc_active <- lihc_active %>%
  relocate(OS, .after = eGLS2)

# Covert from character to factor
lihc_active <- lihc_active |>
  mutate(across(where(is.character), as.factor))

# ===== FIX MISSING VALUES =====
# See a table of number of missing values and percentages by variable
miss_var_summary(lihc_active)

# Estimate ncp to predict the number of components used to predict the missing
# entries
estim_ncpFAMD(lihc_active, sup.var = 5:10)

# Use imputeFAMD to impute (fill in) missing
# values in the dataset before applying FAMD.
# sup.var are GLS1/GLS2 gene expression,
# GLS1/GLS2 CNV and OS
active_impute_famd <- imputeFAMD(lihc_active, sup.var = 4:10, ncp = 4)

# ===== FAMD =====
# Set the maximun number of overlaps for plots using ggrepel globally
# Default is 10
options(ggrepel.max.overlaps = 80)

# Create the FAMD model
lihc_famd <- FAMD(lihc_active,
                  ncp = 10,
                  sup.var = 4:10,
                  tab.disj = active_impute_famd$tab.disj,
                  graph = FALSE)
print(lihc_famd)

# ===== Inspect Principal Components or Dimensions =====

# Extract the eigenvalues/variances retained by each PC (axis).
eig_val <- get_eigenvalue(lihc_famd)
print(eig_val)

# Visualize the eigenvalues as a line
line_eig <- fviz_eig(lihc_famd,
                     choice = "eigenvalue",
                     title = "Eigenvalores de los PCs (LIHC)",
                     xlab = "Componentes Principales",
                     ylab = "Eigenvalores",
                     ylim = c(0, 2.5),
                     geom = c("bar", "line"),
                     barfill = "#95C4CB",
                     barcolor = "#258B9A",
                     linecolor = "#5A7081",
                     addlabels = TRUE,
                     ggtheme = theme_gray())
print(line_eig)

# Visualize PC as percentage of explained variances
screeplot <- fviz_screeplot(lihc_famd,
                            title = "Scree Plot (LIHC)",
                            xlab = "Componentes principales",
                            ylab = "Porcentaje de varianzas explicadas",
                            barfill = "#95C4CB",
                            barcolor = "#258B9A",
                            ggtheme = theme_gray(),
                            addlabels = TRUE)
print(screeplot)

# ===== All variables =====
# used to extract the results for variables. By default, this function returns a
# list containing the coordinates, the cos2 and the contribution
# of all variables
var <- get_famd_var(lihc_famd)
print(var)

# Coordinates of variables
head(var$coord)
# Cos2: quality of representation on the factor map
head(var$cos2)
# Contributions to the  dimensions
head(var$contrib)

# Contribution to the first dimension
contrib_1 <- fviz_contrib(lihc_famd,
                          fill = "#95C4CB",
                          color = "#258B9A",
                          choice = "var",
                          ggtheme = theme_gray(),
                          title = "Contribución de variables al 1er PC",
                          axes = 1)

contrib_1_pc <- ggpar(contrib_1, ylab = "Contribuciones (%)")
print(contrib_1_pc)

# Contribution to the second dimension
contrib_2 <- fviz_contrib(lihc_famd,
                          title = "Contribución de variables al 2do PC",
                          fill = "#95C4CB",
                          color = "#258B9A",
                          choice = "var",
                          ggtheme = theme_gray(),
                          axes = 2)

contrib_2_pc <- ggpar(contrib_2, ylab = "Contribuciones (%)")
print(contrib_2_pc)

# Contribution to both dim
contrib_gen <- fviz_contrib(lihc_famd,
                            title = "Contribución de variables a los PC 1 y 2",
                            fill = "#95C4CB",
                            color = "#258B9A",
                            choice = "var",
                            ggtheme = theme_gray(),
                            axes = c(1, 2))

contrib_gen_pc <- ggpar(contrib_gen, ylab = "Contribuciones (%)")
print(contrib_gen_pc)

# Best variable contribution
# Green : Supplementary variables
contrib_var <- fviz_famd_var(lihc_famd,
                             choice = "var",
                             col.var = "contrib",
                             ggtheme = theme_gray(),
                             axes = c(1, 2),
                             gradient.cols = c("#288B9A",
                                               "#D9AF39",
                                               "#E85B63"),
                             pointsize = 2,
                             show.legend.text = TRUE,
                             title = "Contribución de variables a los PC 1 y 2 - FAMD",
                             repel = TRUE)
print(contrib_var)