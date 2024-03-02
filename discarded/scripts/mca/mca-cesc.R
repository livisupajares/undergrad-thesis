# ~~~~~ MCA Analysis for CESC Cohort ~~~~~

# ===== Load libraries =====
library("ggplot2")
library("tidyverse")
library("naniar")
library("missMDA")
library("FactoMineR")
library("factoextra")
library("haven")
library("haven")
library("corrplot")

# ===== FIXING DATA =====
# Load fixed dataset with empty data as NA
cesc <- read.table("database/fixed/fixed-gdc-tcga-cesc.csv", header = TRUE, sep = ",", na.strings = c(NA, ""))

# Change col names for easier visualization
col.names.cesc <- c("M", "pHPV", "tHPV", "tHPV2", "TM2", "TM.PV", "sncGLS1", "encGLS1", "2encGLS1", "eGLS1", "sncGLS2", "encGLS2", "2encGLS2", "eGLS2", "edad", "etnicidad", "raza", "IMC", "cIMC", "TM", "GH", "EC", "genero", "OS", "CT")

colnames(cesc) <- col.names.cesc 

# Subset necessary variables
cesc.active <- subset(cesc, select = c("pHPV", "tHPV", "TM.PV", "TM", "GH", "EC", "CT"))

# Covert from character to factor
cesc.active <- cesc.active |>
  mutate(across(where(is.character), as.factor))

# ===== FIX MISSING VALUES =====
# See a table of number of missing values and percentages by variable
miss_var_summary(cesc.active)

# Estimate ncp to predict the number of components used to predict the missing entries
estim_ncpMCA(cesc.active)

# Use imputeMCA to impute (fill in) missing values in the dataset before applying MCA.
active.impute.mca <- imputeMCA(cesc.active, ncp = 1)

# ===== MCA =====
# Set the maximun number of overlaps for plots using ggrepel globally
# Default is 10
options(ggrepel.max.overlaps = 80)

# Create the MCA model
cesc.mca <- MCA(cesc.active, tab.disj = active.impute.mca$tab.disj, graph = FALSE)
print(cesc.mca)

# ===== Inspect Principal Components or Dimensions =====

# Extract the eigenvalues/variances retained by each PC (axis).
eig.val <- get_eigenvalue(cesc.mca)
head(eig.val)

# Visualize the eigenvalues as a line
fviz_eig(cesc.mca,  
         choice = 'eigenvalue', 
         geom = 'line')

# Visualize PC as percentage of explained variances
fviz_screeplot(cesc.mca, addlabels = TRUE)

# ===== Biplot (Graph of individuals & variables) =====
fviz_mca_biplot(cesc.mca, # MCA results
                repel = TRUE, # avoid superposition of the graph's labels
                alpha.ind = 0.1, # transparency of rows or  samples
                alpha.var = 1, # transparency of rows or variables
                ggtheme = theme_minimal())

# ===== Describe Principal Components =====
res.desc <- dimdesc(cesc.mca,
                    axes = c(1,2) #PCA 1 y 2
)
res.desc[[1]]
res.desc[[2]]

# ===== Variables =====
var <- get_mca_var(cesc.mca)
print(var)

# coordinates 
head(var$coord)

# cos2: quality on the factor map
head(var$cos2)

# Contribution in factors or dimensions 
head(var$contrib)

# correlation of factors with PC 1 and 2
fviz_mca_var(cesc.mca, choice = "mca.cor", ggtheme = theme_minimal(), repel = TRUE)

# Best variable contribution
fviz_mca_var(cesc.mca, col.var = "contrib",
             gradient.cols = c("#00AFBB", 
                               "#E7B800", 
                               "#FC4E07"),
             axes = c(1,2),
             repel = TRUE)

# Contribution to PC 1 & 2
fviz_contrib(cesc.mca, choice = "var", axes = c(1,2))