# Load libraries
library("ggplot2")
library("dplyr")

# Load fixed dataset with empty data as NA
cesc <- read.table("database/fixed/fixed-gdc-tcga-cesc.csv", header = TRUE, sep = ",", na.strings = c(NA, "")) 

# Fix cesc dataset
## Remove clinical variables, copy number status and segment and reorder columns
## Only leave GLS1 expression, hpv status, hpv type, infection status sample, sample type 2, histologic grade, clinical stage and tobacco smoking history

cesc_gls1 <- subset(cesc, select = -c(sample, sample.type, hpv.type.2, GLS1.copy.number.status.2, GLS2.copy.number.status.2, age.at.initial.pathologic.diagnosis, ethnicity, race, bmi, bmi.category, gender,OS.time.days, GLS1.Copy.Number.Segment, GLS1.copy.number.status, GLS2.Copy.Number.Segment, GLS2.copy.number.status, GLS2.expression.fpkm))

## Move sample.type.2 to the front
cesc_gls1 <- cesc_gls1 %>% relocate(sample.type.2)

## Convert characters to factors, see: https://community.rstudio.com/t/estim-ncpfdma-error-message-no-defined-colums/159382/3
cesc_gls1 <- cesc_gls1 |>
  mutate(across(where(is.character), as.factor))

## See the structure of the data
str(cesc_gls1)

# Fix missing values

## See a table of number of missing values and percentages by variable
library("naniar")
miss_var_summary(cesc_gls1)

## Use estim_ncpFAMD to estimate the number of dimensions (components or factors) in Factor Analysis of Mixed Data (FAMD).
library("missMDA") # use library
ncp <- estim_ncpFAMD(cesc_gls1)

## Use imputeFAMD to impute (fill in) missing values in the dataset before applying FAMD.
impute_famd <- imputeFAMD(cesc_gls1, ncp = 2) # result: 1

### Visualize the new dataset in a matrix of 8 by 8
impute_famd$completeObs[1:8,1:8]

# FAMD
## Set the maximun number of overlaps for plots using ggrepel globally
options(ggrepel.max.overlaps = 20)

## Compute FAMD with the transformed dataset
library("FactoMineR") ## install package
cesc_gls1_famd <- FAMD(cesc_gls1, ncp = 2, tab.disj = impute_famd$tab.disj) ## graph of quantitative variables
print(cesc_gls1_famd)

## Visualization and interpretation
### Use factoextra functions
library("factoextra") # install library
### Extract the eigenvalues/variances retained by each dimension (axis).
eig_val <- get_eigenvalue(cesc_gls1_famd)
print(eig_val)

### Visualize the eigenvalues/variances.
fviz_screeplot(cesc_gls1_famd)
var <- get_famd_var(cesc_famd2_gls1_gen_exp)

# Coordinates of variables
head(var$coord)
# Cos2: quality of representation on the factore map
head(var$cos2)
# Contributions to the  dimensions
head(var$contrib)

# Plot of variables
fviz_famd_var(cesc_famd2_gls1_gen_exp, repel = TRUE)
# Contribution to the first dimension
fviz_contrib(cesc_famd2_gls1_gen_exp, "var", axes = 1)
# Contribution to the second dimension
fviz_contrib(cesc_famd2_gls1_gen_exp, "var", axes = 2)

quanti.var <- get_famd_var(cesc_famd2_gls1_gen_exp, "quanti.var")

quali.var <- get_famd_var(cesc_famd2_gls1_gen_exp, "quali.var")

fviz_famd_var(cesc_famd2_gls1_gen_exp, "quali.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

ind <- get_famd_ind(cesc_famd2_gls1_gen_exp)

fviz_famd_ind(cesc_famd2_gls1_gen_exp, col.ind = "cos2", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)

fviz_mfa_ind(cesc_famd2_gls1_gen_exp, 
             habillage = "infection.status.sample", # color by groups 
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, ellipse.type = "confidence", 
             repel = TRUE # Avoid text overlapping
) 

fviz_ellipses(cesc_famd2_gls1_gen_exp, c("sample.type.2", "hpv.status"), repel = TRUE)