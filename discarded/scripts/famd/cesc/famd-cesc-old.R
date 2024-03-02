# Load libraries
library("ggplot2")
library("dplyr")
library("naniar")
library("missMDA")
library("FactoMineR")


# Load fixed dataset with empty data as NA
cesc <- read.table("database/fixed/fixed-gdc-tcga-cesc.csv", header = TRUE, sep = ",", na.strings = c(NA, ""))

# Fix cesc dataset
## Remove clinical variables, copy number status and segment and reorder columns
## Only leave GLS1 expression, hpv status, hpv type, infection status sample, sample type 2, histologic grade, clinical stage and tobacco smoking history

cesc.active <- subset(cesc, select = -c(sample, sample.type, hpv.type.2, GLS1.copy.number.status.2, GLS2.copy.number.status.2, age.at.initial.pathologic.diagnosis, ethnicity, race, bmi, bmi.category, gender, GLS1.copy.number.status, GLS2.copy.number.status))

## Move sample.type.2 to the front
cesc.active <- cesc.active %>% relocate(sample.type.2)

## Convert characters to factors, see: https://community.rstudio.com/t/estim-ncpfdma-error-message-no-defined-colums/159382/3
cesc.active <- cesc.active |>
  mutate(across(where(is.character), as.factor))

## See the structure of the data
str(cesc.active)

# Fix missing values

## See a table of number of missing values and percentages by variable
library("naniar")
miss_var_summary(cesc)

## Use imputeFAMD to impute (fill in) missing values in the dataset before applying FAMD.
library("missMDA") # use library
impute.famd <- imputeFAMD(cesc.active, ncp = 5)

### Visualize the new dataset in a matrix of 12 by 12
impute.famd$completeObs[1:12, 1:12]

# FAMD
## Set the maximun number of overlaps for plots using ggrepel globally
options(ggrepel.max.overlaps = 18)

## Compute FAMD with the transformed dataset
library("FactoMineR") ## install package
cesc.famd <- FAMD(cesc.active, ncp = 5, tab.disj = impute.famd$tab.disj, graph = FALSE) ## graph of quantitative variables
print(cesc.famd)

## Visualization and interpretation
### Use factoextra functions
library("factoextra") # install library
### Extract the eigenvalues/variances retained by each dimension (axis).
eig.val <- get_eigenvalue(cesc.famd)
print(eig.val)

#### Visualize the eigenvalues/variances. It an be used to draw the scree plot (the percentages of inertia explained by each FAMD dimensions):
fviz_screeplot(cesc.famd, addlabels = TRUE)

### All variables
### used to extract the results for variables. By default, this function returns a list containing the coordinates, the cos2 and the contribution of all variables

## Exlude genomic variables
genomic <- c("GLS1.Copy.Number.Segment", "GLS1.expression.fpkm", "GLS2.Copy.Number.Segment", "GLS2.expression.fpkm", "OS.time.days")

cesc.famd.clinical <- cesc.active[, !names(cesc.active) %in% genomic]

var <- get_famd_var(cesc.famd)
print(var)

#### Access the different components
##### Coordinates of variables
head(var$coord)
##### Cos2: quality of representation on the factore map
head(var$cos2)
##### Contributions to the  dimensions
head(var$contrib)

#### Plot of variables
##### correlation between variables - both quantitative and qualitative variables - and the principal dimensions, as well as, the contribution of variables to the dimensions 1 and 2.

###### fviz_famd_var() to plot both quantitative and qualitative variables.
fviz_famd_var(cesc.famd, repel = TRUE)
###### Contribution to the first dimension
fviz_contrib(cesc.famd, choice = "var", axes = 1, select.vars = cesc.famd.clinical)
###### Contribution to the second dimension
fviz_contrib(cesc.famd, "var", axes = 2) ###### The red dashed line on the graph above indicates the expected average value, If the contributions were uniform.

#### Quatitative variables
##### Extract quantitative variables
quanti.var <- get_famd_var(cesc.famd, "quanti.var")
print(quanti.var)
##### Correlation circle by percentage of contribution to each component
fviz_famd_var(cesc.famd, "quanti.var", repel = TRUE, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

##### Correlation circle by percentage of cos2 values representing the quality of representation on the factor map. If a variable is well represented by two dimensions, the sum of the cos2 is closed to one. For some of the items, more than 2 dimensions might be required to perfectly represent the data.
fviz_famd_var(cesc_famd, "quanti.var", col.var = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

#### Qualitative variables
##### Extract qualitative variables
quali.var <- get_famd_var(cesc_famd, "quali.var")
print(quali.var)
##### Visualize qualitative variables categories
fviz_famd_var(cesc.famd, "quali.var", col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

##### Gráfico de individuos (muestras --> índices)
ind <- get_famd_ind(cesc.famd)
print(ind)

# Plot individuals by cos2 values
# Individuals with similar profiles are close to each other on the factor map.
fviz_famd_ind(cesc.famd, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

# it’s possible to color the individuals using any of the qualitative variables in the initial data table.
fviz_mfa_ind(cesc.famd, habillage = "infection.status.sample", palette = c("#00AFBB", "#E7B800", "#FC4E07"), addEllipses = TRUE, ellipse.type = "confidence", repel = TRUE)

# color individuals using multiple categorical variables at the same time
fviz_ellipses(cesc_famd, c("sample.type.2", "hpv.status"), repel = TRUE)
