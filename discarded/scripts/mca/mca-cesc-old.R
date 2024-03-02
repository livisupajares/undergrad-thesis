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

active.var <- subset(cesc, select = -c(sample, sample.type, hpv.type.2, GLS1.copy.number.status.2, GLS2.copy.number.status.2, age.at.initial.pathologic.diagnosis, ethnicity, race, bmi, bmi.category, gender, GLS1.copy.number.status, GLS2.copy.number.status, GLS1.Copy.Number.Segment, GLS2.Copy.Number.Segment, GLS1.expression.fpkm, GLS2.expression.fpkm, OS.time.days))

## Supplementary variables
quanti.sup <- subset(cesc, select = c(GLS1.Copy.Number.Segment,GLS2.Copy.Number.Segment,GLS1.expression.fpkm,GLS2.expression.fpkm,OS.time.days))
quali.sup <- subset(cesc, select = c(GLS1.copy.number.status,GLS2.copy.number.status,ethnicity,race,bmi.category,gender))

## Move sample.type.2 to the front
active.var <- active.var %>% relocate(sample.type.2)

## Change Row names
# Change the row names
col.names.active <- c("TM", "PHPV", "THPV","TM.PV","GH","EC","CT")
colnames(active.var) <- col.names.active
col.names.quanti <- c("SNCgls1", "SNCgls2", "Egls1", "Egls2", "OS")
colnames(quanti.sup) <- col.names.quanti
col.names.cuali <- c("ENCgls1", "ENCgls2", "Etnicidad", "Raza", "IMC", "Género")
colnames(quali.sup) <- col.names.cuali

## Preview your data
head(active.var)
head(quanti.sup)
head(quali.sup)

## Convert characters to factors, see: https://community.rstudio.com/t/estim-ncpfdma-error-message-no-defined-colums/159382/3
active.var <- active.var |>
  mutate(across(where(is.character), as.factor))
quali.sup <- quali.sup |>
  mutate(across(where(is.character), as.factor))

## Data summary of frequency of factors
summary(active.var)

# Fix missing values

## See a table of number of missing values and percentages by variable
library("naniar")
miss_var_summary(active.var)

## Use imputeFAMD to impute (fill in) missing values in the dataset before applying FAMD.
library("missMDA") # use library
impute.mca.active <- imputeMCA(active.var, 
                               ncp = 5)

### Visualize the new dataset in a matrix of 12 by 12
impute.mca.active$completeObs[1:7,1:7]

# MCA
## Set the maximun number of overlaps for plots using ggrepel globally
options(ggrepel.max.overlaps = 10)

## Compute MCA with the transformed dataset
library("FactoMineR") ## install package
cesc.mca <- MCA(active.var, 
                ncp = 5, 
                tab.disj = impute.mca.active$tab.disj,
                graph = FALSE) ## graph of quantitative variables
print(cesc.mca)

## Visualization and interpretation
### Use factoextra functions
library("factoextra") # install library
### Extract the eigenvalues/variances retained by each dimension (axis).
eig.val <- get_eigenvalue(cesc.mca)
print(eig.val)

#### Visualize the eigenvalues/variances. It an be used to draw the scree plot (the percentages of inertia explained by each FAMD dimensions):
fviz_screeplot(cesc.mca, 
               addlabels = TRUE)

## Biplot
fviz_mca_biplot(cesc.mca,
                repel = TRUE, # Avoid text overlapping (slow if many point
                ggtheme = theme_minimal())

### All variables
### used to extract the results for variables. By default, this function returns a list containing the coordinates, the cos2 and the contribution of all variables
var <- get_mca_var(cesc.mca)
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
## correlation
fviz_mca_var(cesc.mca, choice = "mca.cor", ggtheme = theme_minimal(), repel = TRUE)
## variable categories
fviz_mca_var(cesc.mca, ggtheme = theme_minimal(), repel = TRUE)
###### Contribution to the first dimension
fviz_contrib(cesc.mca, "var", axes = 1)
###### Contribution to the second dimension
fviz_contrib(cesc.mca, "var", axes = 2) ###### The red dashed line on the graph above indicates the expected average value, If the contributions were uniform.

#### Quatitative variables
##### Extract quantitative variables
quanti.var <- get_mca_var(cesc_famd, "quanti.var")
print(quanti.var)
##### Correlation circle by percentage of contribution to each component
fviz_famd_var(cesc_famd, "quanti.var", repel = TRUE, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

##### Correlation circle by percentage of cos2 values representing the quality of representation on the factor map. If a variable is well represented by two dimensions, the sum of the cos2 is closed to one. For some of the items, more than 2 dimensions might be required to perfectly represent the data.
fviz_famd_var(cesc_famd, "quanti.var", col.var = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

#### Qualitative variables
##### Extract qualitative variables
quali.var <- get_famd_var(cesc_famd, "quali.var")
print(quali.var)
##### Visualize qualitative variables categories
fviz_famd_var(cesc_famd, "quali.var", col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

##### Gráfico de individuos (muestras --> índices)
ind <- get_famd_ind(cesc_famd)
print(ind)

# Plot individuals by cos2 values
# Individuals with similar profiles are close to each other on the factor map. 
fviz_famd_ind(cesc_famd, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

# it’s possible to color the individuals using any of the qualitative variables in the initial data table. 
fviz_mfa_ind(cesc_famd, habillage = "infection.status.sample", palette = c("#00AFBB", "#E7B800", "#FC4E07"), addEllipses = TRUE, ellipse.type = "confidence", repel = TRUE) 

# color individuals using multiple categorical variables at the same time
fviz_ellipses(cesc_famd, c("sample.type.2", "hpv.status"), repel = TRUE)