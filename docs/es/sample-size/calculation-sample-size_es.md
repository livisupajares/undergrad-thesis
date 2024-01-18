# Cálculo del tamaño muestral
Livisú Pajares Rojas
01-02-2025

# Cálculo del tamaño muestral

## Fórmula de Murray Larry

Para ello utilizamos la calculadora de
[esta](https://www.calculator.net/sample-size-calculator.html?type=1&cl=95&ci=5&pp=50&ps=317&x=Calculate)
página web.

Para hacer el mismo cálculo en R se realizó lo siguiente en el script
[sample-size-murray-larry.R](https://github.com/livisupajares/undergrad-thesis/blob/master/scripts/sample-size/sample-size-murray-larry.R):

Se calculará el tamaño poblacional para una población finita según el
método de *Murray y Larry* (1) .

$$n=\frac{Z^2\sigma^2N}{e^2(N-1)+Z^2\sigma^2}$$

En donde,

$n$ = Es el tamaño de la población que se desea obtener.

$N$ = Es el tamaño total de la población.

$\sigma$ = Desviación estándar de la población. Si no se conoce, usar
0.5 como constante.

$Z$ = Es el valor obtenido mediante los intervalos de confianza. Un
intervalo de confianza del 95% tiene un valor de $Z$ de 1.96.

$e$ = Es el límite aceptable de error de muestreo.

Primero indicamos el tamaño de la población ($N$):

Para CESC colocamos $N=317$, para HNSC $N=612$, y para LIHC $N=469$.
También están el resto de parámetros de la ecuación de arriba.

``` r
population_size_cesc <- 317
population_size_hnsc <- 612
population_size_lihc <- 469
z_score <- 1.96 # Z score for C.I = 0.95 (95%)
population_sd <- 0.5 # Standard Deviation of the population. Use 50% if not sure
margin_of_error <- 0.05 # Margin of error of 5%
```

Ahora calculamos el tamaño muestral:

``` r
# Calculate sample size
sample_size_cesc <- (z_score^2 * population_sd^2 * population_size_cesc) / (((margin_of_error^2) * (population_size_cesc - 1)) + (z_score^2 * population_sd^2))

sample_size_hnsc <- (z_score^2 * population_sd^2 * population_size_hnsc) / (((margin_of_error^2) * (population_size_hnsc - 1)) + (z_score^2 * population_sd^2))

sample_size_lihc <- (z_score^2 * population_sd^2 * population_size_lihc) / (((margin_of_error^2) * (population_size_lihc - 1)) + (z_score^2 * population_sd^2))
```

Finalmente, el tamaño muestral es:

    The sample size in CESC is  173.9298

    The sample size in HNSC is  236.2494

    The sample size in LIHC is  211.4287

## Z-test: Proporciones - Diferencia entre dos proporciones diferentes

Se realizó en el siguiente script
[sample-size-proportions-zhang-hartmann.R](https://github.com/livisupajares/undergrad-thesis/blob/master/scripts/sample-size/sample-size-proportions-zhang-hartmann.R).
Este método sigue el paper de Zhang y Hartmann (2) . Vamos a necesitar
dos paquetes:

``` r
# Load important packages
library(dplyr)
library(pwr)
```

Ahora subiremos los datasets originales:

``` r
# Load original dataset
cesc <- read.table('database/original/original-gdc-tcga-cesc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))
hnsc <- read.table('database/original/original-gdc-tcga-hnsc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))
lihc <- read.table('database/original/original-gdc-tcga-lihc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))
```

Ahora, con el paquete `dplyr` filtraremos la data para tener dos grupos,
tumores (grupo 1) vs controles (grupo 2) en cada cohorte.

``` r
# Filter data by tumor vs control
primary_tumor_cesc <- cesc %>% filter(cesc$sample_type.samples == 'Primary Tumor') # Grupo 1
solid_tissue_normal_cesc <- cesc %>% filter(cesc$sample_type.samples == 'Solid Tissue Normal') # Grupo 2

primary_tumor_hnsc <- hnsc %>% filter(hnsc$sample_type.samples == 'Primary Tumor') # Grupo 1
solid_tissue_normal_hnsc <- hnsc %>% filter(hnsc$sample_type.samples == 'Solid Tissue Normal') # Grupo 2

primary_tumor_lihc <- lihc %>% filter(lihc$sample_type.samples == 'Primary Tumor') # Grupo 1
solid_tissue_normal_lihc <- lihc %>% filter(lihc$sample_type.samples == 'Solid Tissue Normal') # Grupo 2
```

Para calcular las proporciones, tenemos que saber cuantos tejidos
tenemos en cada dataset:

``` r
# Calculate amount of total observations
total_cesc <- nrow(cesc)
total_hnsc <- nrow(hnsc)
total_lihc <- nrow(lihc)

# Calculate amount of observations for tumor
obs_primary_tumor_cesc <- nrow(primary_tumor_cesc)
obs_primary_tumor_hnsc <- nrow(primary_tumor_hnsc)
obs_primary_tumor_lihc <- nrow(primary_tumor_lihc)

# Calculate amount of observations for control
obs_solid_tissue_normal_cesc <- nrow(solid_tissue_normal_cesc)
obs_solid_tissue_normal_hnsc <- nrow(solid_tissue_normal_hnsc)
obs_solid_tissue_normal_lihc <- nrow(solid_tissue_normal_lihc)
```

Calculamos las proporciones:

``` r
# Calculate proportions for tumor (p1)
p1_primary_tumor_cesc <- obs_primary_tumor_cesc/total_cesc
p1_primary_tumor_hnsc <- obs_primary_tumor_hnsc/total_hnsc
p1_primary_tumor_lihc <- obs_primary_tumor_lihc/total_lihc

# Calculate proportions for normal (p2)
p2_solid_tissue_normal_cesc <- obs_solid_tissue_normal_cesc/total_cesc
p2_solid_tissue_normal_hnsc <- obs_solid_tissue_normal_hnsc/total_hnsc
p2_solid_tissue_normal_lihc <- obs_solid_tissue_normal_lihc/total_lihc
```

Ahora con el paquete `pwr`, calculamos el tamaño muestral:

``` r
# Calculate sample size
# CESC
power.prop.test(p1 = p1_primary_tumor_cesc, p2 = p2_solid_tissue_normal_cesc, sig.level=0.05, power=0.95)

# HNSC
power.prop.test(p1 = p1_primary_tumor_hnsc, p2 = p2_solid_tissue_normal_hnsc, sig.level=0.05, power=0.95)

# LIHC
power.prop.test(p1 = p1_primary_tumor_lihc, p2 = p2_solid_tissue_normal_lihc, sig.level=0.05, power=0.95)
```

Finalmente el tamaño muestral será:

    Total sample size in CESC =  7

    Total sample size in HNSC =  18

    Total sample size in LIHC =  28

## T-test : Diferencia entre dos medias independientes (2 grupos)

Se usó el siguiente script
[sample-size-mean-zhang-hartmann.R](https://github.com/livisupajares/undergrad-thesis/blob/master/scripts/sample-size/sample-size-mean-zhang-hartmann.R).
Este método sigue el paper de Zhang y Hartmann (2) . Vamos a necesitar
tres paquetes en R:

``` r
library(dplyr)
library(effectsize)
library(pwr)
```

Primero importamos la data original:

``` r
cesc <- read.table('database/original/original-gdc-tcga-cesc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))

hnsc <- read.table('database/original/original-gdc-tcga-hnsc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))

lihc <- read.table('database/original/original-gdc-tcga-lihc.tsv', header = TRUE, sep = "\t",na.strings=c(NA,''))
```

Cambiaremos de nombre a la variable expresión de GLS1 y GLS2 para que ya
no salga codificado.

``` r
# Rename gene expression variable
names(cesc)[names(cesc) == "ENSG00000115419.11"] <- "gls1_gene_expression"
names(cesc)[names(cesc) == "ENSG00000135423.11"] <- "gls2_gene_expression"

names(hnsc)[names(hnsc) == "ENSG00000115419.11"] <- "gls1_gene_expression"
names(hnsc)[names(hnsc) == "ENSG00000135423.11"] <- "gls2_gene_expression"

names(lihc)[names(lihc) == "ENSG00000115419.11"] <- "gls1_gene_expression"
names(lihc)[names(lihc) == "ENSG00000135423.11"] <- "gls2_gene_expression"
```

Ahora usaremos `dplyr` para filtrar la data para que tengamos dos
dataframes, uno con los datos que pertenecen al tejido tumoral (Grupo 1)
y el segundo al tejido control (Grupo 2).

``` r
# Filter data by tumor vs control
primary_tumor_cesc <- cesc %>% filter(cesc$sample_type.samples == 'Primary Tumor')
solid_tissue_normal_cesc <- cesc %>% filter(cesc$sample_type.samples == 'Solid Tissue Normal')

primary_tumor_hnsc <- hnsc %>% filter(hnsc$sample_type.samples == 'Primary Tumor')
solid_tissue_normal_hnsc <- hnsc %>% filter(hnsc$sample_type.samples == 'Solid Tissue Normal')

primary_tumor_lihc <- lihc %>% filter(lihc$sample_type.samples == 'Primary Tumor')
solid_tissue_normal_lihc <- lihc %>% filter(lihc$sample_type.samples == 'Solid Tissue Normal')
```

Finalmente, removeremos las filas que tengan NA en las variables de
expresión:

``` r
# Remove NAs in gene expression
primary_tumor_cesc <- primary_tumor_cesc[complete.cases(primary_tumor_cesc$gls1_gene_expression), ]
solid_tissue_normal_cesc <- solid_tissue_normal_cesc[complete.cases(solid_tissue_normal_cesc$gls1_gene_expression), ]

primary_tumor_hnsc <- primary_tumor_hnsc[complete.cases(primary_tumor_hnsc$gls1_gene_expression), ]
solid_tissue_normal_hnsc <- solid_tissue_normal_hnsc[complete.cases(solid_tissue_normal_hnsc$gls1_gene_expression), ]

primary_tumor_lihc <- primary_tumor_lihc[complete.cases(primary_tumor_lihc$gls1_gene_expression), ]
solid_tissue_normal_lihc <- solid_tissue_normal_lihc[complete.cases(solid_tissue_normal_lihc$gls1_gene_expression), ]
```

Un t-test de dos muestras se usa para determinar si las medias de dos
poblaciones, $\mu_{1}$ y $\mu_{2}$ son igules (3) .

$$H_{0}:\:\mu_{1}-\mu_{2}=0$$ $$H_{1}:\:\mu_{1}-\mu_{2}\neq0$$ La prueba
de dos colas (“two-sided”) debe utilizarse.

Primero calculamos la media $\mu_{1}$, $\mu_{2}$.

``` r
# Calculamos la media

# CESC
# GLS1
mean_gls1_exp_tumor_cesc <- mean(primary_tumor_cesc$gls1_gene_expression)
mean_gls1_exp_control_cesc <- mean(solid_tissue_normal_cesc$gls1_gene_expression)
# GLS2
mean_gls2_exp_tumor_cesc <- mean(primary_tumor_cesc$gls2_gene_expression)
mean_gls2_exp_control_cesc <- mean(solid_tissue_normal_cesc$gls2_gene_expression)

# HNSC
# GLS1
mean_gls1_exp_tumor_hnsc <- mean(primary_tumor_hnsc$gls1_gene_expression)
mean_gls1_exp_control_hnsc <- mean(solid_tissue_normal_hnsc$gls1_gene_expression)
# GLS2
mean_gls2_exp_tumor_hnsc <- mean(primary_tumor_hnsc$gls2_gene_expression)
mean_gls2_exp_control_hnsc <- mean(solid_tissue_normal_hnsc$gls2_gene_expression)

# LIHC
# GLS1
mean_gls1_exp_tumor_lihc <- mean(primary_tumor_lihc$gls1_gene_expression)
mean_gls1_exp_control_lihc <- mean(solid_tissue_normal_lihc$gls1_gene_expression)
# GLS2
mean_gls2_exp_tumor_lihc <- mean(primary_tumor_lihc$gls2_gene_expression)
mean_gls2_exp_control_lihc <- mean(solid_tissue_normal_lihc$gls2_gene_expression)
```

Ahora vamos a calcular el **efecto del tamaño** calculado por la fórmula
del **d de Cohen**. Este está dado por la siguiente fórmula (4) (3) :

$$d=\frac{\mu_{1}-\mu_{2}}{\sigma}$$

Si los tamaños muestrales son similares ($n_{1}=n_{2}$), entonces:

$$\sigma=\sqrt{\frac{\sigma_{1}^{2}+\sigma_{2}^{2}}{2}}$$ Sin embargo,
si ambos grupos tienen grandes diferencias de tamaño, $n_{1}\neq n_{2}$,
entonces se usa la desviación estándar entre cada grupo (*pooled
standard deviation*) en lugar de $\sigma$.

Como es el caso, calculamos la desviación estándar entre cada grupo
(*pooled standard deviation*), $\sigma$, con el paquete `effectsize`.

``` r
# Calculate sd
# CESC
# GLS1
sd_pooled_gls1_exp_cesc <- sd_pooled(primary_tumor_cesc$gls1_gene_expression, solid_tissue_normal_cesc$gls1_gene_expression)
# GLS2
sd_pooled_gls2_exp_cesc <- sd_pooled(primary_tumor_cesc$gls2_gene_expression, solid_tissue_normal_cesc$gls2_gene_expression)

# HNSC
# GLS1
sd_pooled_gls1_exp_hnsc <- sd_pooled(primary_tumor_hnsc$gls1_gene_expression, solid_tissue_normal_hnsc$gls1_gene_expression)
# GLS2
sd_pooled_gls2_exp_hnsc <- sd_pooled(primary_tumor_hnsc$gls2_gene_expression, solid_tissue_normal_hnsc$gls2_gene_expression)

# LIHC
# GLS1
sd_pooled_gls1_exp_lihc <- sd_pooled(primary_tumor_lihc$gls1_gene_expression, solid_tissue_normal_lihc$gls1_gene_expression)
# GLS2
sd_pooled_gls2_exp_lihc <- sd_pooled(primary_tumor_lihc$gls2_gene_expression, solid_tissue_normal_lihc$gls2_gene_expression)
```

Ahora sí calculamos el **efecto del tamaño** calculado por la fórmula
del **d de Cohen**, la cual se describió más arriba:

$$d=\frac{\mu_{1}-\mu_{2}}{\sigma}$$

``` r
# Calculate effect size d, tumor = group 1, control = group 2
# CESC
# GLS1
d_pooled_gls1_exp_cesc <- abs((mean_gls1_exp_tumor_cesc-mean_gls1_exp_control_cesc)/sd_pooled_gls1_exp_cesc)
# GLS2
d_pooled_gls2_exp_cesc <- abs((mean_gls2_exp_tumor_cesc-mean_gls2_exp_control_cesc)/sd_pooled_gls2_exp_cesc)

# HNSC
# GLS1
d_pooled_gls1_exp_hnsc <- abs((mean_gls1_exp_tumor_hnsc-mean_gls1_exp_control_hnsc)/sd_pooled_gls1_exp_hnsc)
# GLS2
d_pooled_gls2_exp_hnsc <- abs((mean_gls2_exp_tumor_hnsc-mean_gls2_exp_control_hnsc)/sd_pooled_gls2_exp_hnsc)

# LIHC
# GLS1
d_pooled_gls1_exp_lihc <- abs((mean_gls1_exp_tumor_lihc-mean_gls1_exp_control_lihc)/sd_pooled_gls1_exp_lihc)
# GLS2
d_pooled_gls2_exp_lihc <- abs((mean_gls2_exp_tumor_lihc-mean_gls2_exp_control_lihc)/sd_pooled_gls2_exp_lihc)
```

Finalmente, calculamos el tamaño muestral con el paquete `pwr`, en donde
tenemos un **t-test de dos muestras a priori con dos colas**, $d=$ los
valores calculados arriba, la probabilidad de error $\alpha=0.05$, el
“poder” $1-\beta=0.95$, y el ratio de alocación $N_{2}/N_{1}=1$:

``` r
# Calculate sample size

# CESC
# GLS1
pwr.t.test(d=d_pooled_gls1_exp_cesc,power=0.95,sig.level=0.05,type="two.sample",alternative="two.sided")
# GLS2
pwr.t.test(d=d_pooled_gls2_exp_cesc,power = 0.95,sig.level = 0.05,type="two.sample",alternative="two.sided")

# HNSC
# GLS1
pwr.t.test(d=d_pooled_gls1_exp_hnsc,power=0.95,sig.level=0.05,type="two.sample",alternative="two.sided")
# GLS2
pwr.t.test(d=d_pooled_gls2_exp_hnsc,power = 0.95,sig.level = 0.05,type="two.sample",alternative="two.sided")

# LIHC
# GLS1
pwr.t.test(d=d_pooled_gls1_exp_lihc,power=0.95,sig.level=0.05,type="two.sample",alternative="two.sided")
# GLS2
pwr.t.test(d=d_pooled_gls2_exp_lihc,power = 0.95,sig.level = 0.05,type="two.sample",alternative="two.sided")
```

Ahora sí, el tamaño muestra:

    Total sample size GLS1 in CESC =  238

    Total sample size GLS2 in CESC =  74

    Total sample size GLS1 in HNSC =  32

    Total sample size GLS2 in HNSC =  40902

    Total sample size GLS1 in LIHC =  52

    Total sample size GLS2 in LIHC =  31

## Referencias

(1–4)

<div id="refs" class="references csl-bib-body">

<div id="ref-murray2009" class="csl-entry">

<span class="csl-left-margin">1.
</span><span class="csl-right-inline">Murray S, Larry S. Estadística.
4th ed. México, D.F.: McGraw-Hill; 2009. (Schaum). </span>

</div>

<div id="ref-zhang2023" class="csl-entry">

<span class="csl-left-margin">2.
</span><span class="csl-right-inline">Zhang X, Hartmann P. How to
calculate sample size in animal and human studies. Frontiers in Medicine
\[Internet\]. 2023 Aug 17;10:1215927. Available from:
<https://www.frontiersin.org/articles/10.3389/fmed.2023.1215927/full></span>

</div>

<div id="ref-gpower32023" class="csl-entry">

<span class="csl-left-margin">3.
</span><span class="csl-right-inline">GPower 3.1 manua. 2023 Jun 1;
Available from:
<https://www.psychologie.hhu.de/fileadmin/redaktion/Fakultaeten/Mathematisch-Naturwissenschaftliche_Fakultaet/Psychologie/AAP/gpower/GPowerManual.pdf></span>

</div>

<div id="ref-cohen1988" class="csl-entry">

<span class="csl-left-margin">4.
</span><span class="csl-right-inline">Cohen J. The t Test for Means. In:
2nd ed. New York: Lawrence Erlbaum Associates; 1988. p. 20. </span>

</div>

</div>
