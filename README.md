<h1 align="center"> README</h1>


<div align="center">
<img src="https://img.shields.io/github/last-commit/livisupajares/undergrad-thesis?style=for-the-badge&logo=github&color=F5C2E7&logoColor=CDD6F4&labelColor=313244"/>
<img src="https://img.shields.io/github/repo-size/livisupajares/undergrad-thesis?style=for-the-badge&logo=github&color=CBA6F7&logoColor=CDD6F4&labelColor=313244"/>
<img src="https://img.shields.io/github/languages/top/livisupajares/undergrad-thesis?style=for-the-badge&logo=lua&color=94E2D5&logoColor=CDD6F4&labelColor=313244"/>
</div>
<br>
<div align="center">
<table>
    <tr>
        <th>
            <a href="https://github.com/livisupajares/undergrad-thesis/tree/master">Español 🇪🇸</a>
        <th>
        <!-- TODO: Change href url to .md after rendering in quarto -->
            <a href="https://github.com/livisupajares/undergrad-thesis/blob/master/docs/en/sample-size/calculation-sample-size_en.qmd">English 🇬🇧</a>
    </tr>
</table>
</div>

## Tabla de Contenidos

  - [1. Recolección de bases de datos](#1-recolecci%C3%B3n-de-bases-de-datos)
    - [Selección de variables](#selecci%C3%B3n-de-variables)
      - [Variables genómicas](#variables-gen%C3%B3micas)
      - [Variables fenotípicas](#variables-fenot%C3%ADpicas)
      - [Variables Secundarias](#variables-secundarias)
    - [Cálculo del tamaño muestral](#c%C3%A1lculo-del-tama%C3%B1o-muestral)
  - [2. Procesamiento de la data original](#2-procesamiento-de-la-data-original)

## 1. Recolección de bases de datos

Se descargaron bases de datos del [XenaBrowser](https://xenabrowser.net/) y se seleccionaron los siguientes estudios. El número de sujetos disponibles por cada cohorte es el siguiente:

- GDC TCGA Cervical Cancer (CESC): 317
- GDC TCGA Head and Neck Cancer (HNSC): 612
- GDC TCGA Liver Cancer (LIHC): 469

### Selección de variables

<!-- TODO: Make 3 simple markdown tables -->

#### Variables genómicas

|Para las tres cohortes (CESC, HNSC, LIHC)|
|-----------------------------------------|
|*Copy Number Segment* GLS1/GLS2|
|*HTSeq-FPKM (RNASeq)* GLS1/GLS2|

#### Variables fenotípicas

|CESC|HNSC|LIHC|
|---|----|----|
|sample_type.samples|sample_type.samples|sample_type.samples|
|neoplasm_histologic_grade|neoplasm_histologic_grade|neoplasm_histologic_grade|tumor_stage.diagnoses|
|clinical_stage|clinical_stage|
|OS.time|OS.time|OS.time|
|tobacco_smoking_history|tobacco_smoking_history| - |
| - |hpv_status_by_ish_testing| - |
| - |alcohol_history.exposures| - |

#### Variables Secundarias

> 📝 Estas variables sirven para describir la población. No se usarán en los análisis.

|CESC|HNSC|LIHC|
|---|----|----|
|age_at_initial_pathologic_diagnosis|age_at_initial_pathologic_diagnosis|age_at_initial_pathologic_diagnosis|
|gender.demographic|gender.demographic|gender.demographic|
|race.demographic|race.demographic|race.demographic|
|ethnicity.demographic|ethnicity.demographic|ethnicity.demographic|
|bmi.exposures| - |bmi.exposures|

### Cálculo del tamaño muestral

El proceso de cálculo para el tamaño muestral está documentado [aquí](https://github.com/livisupajares/undergrad-thesis/blob/master/docs/es/sample-size/calculation-sample-size_es.md).

