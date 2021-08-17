## Download transcriptomic data from GEO platform


### Introduction

It's very confused for me to get the data from GEO and process it before doing the data analysis, since I am lack of patient for the numerous repetitions on the Excel.

Luckily, the Package from R sovled this terrible situation. The following details would display how to figure out this problem.


### How it works

Firstly, loading GEOquery pacakges for using the getGEO function, and then the object including clinical information and expression profile would be generate. Morever, performing preprocess step on clinical and expression data would extract all the data we need for the subsequant analysis. Finally, transforming phenotype and profile into ExpressionSet object for further data analysis.


### Command line 

The microarray seems have different regulations including *LIMN* and *array* on naming the probes. 

```bash
Rscript Download_GEOData_v1.R -g GSE65858 -p GPL10558 -t ILMN -o ./
Rscript Download_GEOData_v1.R -g GSE55457 -p GPL96 -t array -o ./
```

### Files' Structure

Downloading the Gene expression profile from GEO

```bash
GSE55457
├── GPL96.soft
├── GSE55457_series_matrix.txt.gz
├── process
│   ├── GPL96_probe2gene_table.tsv
│   ├── GSE55457_GeneExprSet.RDS
│   ├── GSE55457_clinical_origin.csv
│   ├── GSE55457_clinical_post.csv
│   ├── GSE55457_profile_origin.tsv
│   └── GSE55457_profile_post.tsv
└── work.sh
```


### Update 8/16/2021

Since there is R package, providing the relationship between probes and gene symbol, I updated the pipeline. The new version download script named Download_GEOData_v2.R is finished and the usage as following:

```bash
Rscript Download_GEOData_v2.R GSE30552 -p GPL6246 -t other -o ./
```
