## Download Clinical information and Omic-data profile from TCGA project

### Introduction

TCGA project had produced multiple omics data in Cancer research, and the public data could be assessed by R program. This pipeline aims to acquire the clinical information and omic-data profile and dissect the information into analysis model.


### how it works

There are four scripts which separated the pipeline into four parts

1. **01.Download.R**: download clinical information and save omics data into SummariedExperimentData

2. **02.Prerpocess_clinical.R**: Process the clinical information

3. **03.Prerpocess_omics.R**: Process the omics data profile

4. **04.ExpressionSet.R**: Convert Phenotype from clincal information and Profile into ExpressionSet

### Summary

The clinical information contains all the patients' phenotype but it doesn't show any Normal samples' information derived from omics data profile. Therefore, the extra information about normal samples would be obtains from the SummariedExperimentData. In the other hand, to distinct the tumor and normal samples is also be cautious, the Number named "01" at the position from 14 to 15 characters of SummarizedExperiment object's rownames is "tumor" and the others are "normal".
