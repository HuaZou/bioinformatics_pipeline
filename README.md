## Bioinformatics pipeline 



### Introduction

An effective, reproducible and  reliable data analysis workflow is based on the state-of-the-art pipeline, using the most up-to-date methods. In order to  facilitate my future work, building the bioinformatics pipeline is necessary. However, how to construct the workflow is still hard for me. In my view, I wanna use the nextflow or snakemake program language to do it. Before doing it, what I need is to describe the workflows in the mindmap which could make me be clear.



### Interests

* amplicon sequencing analysis
* metagenomics sequencing analysis
* bulk-RNA sequencing analysis



### Workflow

The workflow comprises of two parts, one is from raw data to profile, and the other is data analysis (statistical analysis)  

*  the first parts
  * demultiplex sequences; 
  * scan the quality of reads; 
  * filter the low quality reads and remove host DNA sequence
  * align the high quality reads into the reference database
  * obtain the profile whose structure is $M x N$ matrix (M: features' name; N: sampleid)
* the second part
  * statistical analysis such as wilcoxon rank sum test, LDA, PCoA, linear regression analysis and multivariables association analysis
  * machine learning



### Notice

First and foremost, I utilize the perl program language to do some preliminary work and finally convert all the workflows into snakemake or nextflow.





  