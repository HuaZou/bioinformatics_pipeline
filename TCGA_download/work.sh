Rscript 01.Download.R -t TCGA-UCS -o mRNA
Rscript 02.Prerpocess_clinical.R -p Clinical/TCGA-UCS_clinical_origin.csv -t TCGA-UCS
Rscript 03.Prerpocess_omics.R -e Omics/TCGA-UCS_mRNA.RDS -t mRNA -p TCGA-UCS-post
Rscript 04.ExpressionSet.R -p Clinical/TCGA-UCS-post_clinical.csv -e Clean/TCGA-UCS-post_mRNA_profile.tsv -t mRNA -o TCGA-UCS-post
