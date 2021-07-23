Rscript 01.Download.R -t TCGA-PAAD -o CNV
Rscript 01.Download.R -t TCGA-PAAD -o mRNA
Rscript 01.Download.R -t TCGA-PAAD -o miRNA
Rscript 01.Download.R -t TCGA-PAAD -o DNA_Methylation

Rscript 02.Prerpocess_clinical.R -p Clinical/TCGA-PAAD_clinical_origin.csv -t TCGA-PAAD

Rscript 03.Prerpocess_omics.R -e Omics/TCGA-PAAD_CNV.RDS -t CNV -g util/snp6.na35.remap.hg38.subset.txt.gz -p TCGA-PAAD-post
Rscript 03.Prerpocess_omics.R -e Omics/TCGA-PAAD_mRNA.RDS -t mRNA -g util/human_gene_all.tsv -p TCGA-PAAD-post
Rscript 03.Prerpocess_omics.R -e Omics/TCGA-PAAD_miRNA.RDS -t miRNA -p TCGA-PAAD-post
Rscript 03.Prerpocess_omics.R -e Omics/TCGA-PAAD_DNA_Methylation.RDS -t DNA_Methylation -p TCGA-PAAD-post

Rscript 04.ExpressionSet.R -p Clinical/TCGA-PAAD-post_clinical.csv -e Clean/TCGA-PAAD-post_mRNA_profile.tsv -t mRNA -o TCGA-PAAD-post
Rscript 04.ExpressionSet.R -p Clinical/TCGA-PAAD-post_clinical.csv -e Clean/TCGA-PAAD-post_miRNA_PKM.tsv -t miRNA -o TCGA-PAAD-post-PKM
Rscript 04.ExpressionSet.R -p Clinical/TCGA-PAAD-post_clinical.csv -e Clean/TCGA-PAAD-post_miRNA_count.tsv -t miRNA -o TCGA-PAAD-post-count
