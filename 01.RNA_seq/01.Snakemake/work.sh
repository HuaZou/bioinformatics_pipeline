snakemake -np -r --debug-dag result/00.quality/multiqc/fastqc_multiqc_report.html
snakemake -np -r --debug-dag result/01.trimmed/HBR_Rep1_1.fastq.gz
snakemake -np -r --debug-dag result/02.index/chr22_with_ERCC92.1.ht2
snakemake -np -r --debug-dag result/03.alignment/HBR_Rep1.sam
snakemake -np -r --debug-dag result/03.alignment/HBR_Rep1.sorted.bam
snakemake -np -r --debug-dag result/03.alignment/HBR.merged.bam
