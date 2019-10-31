snakemake -np -r --debug-dag result/00.quality/multiqc/fastqc_multiqc_report.html
snakemake -np -r --debug-dag result/01.trimmed/summary_trimmed.fa
snakemake -np -r --debug-dag result/02.deredundancy/sorted.fa
snakemake -np -r --debug-dag result/03.clusterOTU/OTU.cluster.fa
snakemake -np -r --debug-dag result/03.clusterOTU/OTU.unoise3.fa
snakemake -np -r --debug-dag result/04.rechimeras/OTU.rechimera_silva.cluster.fa
snakemake -np -r --debug-dag result/05.OTUtable/OTU_table.silva.cluster.txt
snakemake -np -r --debug-dag result/06.taxonomy/taxonomy.silva.cluster.txt

snakemake --dag result/06.taxonomy/taxonomy.silva.cluster.txt | dot -Tsvg > workflow.svg

snakemake result/06.taxonomy/taxonomy.silva.cluster.txt
snakemake result/00.quality/multiqc/fastqc_multiqc_report.html
