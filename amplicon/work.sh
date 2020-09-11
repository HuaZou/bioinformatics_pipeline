snakemake -np -r --debug-dag result/07.multiqc/multiqc_report.html
snakemake result/07.multiqc/multiqc_report.html --cores 20

snakemake result/05.finaltable/OTU.table.summary.tsv --cores 20 
