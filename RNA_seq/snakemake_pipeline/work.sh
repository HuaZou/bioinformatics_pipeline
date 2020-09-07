# hisat2
snakemake -np -r --debug-dag hisat2_result/06.multiqc/multiqc_report.html --configfile config_hisat2.yaml --snakefile snakefile_hisat2_pe
snakemake --dag hisat2_result/06.multiqc/multiqc_report.html --configfile config_hisat2.yaml --snakefile snakefile_hisat2_pe | dot -Tsvg > hisat2_workflow.svg
snakemake hisat2_result/06.multiqc/multiqc_report.html  --configfile config_hisat2.yaml --snakefile snakefile_hisat2_pe --cores 20

# star
snakemake -np -r --debug-dag star_result/05.multiqc/multiqc_report.html --configfile config_star.yaml --snakefile snakefile_star_pe
snakemake --dag star_result/05.multiqc/multiqc_report.html --configfile config_star.yaml --snakefile snakefile_star_pe | dot -Tsvg > hisat2_workflow.svg
snakemake star_result/05.multiqc/multiqc_report.html  --configfile config_star.yaml --snakefile snakefile_star_pe --cores 20

# kallisto
snakemake -np -r --debug-dag kallisto_result/03.quantify/kallisto/HBR_Rep1/abundance.tsv --configfile config_kallisto.yaml --snakefile snakefile_kallisto_pe
snakemake kallisto_result/03.quantify/kallisto/HBR_Rep1/abundance.tsv --configfile config_kallisto.yaml --snakefile snakefile_kallisto_pe --cores 20
