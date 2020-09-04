# test
snakemake -np -r --debug-dag result/06.multiqc/multiqc_report.html

# workflow graph
snakemake --dag result/06.multiqc/multiqc_report.html | dot -Tsvg > hisat2_workflow.svg

# final result
snakemake result/06.multiqc/multiqc_report.html  --configfile config.yaml --snakefile snakefile
