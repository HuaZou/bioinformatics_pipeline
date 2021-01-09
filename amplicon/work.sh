conda activate qiime2-2020.8
# debug
## summarize_fastq.smk
snakemake -np -r --debug-dag result/02-remove/paired-end-demux-trim-fastq_counts_describe.md  --configfile config.yaml --snakefile Snakefile
## denoise.smk
snakemake -np -r --debug-dag result/03-denoise/table_final.qza --configfile config.yaml --snakefile Snakefile
## taxonomy 
snakemake -np -r --debug-dag result/04-taxonomy/taxonomy_final.tsv--configfile config.yaml --snakefile Snakefile

# wrokflow 
snakemake --dag result/03-denoise/table_final.qza | dot -Tsvg > workflow.svg





# run
snakemake result/03-denoise/table_final.qza --configfile config.yaml --snakefile Snakefile --cores 2
