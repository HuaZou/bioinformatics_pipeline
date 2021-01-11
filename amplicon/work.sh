conda activate qiime2-2020.8
# debug
## summarize_fastq.smk
snakemake -np -r --debug-dag result/02-remove/paired-end-demux-trim-fastq_counts_describe.md  --configfile config.yaml --snakefile Snakefile
## denoise.smk
snakemake -np -r --debug-dag result/03-denoise/repseqs_lengths_describe.md --configfile config.yaml --snakefile Snakefile
## taxonomy 
snakemake -np -r --debug-dag result/04-taxonomy/taxonomy_final.tsv --configfile config.yaml --snakefile Snakefile
## tree 
snakemake -np -r --debug-dag result/05-alignment-tree/rooted_tree.qza --configfile config.yaml --snakefile Snakefile
## alpha
snakemake -np -r --debug-dag result/06-alpha-diversity/evenness_group_significance.qzv --configfile config.yaml --snakefile Snakefile
## report
snakemake -np -r --debug-dag result/08-report/report.html --configfile config.yaml --snakefile Snakefile

# wrokflow 
snakemake --dag result/08-report/report.html | dot -Tsvg > workflow.svg


# run
snakemake result/08-report/report.html --configfile config.yaml --snakefile Snakefile --cores 10
