snakemake -np -r --debug-dag result/01.trimmed/Sample_path.qiime2.tsv
snakemake -np -r --debug-dag result/02.fastq2qza/qiime2.trimmed.qza
snakemake -np -r --debug-dag result/03.table/table.dada2.qzv  #(result/03.table/biom;result/04.tree/rooted_tree.qza)

