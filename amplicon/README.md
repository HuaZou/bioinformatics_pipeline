## Amplicon Sequence Analysis Workflow


Build this pipeline based on the Tourmaline via snakemake langugae. I reorganized the structures and added some other plugins for data analysis. Qiime2 platform and snakemake makes the effective and extensible standard workflow be possible.

![](workflow.svg)

### Installation 

Before run this pipeline, you should install Qiime2, Snakemake, and some other dependent packages or softwares. All the softwares could be installed by conda commands.

#### install qiime2

I tried use `conda search qiime2` command to find the suitable version, but in fact it failed. Using the configure file from the Qiime2 offical website is the best option. (the following commands just copied from Tourmaline)

```bash
wget https://data.qiime2.org/distro/core/qiime2-2020.8-py36-osx-conda.yml
conda env create -n qiime2-2020.8 --file qiime2-2020.8-py36-osx-conda.yml -y
conda activate qiime2-2020.8
conda install -c bioconda snakemake biopython tabulate pandoc tabview -y 
pip install git+https://github.com/biocore/empress.git
qiime dev refresh-cache
```

#### install R packages

Running R in the linux terminal and then using the *install()* function from the `BiocManger` package to install the dependent packages

```R
packages <- c("msa", "odseq")
BiocManager::install(packages)
```

### Download Reference database 

Compared to greengene, silva could be much better for taxonomic annotation. 

```bash
wget https://data.qiime2.org/2020.8/common/silva-138-99-seqs-515-806.qza
wget https://data.qiime2.org/2020.8/common/silva-138-99-tax-515-806.qza
ln -s silva-138-99-seqs-515-806.qza refseqs.qza
ln -s silva-138-99-tax-515-806.qza reftax.qza
```


###  pipeline directory 

the pipeline comprised of three main parts: **Snakemake file, Rules and Scripts**. the following directory structure shows the final results of this pipeline.

```bash
amplicon/
├── config.yaml
├── css
│   ├── github.css
│   ├── gothic.css
│   ├── newsprint.css
│   ├── night.css
│   ├── pixyll.css
│   └── whitey.css
├── manifest.csv
├── metadata.tsv
├── README.md
├── repseqs_to_filter.tsv
├── result
│   ├── 01-import
│   │   ├── paired-end-demux-fastq_counts_describe.md
│   │   ├── paired-end-demux-fastq_counts.tsv
│   │   ├── paired-end-demux-fastq_summary.qzv
│   │   ├── paired-end-demux.qza
│   │   ├── refseqs.qza
│   │   └── reftax.qza
│   ├── 02-remove
│   │   ├── paired-end-demux-trim-fastq_counts_describe.md
│   │   ├── paired-end-demux-trim-fastq_counts.tsv
│   │   ├── paired-end-demux-trim-fastq_summary.qzv
│   │   └── paired-end-demux-trim.qza
│   ├── 03-denoise
│   │   ├── chimeras.qza
│   │   ├── chimeras_stats.qza
│   │   ├── dada2_stats.qza
│   │   ├── nonchimera.qza
│   │   ├── repseq_final.qza
│   │   ├── repseqs_amplicon_type.txt
│   │   ├── repseqs_final.fasta
│   │   ├── repseqs_final.qza
│   │   ├── repseqs_final.qzv
│   │   ├── repseqs_lengths_describe.md
│   │   ├── repseqs_lengths.txt
│   │   ├── repseqs_nonchimera_contam.qza
│   │   ├── repseqs_nonchimera.qza
│   │   ├── repseqs.qza
│   │   ├── table.biom
│   │   ├── table_final.qza
│   │   ├── table_final.qzv
│   │   ├── table_nonchimera.qza
│   │   ├── table.qza
│   │   ├── table_summary_features.txt
│   │   └── table_summary_samples.txt
│   ├── 04-taxonomy
│   │   ├── classifier.qza
│   │   ├── taxa_barplot.qzv
│   │   ├── taxonomy.qza
│   │   ├── taxonomy.qzv
│   │   └── taxonomy.tsv
│   ├── 05-alignment-tree
│   │   ├── aligned_repseqs.fasta
│   │   ├── aligned_repseqs_gaps.txt
│   │   ├── aligned_repseqs_outliers.tsv
│   │   ├── aligned_repseqs.qza
│   │   ├── outliers.qza
│   │   ├── outliers.tsv
│   │   ├── repseqs_properties.pdf
│   │   ├── repseqs_properties.tsv
│   │   ├── repseqs_to_filter_outliers.tsv
│   │   ├── repseqs_to_filter_unassigned.tsv
│   │   ├── rooted_tree.qza
│   │   ├── rooted_tree.qzv
│   │   ├── unmasked_aligned_repseqs.qza
│   │   └── unrooted_tree.qza
│   ├── 06-alpha-diversity
│   │   ├── alpha_rarefaction.qzv
│   │   ├── evenness_group_significance.qzv
│   │   ├── evenness_vector.qza
│   │   ├── faith_pd_group_significance.qzv
│   │   ├── faith_pd_vector.qza
│   │   ├── observed_features_group_significance.qzv
│   │   ├── observed_features_vector.qza
│   │   ├── rarefied_table.qza
│   │   ├── shannon_group_significance.qzv
│   │   └── shannon_vector.qza
│   ├── 07-beta-diversity
│   │   ├── bray_curtis_distance_matrix.qza
│   │   ├── bray_curtis_emperor.qzv
│   │   ├── bray_curtis_group_significance.qzv
│   │   ├── bray_curtis_pcoa_results.qza
│   │   ├── jaccard_distance_matrix.qza
│   │   ├── jaccard_emperor.qzv
│   │   ├── jaccard_group_significance.qzv
│   │   ├── jaccard_pcoa_results.qza
│   │   ├── unweighted_unifrac_distance_matrix.qza
│   │   ├── unweighted_unifrac_emperor.qzv
│   │   ├── unweighted_unifrac_group_significance.qzv
│   │   ├── unweighted_unifrac_pcoa_results.qza
│   │   ├── weighted_unifrac_distance_matrix.qza
│   │   ├── weighted_unifrac_emperor.qzv
│   │   ├── weighted_unifrac_group_significance.qzv
│   │   └── weighted_unifrac_pcoa_results.qza
│   ├── 08-report
│   │   ├── metadata_summary.md
│   │   ├── report.html
│   │   └── report.md
│   └── logs
│       ├── 02-remove
│       ├── 03-denoise
│       ├── 04-taxonomy
│       ├── 05-alignment-tree
│       ├── 06-alpha-diversity
│       ├── 07-beta-diversity
│       └── 08-report
├── rules
│   ├── denoise.smk
│   ├── diversity.smk
│   ├── import.smk
│   ├── remove_chimera_contamination.smk
│   ├── remove_primer.smk
│   ├── report.smk
│   ├── summarize_fastq.smk
│   ├── summarize_table.smk
│   ├── taxonomy.smk
│   └── tree.smk
├── scripts
│   ├── create_manifest_from_fastq_directory.py
│   ├── detect_amplicon_locus.py
│   ├── fastaLengths.pl
│   ├── fastqc_per_base_sequence_quality_dropoff.py
│   ├── match_manifest_to_metadata.py
│   └── run_odseq.R
├── Snakefile
├── workflow.svg
└── work.sh

20 directories, 109 files
```


### How to run 

```bash
conda activate qiime2-2020.8
# debug
snakemake -np -r --debug-dag result/03-denoise/table_final.qza --configfile config.yaml --snakefile Snakefile

# wrokflow figure
snakemake --dag --debug-dag result/03-denoise/table_final.qza | dot -Tsvg > workflow.svg

# run
snakemake result/03-denoise/table_final.qza --configfile config.yaml --snakefile Snakefile --cores 2
```



### Development 

In the next time, I wanna to add some other common or specific procedures in amplicon sequencing or adjust it to be applied for single end reads. 

Other analysis procedures:

* ANCOM: analysis of composition
* LEfse: linear discriminant Effect Size
* Picrust2: Phylogenetic investigating of Communities by Reconstruction of  Unobserved States
* SourceTracter: Tracking the source of microbes with SourceTracker
* iCAMP: Infer Community Assembly Mechanisms by Phylogenetic bin-based null model analysis 
* etc



### Contributors

* **Hua Zou**



### Reference 

Thanks for [The advanced Amplicon Sequence Processing Workflow](https://github.com/lukenoaa/tourmaline).