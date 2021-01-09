rule import_ref_seqs:
    input:
        config["params"]["classifier"]["refseqs"]
    output:
        os.path.join(config["result"]["import"], "refseqs.qza")
    shell:
        '''
        qiime tools import \
            --type 'FeatureData[Sequence]' \
            --input-path {input} \
            --output-path {output}
        '''

rule import_ref_tax:
    input:
        config["params"]["classifier"]["reftax"]
    output:
        os.path.join(config["result"]["import"], "reftax.qza")
    shell:
        '''
        qiime tools import \
            --type 'FeatureData[Taxonomy]' \
            --input-format HeaderlessTSVTaxonomyFormat \
            --input-path {input} \
            --output-path {output}
        '''

rule build_classifier:
    input:
        refseqs = os.path.join(config["result"]["import"], "refseqs.qza"),
        reftax  = os.path.join(config["result"]["import"], "reftax.qza")
    output:
        os.path.join(config["result"]["taxonomy"], "classifier.qza")
    shell:
        '''
        qiime feature-classifier fit-classifier-naive-bayes \
            --i-reference-reads {input.refseqs} \
            --i-reference-taxonomy {input.reftax} \
            --o-classifier {output}
        '''

rule import_fastq_demux_pe:
    input:
        config["manifest"]
    log:
        os.path.join(config["logs"]["import"], "fastq2qza.log")
    output:
        os.path.join(config["result"]["import"], "paired-end-demux.qza")
    shell:
        '''
        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path {input} \
            --output-path {output} \
            --input-format PairedEndFastqManifestPhred33 2>{log}
        '''
