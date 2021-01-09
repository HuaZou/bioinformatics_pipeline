rule feature_classifier:
    input:
        refseqs = os.path.join(config["result"]["import"], "refseqs.qza"),
        reftax  = os.path.join(config["result"]["import"], "reftax.qza"),
        classifier = os.path.join(config["result"]["taxonomy"], "classifier.qza"),
        repseqs    = os.path.join(config["result"]["denoise"], "repseqs.qza")
    params:
        classifymethod = config["params"]["classifier"]["classify_method"]
    threads: 
        config["params"]["classifier"]["feature_classifier_threads"]
    log:
        os.path.join(config["logs"]["taxonomy"], "taxonomy.log")
    output:
        os.path.join(config["result"]["taxonomy"], "taxonomy.qza")
    shell:
        '''
        echo classify_method: {params.classifymethod}; 
        if [ {params.classifymethod} = 'naive-bayes' ]; then \
        qiime feature-classifier classify-sklearn \
            --i-classifier {input.classifier} \
            --i-reads {input.repseqs} \
            --o-classification {output} \
            --p-n-jobs {threads}; 
        elif [ {params.classifymethod} = 'consensus-blast' ]; then \
        qiime feature-classifier classify-consensus-blast \
            --i-reference-reads {input.refseqs} \
            --i-reference-taxonomy {input.reftax} \
            --i-query {input.repseqs} \
            --o-classification {output}; 
        fi 2>{log}
        '''       

rule visualize_taxonomy:
    input:
        os.path.join(config["result"]["taxonomy"], "taxonomy.qza")
    output:
        os.path.join(config["result"]["taxonomy"], "taxonomy.qzv")
    shell:
        '''
        qiime metadata tabulate \
            --m-input-file {input} \
            --o-visualization {output}
        '''

rule taxa_barplot:
    input:
        table    = os.path.join(config["result"]["denoise"], "table.qza"),
        taxonomy = os.path.join(config["result"]["taxonomy"], "taxonomy.qza"),
        metadata = config["metadata"]
    output:
        os.path.join(config["result"]["taxonomy"], "taxa_barplot.qzv")
    shell:
        '''
        qiime taxa barplot \
            --i-table {input.table} \
            --i-taxonomy {input.taxonomy} \
            --m-metadata-file {input.metadata} \
            --o-visualization {output}
        '''

rule unzip_taxonomy_to_tsv:
    input:
        os.path.join(config["result"]["taxonomy"], "taxonomy.qza")
    output:
        os.path.join(config["result"]["taxonomy"], "taxonomy.tsv")
    shell:
        '''
        unzip -qq -o {input} -d temp; 
        mv temp/*/data/taxonomy.tsv {output}; 
        rm -r temp
        '''

rule import_taxonomy_to_qza:
    input:
        os.path.join(config["result"]["taxonomy"], "taxonomy.tsv")
    output:
        os.path.join(config["result"]["taxonomy"], "taxonomy.qza")
    shell:
        '''
        qiime tools import \
            --type 'FeatureData[Taxonomy]' \
            --input-format TSVTaxonomyFormat \
            --input-path {input} \
            --output-path {output} 
        '''
