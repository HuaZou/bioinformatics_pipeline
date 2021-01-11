rule diversity_alpha_rarefaction:
    input:
        table     = os.path.join(config["result"]["denoise"], "table_final.qza"),
        phylogeny = os.path.join(config["result"]["tree"], "rooted_tree.qza"),
        metadata  = config["metadata"]
    params:
        maxdepth  = config["params"]["diversity"]["alpha_max_depth"]
    log:
        os.path.join(config["logs"]["alpha"], "alpha_rarefaction.log")    
    output:
        os.path.join(config["result"]["alpha"], "alpha_rarefaction.qzv")
    shell:
        '''
        qiime diversity alpha-rarefaction \
            --i-table {input.table} \
            --i-phylogeny {input.phylogeny} \
            --p-max-depth {params.maxdepth} \
            --p-metrics faith_pd \
            --p-metrics shannon \
            --p-metrics observed_features \
            --p-metrics chao1 \
            --m-metadata-file {input.metadata} \
            --o-visualization {output}  2>{log}
        '''

rule diversity_core_metrics_phylogenetic:
    input:
        table     = os.path.join(config["result"]["denoise"], "table_final.qza"),
        phylogeny = os.path.join(config["result"]["tree"], "rooted_tree.qza"),
        metadata  = config["metadata"]
    params:
        samplingdepth = config["params"]["diversity"]["core_sampling_depth"]
    threads: 
        config["params"]["diversity"]["diversity_core_metrics_phylogenetic_threads"]
    log:
        os.path.join(config["logs"]["alpha"], "diversity_core_metrics_phylogenetic.log") 
    output: 
        rarefiedtable = os.path.join(config["result"]["alpha"], "rarefied_table.qza"),
        faithpdvector = os.path.join(config["result"]["alpha"], "faith_pd_vector.qza"),
        observedfeaturesvector = os.path.join(config["result"]["alpha"], "observed_features_vector.qza"),
        shannonvector  = os.path.join(config["result"]["alpha"], "shannon_vector.qza"),
        evennessvector = os.path.join(config["result"]["alpha"], "evenness_vector.qza"),
        unweightedunifracdistancematrix = os.path.join(config["result"]["beta"], "unweighted_unifrac_distance_matrix.qza"),
        weightedunifracdistancematrix   = os.path.join(config["result"]["beta"], "weighted_unifrac_distance_matrix.qza"),
        jaccarddistancematrix    = os.path.join(config["result"]["beta"], "jaccard_distance_matrix.qza"),
        braycurtisdistancematrix = os.path.join(config["result"]["beta"], "bray_curtis_distance_matrix.qza"),
        unweightedunifracpcoaresults = os.path.join(config["result"]["beta"], "unweighted_unifrac_pcoa_results.qza"),
        weightedunifracpcoaresults = os.path.join(config["result"]["beta"], "weighted_unifrac_pcoa_results.qza"),
        jaccardpcoaresults    = os.path.join(config["result"]["beta"], "jaccard_pcoa_results.qza"),
        braycurtispcoaresults = os.path.join(config["result"]["beta"], "bray_curtis_pcoa_results.qza"),
        unweightedunifracemperor = os.path.join(config["result"]["beta"], "unweighted_unifrac_emperor.qzv"),
        weightedunifracemperor   = os.path.join(config["result"]["beta"], "weighted_unifrac_emperor.qzv"),
        jaccardemperor    = os.path.join(config["result"]["beta"], "jaccard_emperor.qzv"),
        braycurtisemperor = os.path.join(config["result"]["beta"], "bray_curtis_emperor.qzv")
    shell:
        '''
        qiime diversity core-metrics-phylogenetic \
            --i-table {input.table} \
            --i-phylogeny {input.phylogeny} \
            --p-sampling-depth {params.samplingdepth} \
            --m-metadata-file {input.metadata} \
            --o-rarefied-table {output.rarefiedtable} \
            --o-faith-pd-vector {output.faithpdvector} \
            --o-observed-features-vector {output.observedfeaturesvector} \
            --o-shannon-vector {output.shannonvector} \
            --o-evenness-vector {output.evennessvector} \
            --o-unweighted-unifrac-distance-matrix {output.unweightedunifracdistancematrix} \
            --o-weighted-unifrac-distance-matrix {output.weightedunifracdistancematrix} \
            --o-jaccard-distance-matrix {output.jaccarddistancematrix} \
            --o-bray-curtis-distance-matrix {output.braycurtisdistancematrix} \
            --o-unweighted-unifrac-pcoa-results {output.unweightedunifracpcoaresults} \
            --o-weighted-unifrac-pcoa-results {output.weightedunifracpcoaresults} \
            --o-jaccard-pcoa-results {output.jaccardpcoaresults} \
            --o-bray-curtis-pcoa-results {output.braycurtispcoaresults} \
            --o-unweighted-unifrac-emperor {output.unweightedunifracemperor} \
            --o-weighted-unifrac-emperor {output.weightedunifracemperor} \
            --o-jaccard-emperor {output.jaccardemperor} \
            --o-bray-curtis-emperor {output.braycurtisemperor} \
            --p-n-jobs-or-threads {threads} 2>{log}
        '''

rule diversity_alpha_group_significance:
    input:
        alphadiversity = os.path.join(config["result"]["alpha"], "{metric}_vector.qza"),
        metadata       = config["metadata"]
    output:
        os.path.join(config["result"]["alpha"], "{metric}_group_significance.qzv")
    shell:
        '''
        qiime diversity alpha-group-significance \
            --i-alpha-diversity {input.alphadiversity} \
            --m-metadata-file {input.metadata} \
            --o-visualization {output}
        '''

rule diversity_beta_group_significance:
    input:
        distancematrix = os.path.join(config["result"]["beta"], "{metric}_distance_matrix.qza"),
        metadata       = config["metadata"]
    params:
        column = config["params"]["diversity"]["beta_group_column"]
    log:
        os.path.join(config["logs"]["beta"], "{metric}_group_significance.log") 
    output:
        os.path.join(config["result"]["beta"], "{metric}_group_significance.qzv")
    shell:
        '''
        qiime diversity beta-group-significance \
            --i-distance-matrix {input.distancematrix} \
            --m-metadata-file {input.metadata} \
            --m-metadata-column {params.column} \
            --o-visualization {output} \
            --p-pairwise 2>{log}
        '''
