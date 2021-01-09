rule denoise_dada2_pe:
    input:
        os.path.join(config["result"]["remove"], "paired-end-demux-trim.qza")
    params:
        trunclenf = config["params"]["denoise"]["dada2pe_trunc_len_f"],
        trunclenr = config["params"]["denoise"]["dada2pe_trunc_len_r"],
        trimleftf = config["params"]["denoise"]["dada2pe_trim_left_f"],
        trimleftr = config["params"]["denoise"]["dada2pe_trim_left_r"],
        maxeef = config["params"]["denoise"]["dada2pe_max_ee_f"],
        maxeer = config["params"]["denoise"]["dada2pe_max_ee_r"],
        truncq = config["params"]["denoise"]["dada2pe_trunc_q"],
        poolingmethod = config["params"]["denoise"]["dada2pe_pooling_method"],        
        chimeramethod = config["params"]["denoise"]["dada2pe_chimera_method"],
        minfoldparentoverabundance = config["params"]["denoise"]["dada2pe_min_fold_parent_over_abundance"],
        nreadslearn = config["params"]["denoise"]["dada2pe_n_reads_learn"],
        hashedfeatureids = config["params"]["denoise"]["dada2pe_hashed_feature_ids"]
    threads: 
        config["params"]["denoise"]["dada2pe_threads"]
    log:
        os.path.join(config["logs"]["denoise"], "denoise_dada2_pe.log")
    output:
        table   = os.path.join(config["result"]["denoise"], "table.qza"),
        repseqs = os.path.join(config["result"]["denoise"], "repseqs.qza"),
        stats   = os.path.join(config["result"]["denoise"], "dada2_stats.qza")
    shell:
        '''
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs {input} \
            --p-trunc-len-f {params.trunclenf} \
            --p-trunc-len-r {params.trunclenr} \
            --p-trim-left-f {params.trimleftf} \
            --p-trim-left-r {params.trimleftr} \
            --p-max-ee-f {params.maxeef} \
            --p-max-ee-r {params.maxeer} \
            --p-trunc-q {params.truncq} \
            --p-pooling-method {params.poolingmethod} \
            --p-chimera-method {params.chimeramethod} \
            --p-min-fold-parent-over-abundance {params.minfoldparentoverabundance} \
            --p-n-reads-learn {params.nreadslearn} {params.hashedfeatureids} \
            --p-n-threads {threads} \
            --o-table {output.table} \
            --o-representative-sequences {output.repseqs} \
            --o-denoising-stats {output.stats} \
            --verbose 2>{log}
        '''
