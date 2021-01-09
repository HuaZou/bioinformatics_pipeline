rule remove_chimera:
    input:
        repseqs = os.path.join(config["result"]["denoise"], "repseqs.qza"),
        table   = os.path.join(config["result"]["denoise"], "table.qza")
    log:
        os.path.join(config["logs"]["denoise"], "remove_chimera.log")
    output:
        chimeras_seq    = os.path.join(config["result"]["denoise"], "chimeras.qza"),
        nonchimeras_seq = os.path.join(config["result"]["denoise"], "nonchimera.qza"),
        chimeras_stats  = os.path.join(config["result"]["denoise"], "chimeras_stats.qza")
    shell:
        '''
        qiime vsearch uchime-denovo \
            --i-sequences {input.repseqs} \
            --i-table {input.table} \
            --o-chimeras {output.chimeras_seq} \
            --o-nonchimeras {output.nonchimeras_seq} \
            --o-stats {output.chimeras_stats} \
            --verbose 2>{log}
        '''

rule filter_table_req_by_chimera:
    input:
        repseqs = os.path.join(config["result"]["denoise"], "repseqs.qza"),
        table   = os.path.join(config["result"]["denoise"], "table.qza"),
        nonchimeras_seq = os.path.join(config["result"]["denoise"], "nonchimera.qza")
    output:
        filteredseqs = os.path.join(config["result"]["denoise"], "repseqs_nonchimera.qza"),
        filteredtable  = os.path.join(config["result"]["denoise"], "table_nonchimera.qza")
    shell:
        '''
        qiime feature-table filter-seqs \
            --i-data {input.repseqs} \
            --m-metadata-file {input.nonchimeras_seq} \
            --o-filtered-data {output.filteredseqs};
        qiime feature-table filter-features \
            --i-table {input.table} \
            --m-metadata-file {input.nonchimeras_seq} \
            --o-filtered-table {output.filteredtable}        
        '''

rule filter_sequences_by_taxonomy:
    input:
        repseqs  = os.path.join(config["result"]["denoise"], "repseqs_nonchimera.qza"),
        taxonomy = os.path.join(config["result"]["taxonomy"], "taxonomy.qza")
    params:
        excludeterms = config["params"]["filter_repseq"]["exclude_terms"]
    log:
        os.path.join(config["logs"]["denoise"], "remove_contamination.log")
    output:
        os.path.join(config["result"]["denoise"], "repseqs_nonchimera_contam.qza")
    shell:
        '''
        qiime taxa filter-seqs \
            --i-sequences {input.repseqs} \
            --i-taxonomy {input.taxonomy} \
            --p-exclude {params.excludeterms} \
            --o-filtered-sequences {output} 2>{log}
        '''

"""
rule filter_table_req_by_contamination:
    input:
        repseqs = os.path.join(config["result"]["denoise"], "repseqs_nonchimera.qza"),
        table   = os.path.join(config["result"]["denoise"], "table_nonchimera.qza"),
        noncontam_seq = os.path.join(config["result"]["denoise"], "repseqs_nonchimera_contam.qza")
    output:
        filteredseqs = os.path.join(config["result"]["denoise"], "repseq_final.qza"),
        filteredtable  = os.path.join(config["result"]["denoise"], "table_final.qza")
    shell:
        '''
        qiime feature-table filter-seqs \
            --i-data {input.repseqs} \
            --m-metadata-file {input.noncontam_seq} \
            --o-filtered-data {output.filteredseqs};
        qiime feature-table filter-features \
            --i-table {input.table} \
            --m-metadata-file {input.noncontam_seq} \
            --o-filtered-data {output.filteredtable}        
        '''
"""

rule filter_sequences_by_id:
    input:
        repseqs         = os.path.join(config["result"]["denoise"], "repseqs_nonchimera_contam.qza"),
        repseqstofilter = config["params"]["filter_repseq"]["repseqs_to_filter"]
    output:
        os.path.join(config["result"]["denoise"], "repseq_final.qza")
    shell:
        '''
        qiime feature-table filter-seqs \
            --i-data {input.repseqs} \
            --m-metadata-file {input.repseqstofilter} \
            --p-exclude-ids \
            --o-filtered-data {output}
        '''

rule filter_table_by_contam:
    input:
        table        = os.path.join(config["result"]["denoise"], "table_nonchimera.qza"),
        filteredseqs = os.path.join(config["result"]["denoise"], "repseq_final.qza")
    output:
        os.path.join(config["result"]["denoise"], "table_final.qza")
    shell:
        '''
        qiime feature-table filter-features \
            --i-table {input.table} \
            --m-metadata-file {input.filteredseqs} \
            --o-filtered-table {output}
        '''

rule filter_taxonomy_by_contam:
    input:
        taxonomy        = os.path.join(config["result"]["taxonomy"], "taxonomy.tsv"),
        repseqstofilter = config["params"]["filter_repseq"]["repseqs_to_filter"]
    params:
        excludeterms = config["params"]["filter_repseq"]["exclude_terms"]
    output:
        taxonomy = os.path.join(config["result"]["taxonomy"], "taxonomy_final.tsv")
    run:
        df_taxonomy = pd.read_csv(input['taxonomy'], sep='\t', index_col=0)
        exclude_terms = params['excludeterms'].split(',')
        exclude_ids_taxa = df_taxonomy.index[[any(x.lower() in y.lower() for x in exclude_terms) for y in df_taxonomy['Taxon']]]
        df_repseqs_to_filter = pd.read_csv(input['repseqstofilter'], sep='\t')
        exclude_ids_direct = df_repseqs_to_filter['featureid'].values
        keep_ids = set(df_taxonomy.index) - set(exclude_ids_taxa) - set(exclude_ids_direct)
        df_taxonomy_filtered = df_taxonomy.loc[list(keep_ids)]
        df_taxonomy_filtered.to_csv(output['taxonomy'], sep='\t')
