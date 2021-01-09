rule alignment_mafft:
    input:
        os.path.join(config["result"]["denoise"], "repseqs_final.qza")
    threads: 
        config["params"]["diversity"]["alignment_mafft_threads"]
    log:
        os.path.join(config["logs"]["tree"], "alignment_mafft.log")
    output:
        os.path.join(config["result"]["tree"], "unmasked_aligned_repseqs.qza")
    shell:
        '''
        qiime alignment mafft \
            --i-sequences {input} \
            --o-alignment {output} \
            --p-n-threads {threads} 2>{log}
        '''

rule alignment_mask:
    input:
        os.path.join(config["result"]["tree"], "unmasked_aligned_repseqs.qza")
    log:
        os.path.join(config["logs"]["tree"], "alignment_mask.log")
    output:
        os.path.join(config["result"]["tree"], "aligned_repseqs.qza")
    shell:
        '''
        qiime alignment mask \
            --i-alignment {input} \
            --o-masked-alignment {output} 2>{log}
        '''

rule phylogeny_fasttree:
    input:
        os.path.join(config["result"]["tree"], "aligned_repseqs.qza")
    threads: 
        config["params"]["diversity"]["phylogeny_fasttree_threads"]
    log:
        os.path.join(config["logs"]["tree"], "phylogeny_fasttree.log")
    output:
        os.path.join(config["result"]["tree"], "unrooted_tree.qza")
    shell:
        '''
        qiime phylogeny fasttree \
            --i-alignment {input} \
            --o-tree {output} \
            --p-n-threads {threads} 2>{log}
        '''

rule phylogeny_midpoint_root:
    input:
        os.path.join(config["result"]["tree"], "unrooted_tree.qza")
    output:
        os.path.join(config["result"]["tree"], "rooted_tree.qza")
    shell:
        '''
        qiime phylogeny midpoint-root \
            --i-tree {input} \
            --o-rooted-tree {output}
        '''

rule visualize_tree:
    input:
        tree     = os.path.join(config["result"]["tree"], "rooted_tree.qza"),
        table    = os.path.join(config["result"]["denoise"], "table_final.qza"),
        metadata = config["metadata"],
        taxonomy = os.path.join(config["result"]["taxonomy"], "taxonomy.qza"),
        outliers = os.path.join(config["result"]["tree"], "outliers.qza")
    log:
        os.path.join(config["logs"]["tree"], "visualize_tree.log")
    output:
         os.path.join(config["result"]["tree"], "rooted_tree.qzv")
    shell:
        '''
        qiime empress community-plot \
            --i-tree {input.tree} \
            --i-feature-table {input.table} \
            --m-sample-metadata-file {input.metadata} \
            --m-feature-metadata-file {input.taxonomy} \
            --m-feature-metadata-file {input.outliers} \
            --o-visualization {output} 2>{log}
        '''

rule unzip_alignment_to_fasta:
    input:
        os.path.join(config["result"]["tree"], "aligned_repseqs.qza")
    output:
        os.path.join(config["result"]["tree"], "aligned_repseqs.fasta")
    shell:
        '''
        unzip -qq -o {input} -d temp; 
        mv temp/*/data/aligned-dna-sequences.fasta {output}; 
        rm -r temp
        '''

rule alignment_count_gaps:
    input:
        os.path.join(config["result"]["tree"], "aligned_repseqs.fasta")
    output:
        os.path.join(config["result"]["tree"], "aligned_repseqs_gaps.txt")
    shell:
        '''
        cat {input} | grep -v '>' | sed 's/[^-]//g' | awk '{{ print length }}' > {output}
        '''

rule alignment_gaps_describe:
    input:
        os.path.join(config["result"]["tree"], "aligned_repseqs_gaps.txt")
    output:
        os.path.join(config["result"]["tree"], "aligned_repseqs_gaps_describe.md")
    run:
        gaps = pd.read_csv(input[0], header=None)
        t = gaps.describe()
        outstr = tabulate(t.iloc[1:], tablefmt="pipe", headers=['Statistic (n=%s)' % t.iloc[0].values[0].astype(int), 'Alignment gaps per sequence'])
        with open(output[0], 'w') as target:
            target.write(outstr)
            target.write('\n')

rule alignment_detect_outliers:
    input:
        os.path.join(config["result"]["tree"], "aligned_repseqs.fasta")
    params:
        metric     = config["params"]["diversity"]["odseq_distance_metric"],
        replicates = config["params"]["diversity"]["odseq_bootstrap_replicates"],
        threshold  = config["params"]["diversity"]["odseq_threshold"]
    output:
        os.path.join(config["result"]["tree"], "aligned_repseqs_outliers.tsv")
    shell:
        "Rscript --vanilla scripts/run_odseq.R {input} {params.metric} {params.replicates} {params.threshold} {output}"

rule tabulate_plot_repseq_properties:
    input:
        lengths  = os.path.join(config["result"]["denoise"], "repseqs_lengths.txt"),                 
        gaps     = os.path.join(config["result"]["tree"], "aligned_repseqs_gaps.txt"),
        outliers = os.path.join(config["result"]["tree"], "aligned_repseqs_outliers.tsv"),
        taxonomy = os.path.join(config["result"]["taxonomy"], "taxonomy.qza"),
        table    = os.path.join(config["result"]["denoise"], "table.qza")
    output:
        proptsv  = os.path.join(config["result"]["tree"], "repseqs_properties.tsv"),
        proppdf  = os.path.join(config["result"]["tree"], "repseqs_properties.pdf"),
        outliersforqza = os.path.join(config["result"]["tree"], "outliers.tsv")
    run:
        lengths = pd.read_csv(input['lengths'], header=None)
        gaps = pd.read_csv(input['gaps'], header=None)
        outliers = pd.read_csv(input['outliers'], header=None, sep='\t')
        taxonomy = Artifact.load(input['taxonomy'])
        taxonomydf = taxonomy.view(view_type=pd.DataFrame)
        taxonomydf['level_1'] = [x.split(';')[0] for x in taxonomydf['Taxon']]
        table = Artifact.load(input['table'])
        tabledf = table.view(view_type=pd.DataFrame)
        merged = pd.concat([lengths, gaps, outliers[1], taxonomydf['Taxon'].reset_index(drop=True), taxonomydf['level_1'].reset_index(drop=True), tabledf.sum().reset_index(drop=True)],
                           axis=1, ignore_index=True)
        merged.columns = ['featureid', 'length', 'gaps', 'outlier', 'taxonomy', 'taxonomy_level_1', 'observations']
        merged['log10(observations)'] = [np.log10(x) for x in merged['observations']]
        merged.sort_values('log10(observations)', ascending=False, inplace=True)
        merged.to_csv(output['proptsv'], index=False, sep='\t')
        g = sns.relplot(data=merged, x='length', y='gaps', col='outlier', hue='taxonomy_level_1', size='log10(observations)', sizes=(1,500), edgecolor = 'none', alpha=0.7)
        g.set_axis_labels('length (bp) not including gaps', 'gaps (bp) in masked multiple sequence alignment')
        plt.savefig(output['proppdf'], bbox_inches='tight')
        outliers.columns = ['Feature ID', 'Outlier']
        outliers = outliers*1
        outliers.to_csv(output['outliersforqza'], index=False, sep='\t')

rule import_outliers_to_qza:
    input:
        outliersforqza = os.path.join(config["result"]["tree"], "outliers.tsv")
    output:
        outliersforqza = os.path.join(config["result"]["tree"], "outliers.qza")
    shell:
        '''
        qiime tools import \
            --type 'FeatureData[Importance]' \
            --input-path {input} \
            --output-path {output}
        '''

rule tabulate_repseqs_to_filter:
    input:
        proptsv  = os.path.join(config["result"]["tree"], "repseqs_properties.tsv")
    output:
        outliers   = os.path.join(config["result"]["tree"], "repseqs_to_filter_outliers.tsv"),
        unassigned = os.path.join(config["result"]["tree"], "repseqs_to_filter_unassigned.tsv")
    shell:
        "cat {input.proptsv} | grep -i 'outlier\|true' | cut -f1,4 > {output.outliers}; "
        "cat {input.proptsv} | grep -i 'taxonomy\|unassigned' | cut -f1,5 > {output.unassigned}"
        