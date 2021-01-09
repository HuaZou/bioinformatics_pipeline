rule summarize_feature_table:
    input:
        table = os.path.join(config["result"]["denoise"], "table_final.qza"),
        metadata = config["metadata"]
    output:
        os.path.join(config["result"]["denoise"], "table_final.qzv")
    shell:
        '''
        qiime feature-table summarize \
            --i-table {input.table} \
            --m-sample-metadata-file {input.metadata} \
            --o-visualization {output}
        '''

rule unzip_table_to_biom:
    input:
        os.path.join(config["result"]["denoise"], "table_final.qza")
    output:
        os.path.join(config["result"]["denoise"], "table.biom")
    shell:
        '''
        unzip -qq -o {input} -d temp;
        mv temp/*/data/feature-table.biom {output};
        rm -r temp
        '''

rule summarize_biom_samples:
    input:
        os.path.join(config["result"]["denoise"], "table.biom")
    output:
        os.path.join(config["result"]["denoise"], "table_summary_samples.txt")
    shell:
        '''
        biom summarize-table \
            --input-fp {input} \
            --output-fp {output}; 
        cat {output} | sed 's/observation/feature/g' | sed 's/.000$//' > temp; 
        mv temp {output}
        '''

rule summarize_biom_features:
    input:
        os.path.join(config["result"]["denoise"], "table.biom")
    output:
        os.path.join(config["result"]["denoise"], "table_summary_features.txt")
    shell:
        '''
        biom summarize-table \
            --observations \
            --input-fp {input} \
            --output-fp {output}; 
        cat {output} | sed 's/observation/feature/g' | sed 's/Counts\/sample/Counts\/feature/g' | sed 's/.000$//' > temp; 
        mv temp {output}
        '''

rule visualize_repseqs:
    input:
        os.path.join(config["result"]["denoise"], "repseqs_final.qza")
    output:
        os.path.join(config["result"]["denoise"], "repseqs_final.qzv")
    shell:
        '''
        qiime feature-table tabulate-seqs \
            --i-data {input} \
            --o-visualization {output}
        '''

rule unzip_repseqs_to_fasta:
    input:
        os.path.join(config["result"]["denoise"], "repseqs_final.qza")
    output:
        os.path.join(config["result"]["denoise"], "repseqs_final.fasta")
    shell:
        '''
        unzip -qq -o {input} -d temp; 
        mv temp/*/data/dna-sequences.fasta {output}; 
        rm -r temp
        '''

rule repseqs_detect_amplicon_locus:
    input:
        os.path.join(config["result"]["denoise"], "repseqs_final.fasta")
    output:
        os.path.join(config["result"]["denoise"], "repseqs_amplicon_type.txt")
    shell:
        '''
        python scripts/detect_amplicon_locus.py -i {input} > {output}
        '''

rule repseqs_lengths:
    input:
        os.path.join(config["result"]["denoise"], "repseqs_final.fasta")
    output:
        os.path.join(config["result"]["denoise"], "repseqs_lengths.txt")
    shell:
        '''
        perl scripts/fastaLengths.pl {input} > {output}
        '''

rule repseqs_lengths_describe:
    input:
        os.path.join(config["result"]["denoise"], "repseqs_lengths.txt")
    output:
        os.path.join(config["result"]["denoise"], "repseqs_lengths_describe.md")
    run:
        s = pd.read_csv(input[0], header=None, index_col=0)
        t = s.describe()
        outstr = tabulate(t.iloc[1:], tablefmt="pipe", headers=['Statistic (n=%s)' % t.iloc[0].values[0].astype(int), 'Sequence length'])
        with open(output[0], 'w') as target:
            target.write(outstr)
            target.write('\n')
