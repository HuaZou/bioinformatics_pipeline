rule summarize_fastq_demux_pe:
    input:
        qza = os.path.join(config["result"]["import"], "paired-end-demux.qza"),
        trim_qza = os.path.join(config["result"]["remove"], "paired-end-demux-trim.qza")
    output:
        fastq_summary = os.path.join(config["result"]["import"], "paired-end-demux-fastq_summary.qzv"),
        trim_fastq_summary = os.path.join(config["result"]["remove"], "paired-end-demux-trim-fastq_summary.qzv")
    shell:
        '''
        qiime demux summarize \
            --i-data {input.qza} \
            --o-visualization {output.fastq_summary}; 
        qiime demux summarize \
            --i-data {input.trim_qza} \
            --o-visualization {output.trim_fastq_summary}      
        '''

rule unzip_fastq_summary:
    input:
        fastq_summary = os.path.join(config["result"]["import"], "paired-end-demux-fastq_summary.qzv"),
        trim_fastq_summary = os.path.join(config["result"]["remove"], "paired-end-demux-trim-fastq_summary.qzv")
    output:
        fastq_summary_tsv = os.path.join(config["result"]["import"], "paired-end-demux-fastq_counts.tsv"),
        trim_fastq_summary_tsv = os.path.join(config["result"]["remove"], "paired-end-demux-trim-fastq_counts.tsv")
    shell:
        '''
        unzip -qq -o {input.fastq_summary} -d temp;
        mv temp/*/data/per-sample-fastq-counts.tsv {output.fastq_summary_tsv};
        rm -r temp;
        unzip -qq -o {input.trim_fastq_summary} -d temp;
        mv temp/*/data/per-sample-fastq-counts.tsv {output.trim_fastq_summary_tsv};
        rm -r temp
        '''

rule describe_fastq_counts:
    input:
        fastq_summary_tsv = os.path.join(config["result"]["import"], "paired-end-demux-fastq_counts.tsv"),
        trim_fastq_summary_tsv = os.path.join(config["result"]["remove"], "paired-end-demux-trim-fastq_counts.tsv")
    output:
        fastq_describe_md = os.path.join(config["result"]["import"], "paired-end-demux-fastq_counts_describe.md"),
        trim_fastq_describe_md = os.path.join(config["result"]["remove"], "paired-end-demux-trim-fastq_counts_describe.md")
    run:
        def get_describe(infile, outfile):
            s = pd.read_csv(infile, sep='\t', index_col=0)
            t = s.describe()
            outstr = tabulate(pd.DataFrame(t.iloc[1:,0]), tablefmt="pipe", headers=['Statistic (n=%s)' % t.iloc[0,0].astype(int), 'Fastq sequences per sample'])
            with open(outfile, 'w') as target:
                target.write(outstr)
                target.write('\n')
        
        get_describe(input[0], output[0])
        get_describe(input[1], output[1])

        