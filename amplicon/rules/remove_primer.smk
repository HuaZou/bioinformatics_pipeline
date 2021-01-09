rule remove_primer_pe:
    input:
        os.path.join(config["result"]["import"], "paired-end-demux.qza")
    params:
        forward_seq = config["params"]["remove"]["forward_sequence_515F"],
        reverse_seq = config["params"]["remove"]["reverse_sequence_806R"]
    threads: 
        config["params"]["remove"]["remove_primer_threads"]
    log:
        os.path.join(config["logs"]["remove"], "remove_primer.log")
    output:
        os.path.join(config["result"]["remove"], "paired-end-demux-trim.qza")
    shell:
        '''
        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences {input} \
            --p-front-f {params.forward_seq} \
            --p-front-r {params.reverse_seq} \
            --p-cores {threads} \
            --o-trimmed-sequences {output} 2>{log}
        '''
