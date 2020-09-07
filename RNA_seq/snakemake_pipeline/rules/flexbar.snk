def get_fastq(wildcards, samples, read_pair="fq1"):
    return samples.loc[wildcards.sample, [read_pair]]

rule flexbar_trim:
    input:
        read1 = lambda wildcards: get_fastq(wildcards, samples, "fq1"),
        read2 = lambda wildcards: get_fastq(wildcards, samples, "fq2")
    output:
        outfile = expand("{trim}/{{sample}}_{read}.fastq.gz",
                        trim=config["result"]["trim"],
                        read=["1", "2"])
    params:
        prefix  = os.path.join(config["result"]["trim"], "{sample}"),
        ao      = config["params"]["trim"]["ao"],
        quality = config["params"]["trim"]["quality"],
        min_length  = config["params"]["trim"]["min_length"],
        threads     = config["params"]["trim"]["cpu"],
        max_uncall  = config["params"]["trim"]["max_uncalled"],
        adapter     = config["params"]["trim"]["adapter"]
    log:
        os.path.join(config["logs"]["trim"], "trim.{sample}.log")
    shell:
        '''
        flexbar -a {params.adapter}\
            -ao {params.ao} -qt {params.quality} \
            --min-read-length {params.min_length} \
            --threads {params.threads} \
            --zip-output GZ \
            --max-uncalled {params.max_uncall} \
            -r {input.read1} \
            -p {input.read2}  \
            --target {params.prefix}
        '''
