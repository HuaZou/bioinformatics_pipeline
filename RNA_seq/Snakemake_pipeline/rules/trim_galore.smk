def get_fastq(wildcards, samples, read_pair="fq1"):
    return samples.loc[wildcards.sample, [read_pair]]

rule trim_galore:
    input:
        read1 = lambda wildcards: get_fastq(wildcards, samples, "fq1"),
        read2 = lambda wildcards: get_fastq(wildcards, samples, "fq2")
    output:
        R1 = os.path.join(config["result"]["trim"], "{sample}.r1_val_1.fq.gz"),
        R2 = os.path.join(config["result"]["trim"], "{sample}.r2_val_2.fq.gz")
    params:
        quality     = config["params"]["trim"]["quality"],
        length      = config["params"]["trim"]["length"],
        error       = config["params"]["trim"]["error"],
        stringency  = config["params"]["trim"]["stringency"],
        adapter1    = config["params"]["trim"]["adapter1"],
        adapter2    = config["params"]["trim"]["adapter2"],
        outdir      = config["result"]["trim"]
    log:
        os.path.join(config["logs"]["trim"], "trim.{sample}.log")
    shell:
        '''
        trim_galore --paired -o {params.outdir}\
            {input.read1} {input.read2} \
            --quality {params.quality} \
            --length {params.length} \
            -e {params.error} \
            --stringency {params.stringency} \
            -a {params.adapter1} \
            -a2 {params.adapter2} \
            --fastqc
        '''
