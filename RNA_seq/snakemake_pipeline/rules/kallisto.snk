rule kallisto_quantify:
    input:
        read1 = expand("{trim}/{{sample}}_1.fastq.gz",
                        trim=config["result"]["trim"]),
        read2 = expand("{trim}/{{sample}}_2.fastq.gz",
                        trim=config["result"]["trim"])
    output:
        tsv = expand("{quantify}/{{sample}}/abundance.tsv",
                       quantify=config["result"]["quantify"]["kallisto"])
    params:
        kallisto_index   = config["params"]["quantify"]["kallisto_index"],
        gtf     = config["params"]["quantify"]["gtf"],
        threads = config["params"]["quantify"]["threads"],
        prefix = os.path.join(config["result"]["quantify"]["kallisto"], "{sample}")
    log:
        os.path.join(config["logs"]["quantify"], "kallisto.quantify.{sample}.log")
    shell:
        '''
        kallisto quant -i {params.kallisto_index} \
        -g {params.gtf} \
        -t {params.threads} \
        -o {params.prefix} \
        -b 30 \
        -l 200 \
        -s 30 \
        {input.read1} {input.read2} 2>{log}
        '''


rule salmon_quantify:
    input:
        read1 = expand("{trim}/{{sample}}_1.fastq.gz",
                        trim=config["result"]["trim"]),
        read2 = expand("{trim}/{{sample}}_2.fastq.gz",
                        trim=config["result"]["trim"])
    output:
        sf = expand("{quantify}/{{sample}}/abundance.sf",
                       quantify=config["result"]["quantify"]["salmon"])
    params:
        salmon_index   = config["params"]["quantify"]["salmon_index"],
        gtf     = config["params"]["quantify"]["gtf"],
        threads = config["params"]["quantify"]["threads"],
        numBootstraps = config["params"]["quantify"]["numBootstraps"],
        prefix = os.path.join(config["result"]["quantify"]["salmon"], "{sample}")
    log:
        os.path.join(config["logs"]["quantify"], "salmon.quantify.{sample}.log")
    shell:
        '''
        salmon quant -i {params.salmon_index} \
        -g {params.gtf} \
        -p {params.threads} \
        -o {params.prefix} \
        -l IU -p 10 \ 
        -1 {input.read1} -2 {input.read2} \
        --numBootstraps {params.numBootstraps} 2>{log}
        '''
