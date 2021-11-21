rule star_alignment:
    input:
        read1 = expand("{trim}/{{sample}}.r1_val_1.fq.gz",
                        trim=config["result"]["trim"]),
        read2 = expand("{trim}/{{sample}}.r2_val_2.fq.gz",
                        trim=config["result"]["trim"])
    output:
        bam = expand("{alignment}/{{sample}}.bam",
                       alignment=config["result"]["alignment"])
    params:
        prefix   = os.path.join(config["result"]["alignment"], "{sample}"),
        starlogs = os.path.join(config["result"]["alignment"], "starlogs"),
        index    = config["params"]["alignment"]["index"],
        threads  = config["params"]["alignment"]["cpu"],
        matchNmin   = config["params"]["alignment"]["matchNmin"],
        misNmax   = config["params"]["alignment"]["misNmax"]
    log:
        os.path.join(config["logs"]["alignment"], "alignment.{sample}.log")
    shell:
        '''
        STAR --runThreadN {params.threads}\
            --genomeDir {params.index}\
            --outFileNamePrefix {params.prefix} --readFilesIn {input.read1} {input.read2}\
            --outSAMtype BAM SortedByCoordinate\
            --outFilterMatchNmin {params.matchNmin}\
            --outFilterMismatchNmax {params.misNmax}\
            --quantMode GeneCounts \
            --readFilesCommand zcat && mv {params.prefix}Aligned.sortedByCoord.out.bam {output} && mkdir -p {params.starlogs} && mv {params.prefix}Log.final.out {params.prefix}Log.out {params.prefix}Log.progress.out {params.starlogs} 2>{log} 
        '''
