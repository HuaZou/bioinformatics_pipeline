rule hisat2_alignment:
    input:
        read1 = expand("{trim}/{{sample}}_1.fastq.gz",
                        trim=config["result"]["trim"]),
        read2 = expand("{trim}/{{sample}}_2.fastq.gz",
                        trim=config["result"]["trim"])
    output:
        sam = temp(expand("{alignment}/{{sample}}.sam",
                    alignment=config["result"]["alignment"])),
        bam = temp(expand("{alignment}/{{sample}}.bam",
                       alignment=config["result"]["alignment"])),
        sortbam = expand("{alignment}/{{sample}}.sorted.bam",
                       alignment=config["result"]["alignment"])
    params:
        prefix = "{sample}",
        index  = config["params"]["alignment"]["index"],
        sortbam_prefix = os.path.join(config["result"]["alignment"], "{sample}.sorted"),
        cpu_hisat2     = config["params"]["alignment"]["cpu_hisat2"],
        cpu_samtools   = config["params"]["alignment"]["cpu_samtools"]
    log:
        os.path.join(config["logs"]["alignment"], "alignment.{sample}.log")
    shell:
        '''
        hisat2 -x {params.index} \
            -p {params.cpu_hisat2} \
            -1 {input.read1} \
            -2 {input.read2} \
            -S {output.sam} \
            --rg-id={params.prefix} 2>{log}
        samtools view -S -b {output.sam} > {output.bam} 2>>{log}
        samtools sort -@ {params.cpu_samtools} {output.bam} {params.sortbam_prefix} 2>>{log}
        samtools index {output.sortbam} 2>>{log}
        '''
