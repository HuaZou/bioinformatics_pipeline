rule sort_by_name:
    input:
        expand("{alignment}/{{sample}}.bam",
                alignment=config["result"]["alignment"])
    output:
        expand("{alignment}/{{sample}}.sortedByName.bam",
                alignment=config["result"]["alignment"])
    threads: 10
    log:
        os.path.join(config["logs"]["alignment"], "{sample}.SortByName.log")
    shell:
        '''
        samtools sort -@ {threads} -on {input} -T /tmp/ -o {output} 2>{log}
        '''

rule htseq_count:
    input:
        expand("{alignment}/{{sample}}.sortedByName.bam",
                alignment=config["result"]["alignment"])
    output: 
        expand("{genecount}/{{sample}}.counts.tsv",
                genecount=config["result"]["genecount"])
    params:
        annnotation  = config["params"]["genecount"]["gtf"],
        phred_cutoff = config["params"]["genecount"]["phred_cutoff"]
    shell:
        '''
        conda activate py27 && htseq-count --order=name --format=bam --mode=intersection-strict --stranded=no --minaqual={params.phred_cutoff} --type=exon --idattr=gene_id {input} {params.annotation} > {output} && conda deactivate
        '''
