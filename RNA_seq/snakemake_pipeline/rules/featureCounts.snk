rule featureCounts_tsv:
    input:
        expand("{alignment}/{sample}.bam",
                alignment=config["result"]["alignment"],
                sample=samples.index)
    output:
        os.path.join(config["result"]["genecount"], "final_count.tsv")
    params:
        annnote  = config["params"]["genecount"]["gtf"],
        cpu      = config["params"]["genecount"]["cpu"],
        gene     = config["params"]["genecount"]["gene"],
        type     = config["params"]["genecount"]["type"]
    log:
        os.path.join(config["logs"]["genecount"], "genecount.featurecounts.log")
    shell:
        '''
        featureCounts \
            -a {params.annnote} \
            -o {output} \
            -t {params.type} \
            -g {params.gene} \
            --primary \
            -T {params.cpu} \
            {input} 2>{log}
        '''
