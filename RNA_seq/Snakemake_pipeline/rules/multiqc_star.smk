rule multiqc_result:
    input:
        datatsv = os.path.join(config["result"]["genecount"], "final_count.tsv"),
        html = os.path.join(config["result"]["fastqc"], "fastqc_multiqc_report.html")
    output:
        os.path.join(config["result"]["multiqc"], "multiqc_report.html")
    params:
        fqc  = config["result"]["fastqc"],
        trim = config["result"]["trim"],
        alignment  = config["result"]["alignment"],
        genecount  = config["result"]["genecount"],
        dir  = config["result"]["multiqc"]
    log:
        os.path.join(config["logs"]["multiqc"], "multiqc.log")
    shell:
        '''
        multiqc {params.fqc} {params.trim} {params.alignment} {params.genecount} --outdir {params.dir} 2>{log}
        '''
