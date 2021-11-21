rule multiqc_result:
    input:
        files = expand("{finaltable}/stringtie.{type}.{out}",
                        finaltable=config["result"]["finaltable"],
                        type=["FPKM", "TPM"],
                        out=["tsv", "cover", "count"]), 
        files2 = expand("{finaltable}/stringtie.{type2}.csv",
                        finaltable=config["result"]["finaltable"],
                        type2=["GeneCount", "TranscriptCount"]),
        html = os.path.join(config["result"]["fastqc"], "fastqc_multiqc_report.html")
    output:
        os.path.join(config["result"]["multiqc"], "multiqc_report.html")
    params:
        fqc  = config["result"]["fastqc"],
        trim = config["result"]["trim"],
        alignment  = config["result"]["alignment"],
        genecount  = config["result"]["genecount"],
        finaltable = config["result"]["finaltable"],
        dir  = config["result"]["multiqc"]
    log:
        os.path.join(config["logs"]["multiqc"], "multiqc.log")
    shell:
        '''
        multiqc {params.fqc} {params.trim} {params.alignment} {params.genecount} {params.finaltable} --outdir {params.dir} 2>{log}
        '''
