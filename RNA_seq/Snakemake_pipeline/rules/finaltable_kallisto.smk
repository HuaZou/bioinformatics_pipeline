rule merge_table:
    input:
        expand("{genecount}/{sampleID}.tsv",
                        genecount=config["result"]["genecount"],
                        sampleID=samples.index)
    output:
        files = expand("{finaltable}/stringtie.{type}.{out}",
                        finaltable=config["result"]["finaltable"],
                        type=["FPKM", "TPM"],
                        out=["tsv", "cover", "count"])
    params:
        count_dir    = config["result"]["genecount"],
        prefix_FPKM  = os.path.join(config["result"]["finaltable"], "stringtie.FPKM"),
        prefix_TPM   = os.path.join(config["result"]["finaltable"], "stringtie.TPM"),
        FPKM    = config["params"]["finaltable"]["FPKM"],
        TPM     = config["params"]["finaltable"]["TPM"],
        script  = config["params"]["finaltable"]["script"]
    log:
        os.path.join(config["logs"]["finaltable"], "finaltable.log")
    shell:
        '''
        {params.script} -d {params.count_dir} -r {params.FPKM} -a 1 -b 1 -o {params.prefix_FPKM} 2>{log}
        {params.script} -d {params.count_dir} -r {params.TPM} -a 1 -b 1 -o {params.prefix_TPM} 2>>{log}
        '''

rule get_gtfpath:
    input:
        expand("{genecount}/{sampleID}.gtf",
                        genecount=config["result"]["genecount"],
                        sampleID=samples.index)
    output:
        os.path.join(config["result"]["finaltable"], "stringtie.gtf.lst"),
    run:
        import re
        import os
        import sys

        out_f = open(output[0], "w")
        for gtf_files in input:
            group = re.findall('(\S+)\.gtf', gtf_files)
            name_path = "\t".join([group[0], gtf_files])
            out_f.write(name_path + "\n")
        out_f.close()

rule stringtie_count:
    input:
        gtf_list = os.path.join(config["result"]["finaltable"], "stringtie.gtf.lst")
    output:
        gene       = os.path.join(config["result"]["finaltable"], "stringtie.GeneCount.csv"),
        transcript = os.path.join(config["result"]["finaltable"], "stringtie.TranscriptCount.csv")
    params:
        script     = config["params"]["finaltable"]["script_count"]
    log:
        os.path.join(config["logs"]["finaltable"], "stringtie_count.log") 
    shell:
        '''
        source activate py27 && python {params.script} -i {input.gtf_list} -g {output.gene} -t {output.transcript} && conda deactivate 2>{log}
        '''
        