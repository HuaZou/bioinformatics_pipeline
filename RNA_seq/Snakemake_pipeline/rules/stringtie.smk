rule stringtie_tsv:
    input:
        bam = expand("{alignment}/{{sample}}.sorted.bam",
                        alignment=config["result"]["alignment"])
    output:
        gtf = expand("{genecount}/{{sample}}.gtf",
                        genecount=config["result"]["genecount"]),
        tsv = expand("{genecount}/{{sample}}.tsv",
                        genecount=config["result"]["genecount"])
    params:
        annnote  = config["params"]["genecount"]["gtf"],
        cpu      = config["params"]["genecount"]["cpu"]
    log:
        os.path.join(config["logs"]["genecount"], "genecount.{sample}.log")
    shell:
        '''
        stringtie --rf -p {params.cpu} -G {params.annnote} -e -B \
            -o {output.gtf} -A {output.tsv} {input.bam} 2>{log}
        '''
"""
rule stringtie_count:
    input:
        gtf = expand("{genecount}/{sampleID}.gtf",
                        genecount=config["result"]["genecount"],
                        sampleID=samples.index)
    output:
        gtf_list   = temp(os.path.join(config["result"]["genecount"], "stringtie.gtf.lst")),
        gene       = os.path.join(config["result"]["genecount"], "stringtie_gene_counts.csv"),
        transcript = os.path.join(config["result"]["genecount"], "stringtie_transcript_counts.csv")
    params:
        count_dir    = config["result"]["genecount"]
    log:
        os.path.join(config["logs"]["genecount"], "tringtie_count.log")
    shell:
        '''
        find {params.count_dir} -name "*.gtf" | perl -ne 'chomp; $name=(split("\/", $_))[-1]; $name=~s/\.gtf//g; print "$name\t$_\n"' > {output.gtf_list}
        conda activate py27 && python /disk/user/zouhua/script/preDE.py -i {output.gtf_list} -g {output.gene} -t {output.transcript} 2>{log} && conda deactivate
        '''
"""