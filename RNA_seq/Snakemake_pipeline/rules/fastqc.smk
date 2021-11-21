def get_fastq(wildcards, samples, read_pair="fq1"):
    return samples.loc[wildcards.sample, [read_pair]]

rule fastqc:
    input:
        read1 = lambda wildcards: get_fastq(wildcards, samples, "fq1"),
        read2 = lambda wildcards: get_fastq(wildcards, samples, "fq2")
    output:
        outfile = expand("{fastqc}/{{sample}}.{read}_fastqc.{out}",
                        fastqc=config["result"]["fastqc"],
                        read=["r1", "r2"],
                        out=["html", "zip"])
    params:
        outdir = config["result"]["fastqc"]
    log:
        os.path.join(config["logs"]["fastqc"], "fastqc.{sample}.log")
    shell:
        '''
        fastqc --noextract -f fastq {input.read1} {input.read2} -o {params.outdir} 2>{log}
        '''

rule multiqc:
    input:
        expand("{fastqc}/{sample}.{read}_fastqc.zip",
                fastqc=config["result"]["fastqc"],
                sample=samples.index,
                read=["r1", "r2"])
    output:
        html = os.path.join(config["result"]["fastqc"], "fastqc_multiqc_report.html"),
        data_dir = directory(os.path.join(config["result"]["fastqc"], "fastqc_multiqc_report_data"))
    params:
        outdir = config["result"]["fastqc"]
    log:
        os.path.join(config["logs"]["fastqc"], "fastqc__multiqc.log")
    shell:
        '''
        multiqc --title fastqc --module fastqc {input} --outdir {params.outdir}  2> {log}
        '''
