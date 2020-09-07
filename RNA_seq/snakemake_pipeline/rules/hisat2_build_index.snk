rules hisat2_build:
    input:
        gtf = config[params.gtf],
        fa  = config[params.fa]
    output:
        exon   = temp(expand({index}/Mus_musculus.GRCm38.100.exons.gtf,
                    index=config["result"]["index"])),
        splice = temp(expand({index}/Mus_musculus.GRCm38.100.splicesites.gtf,
                    index=config["result"]["index"])),                    
        index  = expand("{index}/Mus_musculus.GRCm38.100.exon_splice.{number}.ht2l",
                        index=config["result"]["index"],
                        number=["1", "2"])
    params:
        prefix = os.path.join(config["result"]["index"], "Mus_musculus.GRCm38.100.exon_splice"),
        cpu    = config["parsms"]["index"]["cpu"]
    log:
        os.path.join(config["logs"]["index"], "index.log")
    
    shell:
        '''
        hisat2_extract_exons.py {input.gtf} > {output.exon}
        hisat2_extract_splice_sites.py {input.gtf} > {output.splice}
        hisat2-build -p 30 --ss {output.splice} --exon {output.exon} {input.fa} {params.prefix} 2>{log}
        '''
        