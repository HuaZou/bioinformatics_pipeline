
# how to get transcript id and gene id 
# devtools::install_github("BUStools/BUSpaRse")
library(BUSpaRse)
tr2g <- transcript2gene(
            c("Homo sapiens", "Mus musculus"), 
            type = "vertebrate",
            ensembl_version = 100, 
            kallisto_out_path = "./")

# get geneid and length on Homo_sapiens.GRCh38.101.gtf
#cat Homo_sapiens.GRCh38.101.gtf | awk -F'\t' '{if($3=="gene") {split($9,a,";"); print a[1]"\t"$5-$4};}' | sed 's/[gene_id |"|]//g' | sort -u > Homo_sapiens.GRCh38.101.genelength.tsv

library(biomaRt)
library(curl)

genelist <- read.table("Homo_sapiens.GRCh38.101.genelength.tsv", header = T)
human_mart <- useMart(host="www.ensembl.org", 
                     biomart="ENSEMBL_MART_ENSEMBL", 
                     dataset = "hsapiens_gene_ensembl")
human_gene_all <- getBM(attributes=c("ensembl_gene_id", 
                                     "entrezgene_id",
                                     "external_gene_name", 
                                     "ensembl_transcript_id", 
                                     "ensembl_transcript_id_version",
                                     "transcript_biotype", 
                                     "description"),
                            filters="ensembl_gene_id",
                            values = genelist$Geneid,
                            mart=human_mart)
