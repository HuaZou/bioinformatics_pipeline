# step1 remove adapter
#find /disk/user/zouhua/pipeline/miRNA_seq/00.RawData/ -name "*.fq.gz" | perl -e 'print"SampleID\tLaneID\tPath\n";while(<>){chomp; $name=(split("\/", $_))[-1]; $name1=$name; $name2=$name; $name1=~s/_miRNA.*//g; $name2=~s/.fq.gz//; print "$name1\t$name2\t$_\n";}' > miRNA.samples.path.tsv 
#ln -s ../dataset/01.miRNA_seq/4samples.path.tsv  miRNA.samples.path.tsv 
perl batch_fqc.pl -f miRNA.samples.path.tsv -o 01.fastqc.sh
perl batch_trim.pl -f miRNA.samples.path.tsv -a AACTGTAGGCACCATCAAT -a2 AGATCGGAAGAG -o 02.trim.sh
#rm miRNA.samples.path.tsv

# step2 2nd fastqc
find /disk/user/zouhua/pipeline/miRNA_seq/result/02.trim/trimmomatic/ -name "*.fq.gz" | perl -e 'print"SampleID\tLaneID\tPath\n";while(<>){chomp; $name=(split("\/", $_))[-1]; $name1=$name; $name2=$name; $name1=~s/.clean.*//g; $name2=~s/.fq.gz//; print "$name1\t$name2\t$_\n";}'  > miRNA_clean.samples.path.tsv
perl batch_fqc_2nd.pl -f miRNA_clean.samples.path.tsv -o 02.trim_fastqc.sh
perl batch_collapse.pl -f miRNA_clean.samples.path.tsv -o 03.collapse.sh
rm miRNA_clean.samples.path.tsv 

# step3 remove ncRNA
find /disk/user/zouhua/pipeline/miRNA_seq/result/03.collapse/collapse/ -name "*.fa" | perl -e 'print"SampleID\tLaneID\tPath\n";while(<>){chomp; $name=(split("\/", $_))[-1]; $name1=$name; $name1 =~ s/\.collapse\.fa//; print "$name1\t$name1\t$_\n";}' > miRNA_collapse.samples.path.tsv
perl batch_blast.pl -f miRNA_collapse.samples.path.tsv -d /disk/share/database/Mus_Rfam_ncRNA/blast_index -o 04.blast.sh
rm miRNA_collapse.samples.path.tsv

# step4 remove miRNA in exon
find /disk/user/zouhua/pipeline/miRNA_seq/result/04.blast/ -name "*.fa" | perl -e 'print"SampleID\tLaneID\tPath\n";while(<>){chomp; $name=(split("\/", $_))[-1]; $name1=$name; $name1 =~ s/\.remained\.fa//; print "$name1\t$name1\t$_\n";}' > miRNA_remained.samples.path.tsv
perl batch_mapped_intron.pl -f miRNA_remained.samples.path.tsv -d /disk/share/database/Musculus_UCSC/bowtie_index -o 05.exon_intron.sh
rm miRNA_remained.samples.path.tsv

# step5 mirDeep2
find /disk/user/zouhua/pipeline/miRNA_seq/result/05.exon_intron/ -name "*.fa" | perl -e 'print"SampleID\tLaneID\tPath\n";while(<>){chomp; $name=(split("\/", $_))[-1]; $name1=$name; $name1 =~ s/\_exonintron\_filtered\.fa//; print "$name1\t$name1\t$_\n";}' > miRNA_exon.samples.path.tsv
perl batch_mirdeep2.pl -f miRNA_exon.samples.path.tsv -d /disk/share/database/Mus_musculus.GRCm38_release100/bowtie_remove_whitespace_index -fa /disk/share/database/Mus_musculus.GRCm38_release100/Mus_musculus.GRCm38.dna_sm.primary_assembly_remove_whitespace.fa -m /disk/share/database/mirbase -o 06.mirdeep2.sh
rm miRNA_exon.samples.path.tsv
