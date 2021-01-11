## Metagenomics Sequencing Analysis Workflow



To obtain the metagenome composition and function annotation profile via the popular software such as metaphlan and humman version 3 is very helpful when we do some research focusing on the known microorganisms in the environment like the gut. 

<img src="./mindmap_workflow.jpg" width="1000" height="600">



### Installation

This pipeline is based on the metaphlan and humann version 3, so please install the two software before running it. In addition, removing the host DNA sequences is used the kneaddata and the user should download the dependent files. All the softwares could be installed by conda

```bash   
# install software 
conda create --name metagenome metaphlan humann kneaddata fastqc multiqc fastp -y
```



### Download Referecne database 

After installing the software, we should configure all the reference

```bash
# configure the database
# metaphlan 
# humann 
# kneaddata
```



### How  to run

#### input file

the fastqpath file

```bash
find /RawData/ -name "*fq.gz" |  sort | perl -e 'print "SampleID\tLaneID\tPath\n"; while(<>){chomp; $fq=(split("\/", $_))[-1]; $sampleid=$fq; $laneid=$fq; $sampleid=~s/\_R[1|2]\.fq.gz//g; $laneid=~s/\.fq.gz//g;print "$sampleid\t$laneid\t$_\n";}' > samples.fqpath.tsv
```



| SampleID | LaneID | Path                 |
| -------- | ------ | -------------------- |
| ND2      | ND2_R1 | RawData/ND2_R1.fq.gz |
| ND2      | ND2_R2 | RawData/ND2_R2.fq.gz |

#### command line 

```bash
perl main.pl -f test.tsv -a TruSeq2-PE.fa -o Run.all.sh
```



### Final directory structure

```bash
./MetaGenomics/
├── bin
│   ├── convert2matrix.pl
│   ├── humann.pl
│   ├── kneaddata.pl
│   ├── merge.pl
│   ├── metabolic_pathway.pl
│   ├── metaphlan.pl
│   └── qc.pl
├── main.pl
├── result
│   ├── 00.quality
│   │   ├── fastqc
│   │   └── multiqc
│   ├── 01.kneaddata
│   ├── 02.merge
│   ├── 03.humann
│   │   ├── genefamilies
│   │   ├── log
│   │   ├── metaphlan
│   │   ├── pathabundance
│   │   └── pathcoverage
│   ├── 04.metaphlan
│   ├── 05.profile
│   ├── 06.metabolic_pathway
│   │   ├── EggNOG_COGs
│   │   ├── genefamilies
│   │   ├── Gene_Ontology
│   │   ├── KEGG_Orthogroups
│   │   ├── Level4_enzyme
│   │   └── MetaCyc_Reactions
│   ├── Run.s1.qc.sh
│   ├── Run.s2.kneaddata.sh
│   ├── Run.s3.merge.sh
│   ├── Run.s4.humann.sh
│   ├── Run.s5.metaphlan.sh
│   ├── Run.s6.profile.sh
│   ├── Run.s7.metabolic.sh
│   └── script
│       ├── 00.quality
│       ├── 01.kneaddata
│       ├── 02.merge
│       ├── 03.humann
│       ├── 04.metaphlan
│       ├── 05.profile
│       └── 06.metabolic_pathway
├── Run.all.sh
├── test.tsv
├── TruSeq2-PE.fa 
├── util
│   ├── calculate_unifrac.R
│   └── mpa_v30_CHOCOPhlAn_201901_species_tree.nwk
└── work.sh

31 directories, 21 files
```



### Contributors

-   [Hua Zou](https://github.com/zouhua)

