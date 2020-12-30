## Metagenomics Sequencing Analysis Workflow



### Summary of metagenomics workflow

<img src="./mindmap_workflow.jpg" width="1000" height="600">





The workflow consists of two parts. 

* the software involved in Part1:

  * fastqc/multiqc
  * kneaddata
  * fastp 
  * MetaPhlAn version 3
  * HumaAn version 3

  

### Structure

```bash
# tree -L 3 /metagenomics/
/metagenomics/
├── README.md
├── bin
│   ├── humann.pl
│   ├── kneaddata.pl
│   ├── merge.pl
│   ├── metaphlan.pl
│   └── qc.pl
├── example
│   ├── Run.all.sh
│   ├── TruSeq2-PE.fa
│   ├── result
│   │   ├── 00.quality/
│   │   ├── 01.kneaddata/
│   │   ├── 02.merge/
│   │   ├── 03.humann/
│   │   ├── 04.metaphlan/
│   │   ├── Run.s1.qc.sh
│   │   ├── Run.s2.kneaddata.sh
│   │   ├── Run.s3.merge.sh
│   │   ├── Run.s4.humann.sh
│   │   ├── Run.s5.metaphlan.sh
│   │   └── script/
│   ├── test.tsv
│   └── work.sh
├── main.pl
├── mindmap_workflow.emmx
└── mindmap_workflow.jpg
```

