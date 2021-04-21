# Methylation_analysis_scripts
R scripts for Illumina array analyses

## Pre-processing: *Methylation_pre-processing.R*

This script performs Illumina EPIC 850K array and Illumina 450K pre-processing and QC from idat files. 

## Prerequisites
This R script requires the following packages:
- BiocManager
- minfi
- ENmix
- wateRmelon
- MASS
- broom
- RColorBrewer
- IlluminaHumanMethylationEPICmanifest
- IlluminaHumanMethylationEPICanno.ilm10b4.hg19
- IlluminaHumanMethylation450kmanifest
- IlluminaHumanMethylation450kanno.ilm12.hg19

```
install.packages("BiocManager")
BiocManager::install("minfi")
BiocManager::install("ENmix")
BiocManager::install("wateRmelon")
install.packages("MASS")
install.packages("broom")
install.packages("RColorBrewer")
BiocManager::install("IlluminaHumanMethylationEPICmanifest")
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
BiocManager::install("IlluminaHumanMethylation450kmanifest")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
```

### Usage
```bash
Rscript Methylation_pre-processing.R -f input_folder -t targetfile.txt [options]
```

| **PARAMETER** | **DEFAULT** | **DESCRIPTION** |
|-----------|--------------:|-------------| 
*-f* | . | folder with idat files |
*-o* | out | output directory name |
*-t* |  . | target file with sample information : sample_id, plate_id (can be null), sample_well (can be null), sentrix_id (like "3999492078"), sentrix_position (like "R01C01"), sex (can be null) |
*-p* |  \t | target file field separator |
*-c* | NULL | file with list of cross-reactive probes |
*-s* | TRUE | filter SNP-probes (maf >= 0.05) |
*-m* | FALSE | filter multimodal probes |
*-h*    |  | Show help message and exit|

### Details
The script involves 5 steps
- **Reading idat files**
- **Raw data QC** using minfi and ENmix packages
- **Data normalization** using functional normalisation 
- **Data filtering** of poor quality and cross reactive probes, XY chromosome probes and SNP-associated probes
- **Making analysis-ready beta and m tables** 

![DNAmethylationProcess](https://github.com/IARCbioinfo/Methylation_analysis_scripts/blob/master/DNAmethylationProcess.png)

### Output
- RGSet.RData file
- MSet.RData file
- gset.RData file
- sva90.RData file
- PdetectionTables.RData file
- Fun/Fun1/Fun2/Fun3.RData files
- Normalised and filtered beta and m value tables

### QC output plots
| **PLOT NAME** | **DESCRIPTION** |
|-----------|-------------| 
*STAINING Green/Red* | Staining controls are used to examine the efficiency of the staining step in both the red and green channels and are independent of the hybridisation and extension step |
*EXTENSION Green/Red* | Extension controls test the extension efficiency of A, T, C and G nucleotides from a hairpin probe, and are therefore sample-independent. A and T should be monitored in the red channel, C and G should be monitored in the green channel. |
*TARGET REMOVAL Green/Red* | Target removal controls test the efficiency of the stripping step after the extension reaction. The probe sequences are designed such that extension from the probe does not occur. All target removal controls should result in low signal compared to the hybridisation controls, indicating that the targets were removed efficiently after extension. Target removal controls are present in the hybridisation buffer RA1. The performance of the target removal controls should be monitored only in the green channel |
*HYBRIDIZATION Green/Red* | Hybridisation controls test the overall performance of the Infinium assay using synthetic targets instead of amplified DNA. These synthetic targets complement the sequence on the array perfectly, allowing the probe to extend on the synthetic target as a template. Synthetic targets are present in the hybridisation buffer at three levels, monitoring the response from high concentration, medium concentration and low concentration targets. All bead type IDs should result in signal with various intensities, corresponding to the concentrations of the initial synthetic targets. The performance of the hybridisation controls should be monitored only in the green channel.|
*BISULFITE CONVERSION I Green/Red* | These controls use Infinium I probe design and monitor the efficiency of bisulphite conversion. There are no C bases in the probe except for the site in question. If the conversion reaction was successful, the ‘converted’ C probes will match the converted sequence and get extended, if there is unconverted DNA the ‘unconverted’ U probes will get extended. C1, and C2 should be monitored in the green channel, C3, C4 and C5 should be monitored in the red channel. |
*BISULFITE CONVERSION II Green/Red* | These controls use Infinium II probe design and single base extension to monitor efficiency of bisulfite conversion. If the bisulfite conversion reaction was successful, the "A" base will get incorporated and the probe will have intensity in the Red channel. If the sample has unconverted DNA, the "G" base will get incorporated across the unconverted cytosine, and the probe will have elevated signal in the Green channel.|
*SPECIFICITY I Green/Red* | These controls are designed to monitor allele-specific extension for Infinium I probes. In this probe design, the A/T match corresponds to the unmethylated status of the interrogated C, and G/C match corresponds to the methylated status of C. G/T mismatch controls check for non-specific detection of methylation signal over unmethylated background. PM controls correspond to A/T perfect match and should give high signal. MM controls correspond to G/T mismatch and should give low signal. High in green channel: GT mismatch 1 PM, GT mismatch 2 PM and GT mismatch 3 PM. High in red channel: GT mismatch 4 PM, GT mismatch 5 PM, GT mismatch 6 PM |
*SPECIFICITY II Green/Red* | These controls are designed to monitor allele-specific extension for Infinium II probes and check for potential non-specific detection of methylation signal over unmethylated background. Specificity II probes should incorporate the A base across the non-polymorphic T and have intensity in the red channel. In case of nonspecific incorporation of the G base, the probe will have elevated signal in the green channel. |
*NEGATIVE Green/Red* | Negative controls target bisulphite-converted sequences that do not contain CpG dinucleotides. Assay probes are randomly permutated and should not hybridise to the DNA template. The mean signal of these probes defines the system background.|
*NON-POLYMORPHIC Green/Red* | Non-polymorphic controls test the overall performance of the assay by querying a particular base in a non-polymorphic region of the genome. They let you compare assay performance across difference samples. One non-polymorphic control has been designed for each of the four nucleotides. |

