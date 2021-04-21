#Rscript DNA_Methylation_preprocessing.R -f IDAT -t RawData_IDAT.tsv -o result
Rscript DNA_Methylation_preprocessing_v2.R -f IDAT -t RawData_IDAT.tsv -o result
Rscript DNA_Methylation_preprocessing_v3.R -f IDAT -t RawData_IDAT_v2.csv -c 850k -o result_v2
