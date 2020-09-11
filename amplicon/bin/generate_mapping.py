# -*- coding: utf-8 -*-
"""
Created on Wed Apr 01 15:25:06 2015

Try to guess the sample name and read type from the raw Illumina FASTQ files.
The format is according to that used by University of Minnesota Genomic Center (UMGC).
The output is a mapping file which can be used as the input of add_labels.py

This is NOT a Qiime mapping file for add_qiime_labels.py, but should be easy to convert to that file.

Please feel free to contact me for any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""
from __future__ import print_function
from __future__ import division

def main(Namespace):

    def parse_file_name(filename):
        if filename.find('_') != -1:
            sample_name = filename[0:filename.find("_")]
        else:
            sample_name = filename
        if filename.find("R1") != -1:
            read_type = "R1"
        elif filename.find("R2") != -1:
            read_type = "R2"
        else:
            read_type = "unknown"
        return (filename, sample_name, read_type)
    
    
    def write_mapping(folder, mapping_file='mapping.txt'):
        f = File_IO.file_list(folder)
        #mapping = [("Filename", "Sample_name", "Read_type")]
        mapping = [('#SampleID','BarcodeSequence','LinkerPrimerSequence','InputFileName','ReadType','Description')]
        unknown_read = 0
        count_file = 0
        for item in f:
            file_info = parse_file_name(item)
            new_record = (file_info[1],'','',file_info[0],file_info[2],file_info[1])
            mapping.append(new_record)
            count_file += 1
            if mapping[-1][4] == 'unknown':
                unknown_read += 1
        with open(mapping_file, 'w') as m:
            for line in mapping:
                m.write("%s\n" % '\t'.join(line))
        return count_file, unknown_read

    import argparse
    import textwrap
    from lib import File_IO

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                    ------------------------
                                    By Zewei Song
                                    University of Minnesota
                                    Dept. Plant Pathology
                                    songzewei@outlook.com
                                    ------------------------'''), prog='fast.py -generate_mapping')
    parser.add_argument('-i', '--input', help='Name of the folder with the raw data.')
    parser.add_argument('-o', '--output', default='mapping.txt', help='Name of the new mapping file.')
    args = parser.parse_args(Namespace)

    folder_name = args.input
    count_file, unknown_read = write_mapping(folder_name, mapping_file=args.output)

    print('Generated a mapping file with %d files in %s.' % (count_file, args.output))
    if unknown_read > 0:
        print('%s files have unknown read type, please check the mapping file.' % unknown_read)
        
if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
