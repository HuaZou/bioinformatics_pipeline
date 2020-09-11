# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 12:58:05 2015

Merged all sequence files (FASTA or FASTQ) in a folder into a single file.

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
    import argparse
    from lib import File_IO
    import sys
    import os
    import time
    import textwrap
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                    ------------------------
                                    By Zewei Song
                                    University of Minnesota
                                    Dept. Plant Pathology
                                    songzewei@outlook.com
                                    ------------------------'''), prog='fast.py -merge_seqs')
    parser.add_argument('-i', '--input', help='Name of the input folder containing files to be merged')
    parser.add_argument('-o', '--output', help='Name of the merged file')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-fasta', action='store_true', help='Set the file type to FASTA')
    group.add_argument('-fastq', action='store_true', help='Set the file type to FASTQ, this is the default option')
    args = parser.parse_args(Namespace)
    
    input_folder = args.input
    if not input_folder:
        print('please specified an input folder')
        sys.exit()
    output_file = args.output
    if not output_file:
        print('Please specified an output file.')
        sys.exit()
    if os.path.isfile(output_file):
        file_size = round(os.path.getsize(output_file)/1024**2, 0)
        exist = raw_input('%s (%d MB)already exists , do you want to overwrite it? [y/n]' % (output_file, file_size))
        if exist == 'y' or exist == 'Y':
            os.remove(output_file)
        else:
            print('Program stopped.')
            sys.exit()
    file_type = 'fastq'
    if args.fasta:
        file_type = 'fasta'
    
    start = time.time()
    f_list = File_IO.file_list(input_folder)
    f_list.sort()
    print('Found %i files in the folder %s' % (len(f_list), input_folder))
    count = 0
    n = 1
    count_total = 0
    for seq_file in f_list:
        current_file = input_folder+'/'+seq_file
        count = File_IO.write_seqs(File_IO.read_seqs(current_file, file_type), output_file, checker=False, overwrite=False)
        print('%d. Merged %d sequences from %s into the new file.' % (n, count, seq_file))
        n += 1
        count_total += count
    end = time.time()
    used_time = round(end-start, 2)
    print('Spent %s sec to merge %d records in %d files into %s' % (str(used_time), count_total, len(f_list), output_file))

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])