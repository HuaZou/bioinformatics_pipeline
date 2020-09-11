#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

Please feel free to contact me for any question.
Fri Jul 21 11:06:57 2017
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
www.songzewei.org
"""
#%%
from __future__ import print_function
from __future__ import division

def main(name_space):
    import argparse
    import textwrap
    from lib import File_IO
    import os
    import sys

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                        ------------------------
                                        By Zewei Song
                                        University of Minnesota
                                        Dept. Plant Pathology
                                        songzewei@outlook.com
                                        ------------------------'''), prog='fast.py -add_labels')
    parser.add_argument("-i", "--input", help="Name of the input file, merging from multiple sample sequences.")
    parser.add_argument("-o", "--output", help="Name of the output folder.")
    parser.add_argument("-r", "--read", choices = ['r1', 'read1', 'r2', 'read2'], help="Read direction, read1 or read2.")
    args = parser.parse_args(name_space)
    #args = argparse.Namespace(input = 'read1.cut2.fastq', output = 'unmerged', read = 'read1') # This line is for testing purpose

    input_file = args.input
    output_folder = args.output
    read_type = args.read

    if read_type == 'r1' or read_type == 'read1':
        read_type = 'R1'
    elif read_type == 'r2' or read_type == 'read2':
        read_type = 'R2'
    else:
        print('Please specify the correct read type using the -r option.')
        sys.exit()

    os.makedirs(output_folder, exist_ok = True)
    input_seqs = File_IO.read_seqs(input_file)
    output_records = {}

    for record in input_seqs:
        sample_name = record[0][0:record[0].index('_')]
        try:
            output_records[sample_name].append(record)
        except KeyError:
            output_records[sample_name] = [record]

    for key, value in output_records.items():
        output_file = output_folder + '/' + key + '_' + read_type + '.fastq'
        File_IO.write_seqs(value, output_file, checker = False, overwrite = True)

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
    #a = main() # This line is for testing purpose