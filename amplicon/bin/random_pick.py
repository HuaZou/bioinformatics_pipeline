#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 10:18:19 2016



Please feel free to contact me for any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
www.songzewei.org
"""
from __future__ import print_function
def main(name_space):
    from lib import random_subsample as rs
    from lib import File_IO
    import argparse
    import textwrap
    import sys
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                        ------------------------
                                        By Zewei Song
                                        University of Minnesota
                                        Dept. Plant Pathology
                                        songzewei@outlook.com
                                        ------------------------'''), prog='fast.py -random_subsample')
    parser.add_argument('-r1', help='Name of the Read1 file.')
    parser.add_argument('-r2', help='Name of the Read2 file if applicable.')
    parser.add_argument('-size', default = 10000, help='Sampling size for each file, default=10,000.')
    
    args = parser.parse_args(name_space)
    
    read1 = args.r1
    if args.r2:
        read2 = args.r2
    sample_size = int(args.size)
    
    read1_content = File_IO.read_seqs(read1)
    total_size = len(read1_content)
    
    file_type = "fasta"
    if read1_content[0][2] == "+":
        file_type = "fastq"
    
    if sample_size > total_size:
        print('The specified sampling size is larger than the total number of sequences.')
        sys.exit()
    else:
        seq_index = rs.generate_random_index(total_size, sample_size)
    
    # Get sequences in read1 file
    read1_picked = []
    for index in seq_index:
        read1_picked.append(read1_content[index])
    
    # Pick read1 file is the filename is specified
    if args.r2:
        read2_content = File_IO.read_seqs(read2)
        read2_picked = []
        for index in seq_index:        
            read2_picked.append(read2_content[index])

    # write to new files
    read1_output = "R1."+file_type
    read1_count = File_IO.write_seqs(read1_picked, read1_output, checker=False, overwrite=True)
    print('{0} sequences have been randomly picked from {1}, and saved in {2}.'.format(read1_count, read1, read1_output))
    if args.r2:
        read2_output = "R2."+file_type
        read2_count = File_IO.write_seqs(read2_picked, read2_output, checker=False, overwrite=True)
        print('{0} sequences have been randomly picked from {1}, and saved in {2}.'.format(read2_count, read2, read2_output))

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
