#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 03 11:07:15 2016



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
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                        ------------------------
                                        By Zewei Song
                                        University of Minnesota
                                        Dept. Plant Pathology
                                        songzewei@outlook.com
                                        ------------------------'''), prog='fast.py -random_subsample')
    parser.add_argument('-i', '--input', help='Name of the input folder with raw data')
    parser.add_argument('-o', '--output', default = 'random_dataset', help='Name of the output folder with raw data')
    parser.add_argument('-file_number', default = 10, help='Number of file to pick.')
    parser.add_argument('-size', default = 10000, help='Sampling size for each file.')
    
    args = parser.parse_args(name_space)
    
    input_folder = args.input
    output_folder = args.output
    
    file_number = int(args.file_number)
    sample_size = int(args.size)
    
    # Create new folder
    File_IO.mk_dir(output_folder)
    
    # Randomly pick files to be sampled
    input_file_list = File_IO.file_list(input_folder)
    print('Found {0} files in the folder {1}'.format(len(input_file_list), input_folder), end = '\n')
    file_index = rs.generate_random_index(len(input_file_list), file_number)
    file_list = []
    for index in file_index:
        file_list.append(input_file_list[index])
    
    # Randomly pick sequences from each file
    for raw_file in file_list:
        print('\tRandoming sampling {0} for {1} sequences ...'.format(raw_file, sample_size, end='\r'))
        current_content = File_IO.read_seqs(input_folder + '/' + raw_file)
        seq_index = rs.generate_random_index(len(current_content), sample_size)
        sampled_content = []

        for index in seq_index:
            sampled_content.append(current_content[index])
        
        count = File_IO.write_seqs(sampled_content, output_folder + '/' + raw_file)
    
    print('A randomly sampled dataset ({0} files, {1} sequences per file) was generated under the folder {2}'.format(file_number, sample_size, output_folder, end = '\n'))
    
if __name__ == '__main__':
    import sys
    main(sys.argv[1:])