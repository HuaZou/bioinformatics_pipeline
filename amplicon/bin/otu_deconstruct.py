#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 02 09:56:28 2016

This script will output the detail structure of each individual OTU using the information
of a hybrid FAST map.

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
    import argparse
    import textwrap
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                        ------------------------
                                        By Zewei Song
                                        University of Minnesota
                                        Dept. Plant Pathology
                                        songzewei@outlook.com
                                        ------------------------'''), prog = 'fast.py -otu_deconstruct')
    parser.add_argument('-map', help='Name of the FAST-derep map.')
    parser.add_argument('-o', '--output', default = 'otu_deconstruct', help='Name of the output folder')

    args = parser.parse_args(name_space)
    
    input_map_file = args.map
    output_folder = args.output
    
    from lib import ParseOtuMap
    from lib import File_IO
    
    File_IO.mk_dir(output_folder)
    
    input_map = ParseOtuMap.read_fast_output(input_map_file)
    
    input_map = ParseOtuMap.fast_output_parser(input_map)
    
    input_map_size = input_map.unit_count
    print('{0} contains {1} OTUs.'.format(input_map_file, input_map_size))
    
    otu_list = input_map.get_seqs() # get a list of otu with their sequences
    
    for unit in otu_list:
        output_file = output_folder + '/' + unit[0] + '.txt'
        current_otu = input_map.detail_sample_unit(unit[0])
        
        print('\tWriting: {0} ...\r'.format(output_file, end='\r'))
        with open(output_file, 'wb') as f:
            for line in current_otu:
                line = '\t'.join([str(i) for i in line])
                f.write('%s\n' % line)
    
    print('All files wrote to the folder: {0}.'.format(output_folder))

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
