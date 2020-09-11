#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 08 16:38:10 2016



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
    from lib import ParseOtuTable
    import sys
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                    ------------------------
                                    By Zewei Song
                                    University of Minnesota
                                    Dept. Plant Pathology
                                    songzewei@outlook.com
                                    ------------------------'''), prog='fast.py -parse_uc_cluster')
    parser.add_argument('-otu', help='Input OTU table')
    parser.add_argument('-o', '--output', help='Output OTU table')
    parser.add_argument('-sample_list', help='A list of new sample order')
    parser.add_argument('-meta_column', default='taxonomy', help='Name of the meta column')
    
    args = parser.parse_args(name_space)
    input_file = args.otu
    output_file = args.output
    input_names = args.sample_list
    
    # Read in the sample list
    sample_list = []
    with open(input_names ,'rU') as f:
        for line in f:
            line = line.strip('\n')
            sample_list.append(line)
    print('Found {0} samples to be re-ordered.'.format(len(sample_list)))
    
    # Read in the OTU table
    otu_old = ParseOtuTable.parser_otu_table(input_file, meta_col = args.meta_column)
    
    # Get the sample list in the current OTU table
    current_sample_list = otu_old.sample_id
    sample_failed = []
    for sample in sample_list:
        try:
            current_sample_list.index(sample)
        except ValueError:
            print('{0} is not found in current OTU table.'.format(sample))
            sample_failed.append(sample)
    
    if len(sample_failed) > 0:
        print('{0} samples are not in current OTU table, please check your sample list.'.format(len(sample_failed)))        
        sys.exit()
    
    # Output the new OTU table with the new sample order
    ParseOtuTable.write_sample_dict_newlist(otu_old.sample_dict(),otu_old.meta_dict(),otu_old.species_id, output_file, sample_list)
    print('A new OTU table was written to {0}.'.format(output_file))

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
