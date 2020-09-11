#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 02 09:01:07 2016

This script will generate a new FAST style map file. This map file is similar to a Qiime style OTU map,
but also contains the information about sequences. It could be a map for either dereplicated seuqences,
or OTUs. If it is for OTUs, it can also contain information about the detail dereplicated unit within each
OTUs (abudnance and sequence).

The map is saved in JSON format, and can be read using the json module in python.

Please feel free to contact me for any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
www.songzewei.org
"""
from __future__ import print_function
from __future__ import division

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
                                        ------------------------'''), prog='fast.py -generate_fast_map')
    parser.add_argument('-map', help='Name of the Qiime style OTU/Derep map file')
    parser.add_argument('-seq', help='Name of the sequence file corresponding to the Qiime map.')
    parser.add_argument('-o', '--output', help='Name of the output FAST map.')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-derep', action='store_true', help='Indicate the source is a dereplication map.')
    group.add_argument('-otu', action='store_true', help='Indicate the source is an OTU map.')
    parser.add_argument('-separator', default = ';', help='Set the separator for parsing the sequence label.')
    
    args = parser.parse_args(name_space)
    
    input_map_file = args.map
    input_seq_file = args.seq
    output_map_file = args.output
    separator = args.separator
    
    if args.derep:
        real_sample = True
    elif args.otu:
        real_sample = False
    
    from lib import ParseOtuMap
    from lib import File_IO
    
    input_map = ParseOtuMap.read_otu_map(input_map_file)
    input_seq = File_IO.read_seqs(input_seq_file)
    
    output_map = ParseOtuMap.generate_fast_output(input_map, input_seq, real_sample = real_sample, separator=separator)
    
    ParseOtuMap.write_fast_output(output_map, output_map_file)
    
    print('FAST style map file wrote to %s.' %output_map_file)

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])