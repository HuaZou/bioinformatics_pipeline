#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 02 09:34:07 2016

This script will combine the FAST-derep map with the FAST-OTU map.
The new FAST map is a hybrid of these two that contains all sample information.

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
                                        ------------------------'''),prog='fast.py -combine_fast_map')
    parser.add_argument('-derep_map', help='Name of the FAST-derep map.')
    parser.add_argument('-otu_map', help='Name of the FAST-OTU map.')
    parser.add_argument('-o', '--output', help='Name of the output hybrid FAST map.')
    
    #parser.add_argument('-separator', default = ';', help='Set the separator for parsing the sequence label.')

    args = parser.parse_args(name_space)
    
    input_derep_map_file = args.derep_map
    input_otu_map_file = args.otu_map
    
    output_hybrid_fast_file = args.output
    
    from lib import ParseOtuMap
    
    input_derep_map = ParseOtuMap.read_fast_output(input_derep_map_file)
    input_otu_map = ParseOtuMap.read_fast_output(input_otu_map_file)
    
    output_hybrid_fast = ParseOtuMap.merge_fast_output(input_derep_map, input_otu_map)
    
    ParseOtuMap.write_fast_output(output_hybrid_fast, output_hybrid_fast_file)
    
    print('FAST style map file wrote to %s.' %output_hybrid_fast_file)
    
if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
