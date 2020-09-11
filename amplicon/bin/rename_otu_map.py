#!usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create on 4/27/2015.

Rename the OTU in a Qiime or FAST style OTU map using a provided label.

Please feel free to contact me for any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""
from __future__ import print_function
from __future__ import division

def main(name_space):
    import argparse
    import textwrap
    from lib import ParseOtuMap
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                    ------------------------
                                    By Zewei Song
                                    University of Minnesota
                                    Dept. Plant Pathology
                                    songzewei@outlook.com
                                    ------------------------'''), prog='fast.py -rename_otu_map')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-qiime_map', help='Input Qiime style OTU map')
    group.add_argument('-fast_map', help='Input FAST style OTU map')
    parser.add_argument('-o', '--output', help='Output OTU map')
    parser.add_argument('-label', default='OTU_', help='New OTU name label')
    args = parser.parse_args(name_space)
    
    if args.qiime_map != None:
        input_otu = args.qiime_map
        method = 'qiime'
    elif args.fast_map != None:
        input_otu = args.fast_map
        method = 'fast'
    
    if method == 'qiime':
        input_otu = args.qiime_map
        output_otu = args.output
        otu_label = args.label
        
        print('Reading in %s ...' % input_otu)
        old_map = ParseOtuMap.read_otu_map(input_otu)
        otu_temp = old_map.items()
        print('Sorting OTUs by size ...')
        otu_temp.sort(key=lambda x: len(x[1]), reverse=True)
        
        with open(output_otu, 'w') as f:
            count = 1
            for otu in otu_temp:
                new_name = otu_label + str(count)
                count += 1
                line = new_name + '\t' + '\t'.join(otu[1])
                f.write('%s\n' % line)
        print('New OTU map saved in %s.' % output_otu)
    
    elif method == 'fast':
        output_otu = args.output
        fast_otu_map = ParseOtuMap.read_fast_output(input_otu)
        map_parser = ParseOtuMap.fast_output_parser(fast_otu_map)
        seq_list = map_parser.get_seqs()
        
        new_map = {}
        count = 1
        for unit in seq_list:
            new_otu_name = args.label + str(count)
            new_map[new_otu_name] = fast_otu_map[unit[0]]
            count += 1
        ParseOtuMap.write_fast_output(new_map, output_otu)

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])