# -*- coding: utf-8 -*-
"""
Created on Mon Apr 06 14:30:04 2015

Filter a Qiime style OTU map based on minimum size or sequence name.

User can provide 1) a minimum size to be kept, 2) a file contains names in separated lines, or
3) A fasta file contains sequences to be kept.

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
    import textwrap
    from lib import ParseOtuMap
    from lib import File_IO
    import sys
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                    ------------------------
                                    By Zewei Song
                                    University of Minnesota
                                    Dept. Plant Pathology
                                    songzewei@outlook.com
                                    ------------------------'''), prog='fast.py -filter_otu_map')
    parser.add_argument('-i', '--input', help='Input OTU map')
    parser.add_argument('-o', '--output', help='Output OTU map')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-min_size', default=2, help='The minimum size of an OTU to be kept')
    group.add_argument('-name_list', help='A file contains a list of sequence name to be picked')
    group.add_argument('-fasta', help='A FASTA file contains sequence to be picked')
    args = parser.parse_args(Namespace)
    
    map_file = args.input
    output_file = args.output
    if args.min_size:
        min_size = int(args.min_size)
    if args.name_list:
        with open(args.name_list, 'rU') as f:
            pick_list = []
            for line in f:
                pick_list.append(line.strip('\n'))
    if args.fasta:
        seqs = File_IO.read_seqs(args.fasta)
        pick_list = []
        for record in seqs:
            pick_list.append(record[0])
    
    print('Reading in %s ...' % map_file)
    MapDict = ParseOtuMap.read_otu_map(map_file)
    
    # Filter OTU map based on parameters
    if args.min_size:
        print('Filtering OTUs with less than %d sequences ...' % min_size)
        MapDictFiltered = ParseOtuMap.filter_by_size(MapDict, min_size=min_size)
    
    if args.name_list or args.fasta:
        print('Pick OTUs based on the names in %s ...' % args.name_list)
        MapDictFiltered = {}
        for name in pick_list:
            try:
                MapDictFiltered[name] = MapDict[name]
            except KeyError:
                print('Cannot find %s in the OTU map. Program exits.')
                sys.exit()
    
    #  Report comparison of original and filtered maps
    old_map = ParseOtuMap.otu_map_parser(MapDict)
    new_map = ParseOtuMap.otu_map_parser(MapDictFiltered)
    
    print('\n')
    print('Original OTU map:')
    print('\t OTU=%i (Total Sequences=%i, Max=%i, Min=%i, Ave=%i)' % (
        old_map.derep_count, old_map.seqs_count, old_map.max_derep, old_map.min_derep, old_map.ave_derep))
    print('Filtered OTU map:')
    print('\t OTU=%i (Total Sequences=%i, Max=%i, Min=%i, Ave=%i)' % (
        new_map.derep_count, new_map.seqs_count, new_map.max_derep, new_map.min_derep, new_map.ave_derep))
    
    print('Writing new map ...')
    ParseOtuMap.write_otu_map(MapDictFiltered, output_file=output_file)
    print('New map saved in %s.' % output_file)

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])