#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create  on April 25th, 2015.

Convert a UPARSE denovo cluster output file into a Qiime style OTU table.

The second column of the input file should contains 'otu, match, or chimera".

Please feel free to contact me with any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""
from __future__ import print_function
def main(name_space):
    
    import argparse
    import textwrap
    from lib import ParseOtuMap
    import sys
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                    ------------------------
                                    By Zewei Song
                                    University of Minnesota
                                    Dept. Plant Pathology
                                    songzewei@outlook.com
                                    ------------------------'''), prog='fast.py -parse_uc_cluster')
    parser.add_argument('-i', '--input', help='Input uparse denovo cluster output file')
    parser.add_argument('-o', '--output', help='Output Qiime style OTU map')
    parser.add_argument('-separator', default=';', help='Separator between sample name and size annotation')
#    parser.add_argument('-fast', default='', help='Name of FAST style output.')
#    parser.add_argument('-centroid', default = '', help='Name of the FASTA file contains sequences of centroids.')
    
    args = parser.parse_args(name_space)
    input_file = args.input
    output_file = args.output
    separator = args.separator
    
    # Read in the uparse file
    uc_content = []
    count_uc = 0
    with open(input_file, 'rU')  as f:
        for line in f:
            temp = line.strip('\n').split('\t')
            #temp[0] = temp[0][:temp[0].find(separator)]  # Remove size annotation and other extra labels
            uc_content.append(temp)
            count_uc += 1
            #line = '%i line in the uparse file ... \r' % count_uc
            #sys.stderr.write(line)
    
    # Convert uparse output into Qiime style OTU map
    otu_map = {}
    
    for line in uc_content:
        centroid = ""
        current_seq = ""
        if line[0] == 'S':
            centroid = line[8][:line[8].find(separator)]
            otu_map[centroid] = [centroid]
        elif line[0] == 'H':
            current_seq = line[8][:line[8].find(separator)]
            centroid = line[9][:line[9].find(separator)]        
            otu_map[centroid].append(current_seq)
        elif line[0] == 'C':
            pass
    
    # Report basic informatin of the OTU map
    otu_map_parser = ParseOtuMap.otu_map_parser(otu_map)
    print("\n")
    print('The OTU map contains:')
    print('%i OTUs' % otu_map_parser.derep_count)
    print('%i sequences' % otu_map_parser.seqs_count)
    print('Max size=%i, Min size=%i, Ave size=%i' % (
        otu_map_parser.max_derep, otu_map_parser.min_derep, otu_map_parser.ave_derep))
    
    # Write the OTU map to a new file
    ParseOtuMap.write_otu_map(otu_map, output_file)
    print('Converted OTU map saved to %s' % output_file)

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
