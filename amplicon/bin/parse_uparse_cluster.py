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
                                    ------------------------'''), prog='fast.py -parse_uparse_cluster')
    parser.add_argument('-i', '--input', help='Input uparse denovo cluster output file')
    parser.add_argument('-o', '--output', help='Output Qiime style OTU map')
    parser.add_argument('-separator', default=';', help='Separator between sample name and size annotation')
    args = parser.parse_args(name_space)
    input_file = args.input
    output_file = args.output
    separator = args.separator
    
    # Read in the uparse file
    uparse_content = []
    count_uparse = 0
    with open(input_file, 'rU')  as f:
        for line in f:
            temp = line.strip('\n').split('\t')
            temp[0] = temp[0][:temp[0].find(separator)]  # Remove size annotation and other extra labels
            uparse_content.append(temp)
            count_uparse += 1
            #sys.stderr.write('%i line in the uparse file ...' % count_uparse + '\b' * 100,)
    
    # Convert uparse output into Qiime style OTU map
    otu_map = {}
    count_chimera = 0
    for line in uparse_content:
        if line[1] == 'otu':
            otu_map[line[0]] = [line[0]]
        elif line[1] == 'match':
            otu_map[line[4]].append(line[0])
        elif line[1] == 'chimera':
            count_chimera += 1
    
    # Report basic informatin of the OTU map
    otu_map_parser = ParseOtuMap.otu_map_parser(otu_map)
    print '\n'
    print 'The OTU map contains:'
    print '%i OTUs' % otu_map_parser.derep_count
    print '%i sequences' % otu_map_parser.seqs_count
    print 'Max size=%i, Min size=%i, Ave size=%i' % (
        otu_map_parser.max_derep, otu_map_parser.min_derep, otu_map_parser.ave_derep)
    
    # Write the OTU map to a new file
    ParseOtuMap.write_otu_map(otu_map, output_file)
    print 'Converted OTU map saved to %s' % output_file

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])