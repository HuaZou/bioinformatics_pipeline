#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 15:51:12 2016



Please feel free to contact me for any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
www.songzewei.org
"""
from __future__ import print_function
def main(Name_space):
    import argparse
    import textwrap
    from lib import ParseOtuTable

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                    ------------------------
                                    By Zewei Song
                                    University of Minnesota
                                    Dept. Plant Pathology
                                    songzewei@outlook.com
                                    ------------------------'''), prog='fast.py -summary_otu_table')
    parser.add_argument('-otu', help='Name of the input OTU table.')
    parser.add_argument('-o', '--output', default='otu_report.txt', help='Name of the output report.')
    args = parser.parse_args(Name_space)
    
    input_file = args.otu
    output_file = args.output
    
    otu_table = ParseOtuTable.parser_otu_table(input_file)

    sample_id = otu_table.sample_id
    
    print('{0} has {1} samples.'.format(input_file, len(sample_id)))    
    
    sample_dict = otu_table.sample_dict()    
    report_dict = {}    
    for sample in sample_id:
        report_dict[sample] = {'depth':0, 'richness':0}
    
    for sample in sample_id:    
        report_dict[sample]['depth'] = sum(sample_dict[sample].values())
        report_dict[sample]['richness'] = len(sample_dict[sample].values())
    
    with open(output_file, 'w') as f:
        f.write('Sample\tDepth\tRichness\n')
        for sample in sample_id:
            line = ''
            line = sample + '\t' + str(report_dict[sample]['depth']) + '\t' + str(report_dict[sample]['richness'])
            f.write('%s\n' % line)
    
    print('A summary of {0} was wrote to {1}'.format(input_file, output_file))

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])