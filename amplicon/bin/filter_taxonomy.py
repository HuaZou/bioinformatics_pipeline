#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 05 12:18:11 2016



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
    parser.add_argument('-i', '--input', help='Input blast result, or any taxonomy result based on BLASTn output format 6.')
    parser.add_argument('-op', '--output_pass', default = 'pass.txt', help='Output for passed blast result.')
    parser.add_argument('-of', '--output_fail', help='Output for failed blast result.')
    parser.add_argument('-match_length', default = 0.9, help='Threshold for matched length (0.0 - 1.0).')
    parser.add_argument('-pident', default = 97, help='Threshold for pident (0 - 100).')
    parser.add_argument('-qlen_position', default = 3, help='Column position for qlen, length of query sequence in BLAST.')
    parser.add_argument('-length_position', default = 4, help='Column position for length, alignment length in BLAST.')
    parser.add_argument('-pident_position', default = 5, help='Column position for pident, Pident value in BLAST.')
    
    args = parser.parse_args(name_space)

    #%% Read in file
    input_file = args.input
    
    input_content = []
    with open(input_file, 'rU') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            input_content.append(line)
    
    print("{0} contains {1} records.".format(input_file, len(input_content)))
    
    #%% Filter
    pass_content = []
    fail_content = []
    
    match_l = float(args.match_length)
    pident = float(args.pident)
    
    print("Minimum matched length set to {0}, minimum Pident set to {1}.".format(match_l, pident))
    
    qlen_p = int(args.qlen_position) - 1
    length_p = int(args.length_position) - 1
    pident_p = int(args.pident_position) - 1
    
    for record in input_content:
        current_match_l = float(int(record[length_p]))/int(record[qlen_p])
        current_pident = float(record[pident_p])
        if current_match_l >= match_l and current_pident >= pident:
            pass_content.append(record)
        else:
            fail_content.append(record)
    
    #%% Output
    output_pass = args.output_pass
    with open(output_pass, 'w') as f:
        for record in pass_content:
            current_line = '\t'.join(record)
            f.write('%s\n' %current_line)
    
    print("{0} record passed and saved in {1}.".format(len(pass_content), output_pass))
    
    if args.output_fail:
        output_fail = args.output_fail
        with open(output_fail, 'w') as f:
            for record in fail_content:
                current_line = '\t'.join(record)
                f.write('%s\n' %current_line)
                print("{0} record failed and saved in {1}.".format%(len(fail_content), output_pass))

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])