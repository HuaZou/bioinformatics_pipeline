#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 15:22:31 2016

Pull a subset from the hybrid FAST map using a list of OTU names.

It will generate a FAST derep map and the corresponding derep sequences.

This subset can be used to perform OTU clustering.

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
    from lib import File_IO
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                    ------------------------
                                    By Zewei Song
                                    University of Minnesota
                                    Dept. Plant Pathology
                                    songzewei@outlook.com
                                    ------------------------'''), prog='fast.py -subset_fast_hybrid')
    parser.add_argument('-i', '--input', help='Input of FAST hybrid map.')
    parser.add_argument('-o', '--output', help='Output prefix for the derep map and sequence')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-otu_list', help='A list of OTU names seperated by ","')
    group.add_argument('-otu_file', help='A file contains a list of OTU names (no header)')
    args = parser.parse_args(name_space)
    
    print ('Subtracting a FAST hybrid map with provieded OTU names ...')    
    input_file = args.input
    output_derep = args.output + '.txt'
    output_fasta = args.output + '.fasta'
    
    otu_list = []
    if args.otu_list:
        otu_list = args.otu_list.split(',')
    elif args.otu_file:
        otu_list = []
        with open(args.otu_file) as f:
            for line in f:
                otu_list.append(line)
    print ('Found {0} OTU names.'.format(len(otu_list)))

    print ('Reading in the FAST hybrid map: {0} ...'.format(input_file))
    hybrid_map = ParseOtuMap.read_fast_output(input_file)
    fast_derep = {}
    for otu in otu_list:
        fast_derep.update(hybrid_map[otu]['sample'])
    ParseOtuMap.write_fast_output(fast_derep, output_derep)
    print ('A FAST derep map wrote to: {0}.'.format(output_derep))    
    
    derep_seq = []
    for key, value in fast_derep.items():
        current_seq = []
        derep_size = sum(value['sample'].values())
        seq_label = key + ';size=' + str(derep_size)
        current_seq = [derep_size, seq_label, value['seq']]
        derep_seq.append(current_seq) 
    
    derep_seq.sort(reverse=True)
    derep_seq = [i[1:] for i in derep_seq]

    count = File_IO.write_seqs(derep_seq, output_fasta, checker=False)
    print ('A dereplicated FASTA file wrote to {0}, containing {1} sequences with size annotation.'.format(output_fasta, count))
    print ('\n')

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
