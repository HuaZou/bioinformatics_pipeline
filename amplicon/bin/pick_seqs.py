# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 16:41:33 2015

Pick sequences from a given FASTA file using a OTU map or a name list, or randomly sampled the file.

For OTU map, the default method is to pick using name of each OTU.

Using -sequence to specify picking using sequence names within each OTU.

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
    from lib import random_subsample as rs
    from lib import ParseOtuMap
    from lib import Seq_IO
    from lib import File_IO
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                    ------------------------
                                    By Zewei Song
                                    University of Minnesota
                                    Dept. Plant Pathology
                                    songzewei@outlook.com
                                    ------------------------'''), prog='fast.py -pick_seqs')
    parser.add_argument('-i', '--input', help='Input FASTA file to be picked')
    parser.add_argument('-o', '--output', help='Name for output FASTA file.')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-map', help='OTU map file, will pick OTU names in default')
    group.add_argument('-name_list', help='File with names in separated lines')
    group.add_argument('-random_pick', help='Randomly pick the given number of sequences.')
    parser.add_argument('-sequence', action='store_true', help='Indicate to pick by sequence names instead of OTU names')
    parser.add_argument('-sizeout', action='store_true', help='Indicate to output size label.')
    args = parser.parse_args(name_space)
    
    input_fasta = args.input
    output_fasta = args.output
    pick_list = False
    
    print('\n')
    if args.map:
        pick_list = []
        if not args.sequence:
            otu_map = ParseOtuMap.read_otu_map(args.map)
            for key in otu_map:
                pick_list.append(key)
        elif args.sequence:
            otu_map = ParseOtuMap.read_otu_map(args.map)
            for key, value in otu_map.items():
                pick_list += value
        print('Picking sequences from the OTU map: %s.' % (args.map))
        print('Found %i names to be picked.' % len(pick_list))
    if args.name_list:
        pick_list = []
        with open(args.name_list, 'rU') as f:
            for line in f:
                pick_list.append(line.strip('\n'))
        print('Picking sequences from a OTU list.')
        print('Found %i names to be picked.' % len(pick_list))
    if args.random_pick:
        pick_size = int(args.random_pick)
        print('Randomly pick %i sequences.' % (pick_size))
        
    if pick_list == []:
        input_content = File_IO.read_seqs(input_fasta)
        print('Reaing in the original FASTA file: %s ...' % input_fasta)
        print('Randomly sampling %i sequences out of %i ...' %(pick_size, len(input_content)))
        seq_index = rs.generate_random_index(len(input_content), pick_size)
        sampled_content = []

        for index in seq_index:
            sampled_content.append(input_content[index])
        
        count = File_IO.write_seqs(sampled_content, output_fasta, checker=False, overwrite=True)
        print('Picked sequences wrote to %s.' % output_fasta)
    
    else:
        print('Reaing in the original FASTA file: %s ...' % input_fasta)
        input_content = File_IO.read_seqs(input_fasta)
        for record in input_content:
            record[0] = record[0].split(' ')[0]  # OTU name will be cut at the first space
            if record[0].find(';') != -1:
                record[0] = record[0][:record[0].find(';')] # Cut the label at the first ";"
        print('Indexing the original sequence file ...')
        input_dict = Seq_IO.make_dict(input_content)
        
        count_picked = 0
        count_missed = 0
        print('Search name list in the sequence file ...')
        picked_content = []
        
        if args.sizeout:
            print("Output size labels ...")
            size_list = []
            for record in pick_list:
                size_list.append([record, len(otu_map[record])])
            size_list = sorted(size_list, key=lambda x:x[1], reverse=True)
            for record in size_list:
                try:
                    new_label = record[0] + ';size=' + str(record[1]) 
                    picked_content.append([new_label, input_dict[record[0]][0]])
                    count_picked += 1
                except KeyError:
                    count_missed += 1
        
        else:
            for name in pick_list:
                try:
                    picked_content.append([name, input_dict[name][0]])
                    count_picked += 1
                except KeyError:
                    count_missed += 1
      
        print('Finished searching.')
        print('Original sequence=%i' % len(input_content))
        print('Input names=%i' % len(pick_list))
        print('Picked sequences=%i' % count_picked)
        print('Not found sequences=%i' % count_missed)
        
        print('Writing to a new FASTA file ...')
        count = File_IO.write_seqs(picked_content, output_fasta, checker=False, overwrite=True)
        print('Picked sequences wrote to %s.' % output_fasta)

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])