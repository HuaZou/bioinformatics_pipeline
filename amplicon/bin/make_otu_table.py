# -*- coding: utf-8 -*-
"""
Create on April 20th, 2015.

Create an OTU table from a Qiime style OTU map. The OTUs will be sorted by their total abundance.

Please feel free to contact me with any question.
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
    from lib import File_IO
    #import sys
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                    ------------------------
                                    By Zewei Song
                                    University of Minnesota
                                    Dept. Plant Pathology
                                    songzewei@outlook.com
                                    ------------------------'''), prog='fast.py -make_otu_table')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-qiime_map', help='The Qiime style OTU map.')
    group.add_argument('-fast_map', help='The FAST hybrid OTU map.')
    parser.add_argument('-o', '--output', help='Output OTU table.')
    parser.add_argument('-rep', help='Indicate to output a representative sequnce if using FAST method.')
    args = parser.parse_args(name_space)
    
    if args.qiime_map != None:
        input_file = args.qiime_map
        method = 'qiime'
    elif args.fast_map != None:
        input_file = args.fast_map
        method = 'fast'
        
        if args.rep != None:
            output_seq_file = args.rep
            
    output_file = args.output
    
    # Parse OTU map into OTU table dictionary
    
    # Use Qiime style
    if method == 'qiime':
        print('Reading the Qiime style OTU map: {0} ...'.format(input_file))
        otu_map = ParseOtuMap.read_otu_map(input_file)
        sample_list = []
        otu_table_dict = {}
        for key, value in otu_map.items():
            otu_table_dict[key] = {}
            for sample in value:
                treatment = sample[:sample.find('_')]
                if treatment not in sample_list:
                    sample_list.append(treatment)
                try:
                    otu_table_dict[key][treatment] += 1
                except KeyError:
                    otu_table_dict[key][treatment] = 1
    
    # Use FAST style
    if method == 'fast':
        print('Reading the FAST hybrid OTU map: {0} ...'.format(input_file))
        otu_map = ParseOtuMap.read_fast_output(input_file)
#        sample_list = []
#        otu_table_dict = {}
#        for otu, value in otu_map.items():
#            otu_table_dict[otu] = {}
#            for derep_unit, derep_value in value['sample'].items():
#                for sample, abundance in derep_value['sample'].items():
#                    treatment = sample
#                    if treatment not in sample_list:
#                        sample_list.append(treatment)
#                    try:
#                        otu_table_dict[otu][treatment] += abundance
#                    except KeyError:
#                        otu_table_dict[otu][treatment] = abundance
        otu_map_parser = ParseOtuMap.fast_output_parser(otu_map)
        sample_list, otu_table_dict = otu_map_parser.parse_otu_table()
        
        if args.rep != None:
            temp_content = otu_map_parser.get_seqs()
            rep_seq = []
            for item in temp_content:
                rep_seq.append(item[:2])
            rep_seq_count = File_IO.write_seqs(rep_seq, output_seq_file, checker=False, overwrite=True) 
            print('{0} OTUs were wrote to {1}.'.format(rep_seq_count, output_seq_file))
    
    # Convert OTU table dictionary to table
#    otu_abundance = {}
#    for sample in sample_list:
#        otu_abundance[sample] = 0  # Set initial abundance to zero as place holder
    sample_list.sort()
    otu_table = []
    for key, value in otu_table_dict.items():
        current_otu = [key]
        for sample in sample_list:
            try:
                current_otu.append(value[sample])
            except KeyError:
                current_otu.append(0)
        otu_table.append(current_otu)
    otu_table.sort(key=lambda x: sum(map(int, x[1:])), reverse=True)
    otu_table = [['OTU_ID'] + sample_list] + otu_table
    
    # Write OTU table to a new file
    sample_list = ['OTU_ID'] + sample_list
    with open(output_file, 'w') as f:
        for line in otu_table:
            line = [str(i) for i in line]
            f.write('%s\n' % '\t'.join(line))
    
    print('OTU table with {0} samples was saved in {1}.'.format(len(sample_list)-1, output_file))

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
