#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 11:44:11 2015

Assign taxonomy information from  a BLAST output to the OTU table.

Now it only accept following output from BLAST:
@ blastn -db blast\\unite_02.03.2015 -query raw.qc.fasta_rep_set.fasta -max_target_seqs 1 -outfmt "6 qseqid stitle qlen length pident evalue" -out rep.otu.txt

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
    from lib import ParseOtuTable
    import argparse
    import textwrap
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                        ------------------------
                                        By Zewei Song
                                        University of Minnesota
                                        Dept. Plant Pathology
                                        songzewei@outlook.com
                                        ------------------------'''), prog = '-assign_taxonomy')
    parser.add_argument("-otu", help="Input OTU table")
    parser.add_argument("-tax", help="Taxonomy search result in blast6out format: query+target+ql+pairs+id")
    parser.add_argument("-o", "--output", help="Output OTU table")
    parser.add_argument('-scores', action='store_true', help='Indicate to output matching scores following the taxonomy column.')
    parser.add_argument('-keep_all', action='store_true', help='Indicate to keep all OTUs including low confidencial ones.')
    parser.add_argument('-min_match', default = 0.0, help='Threshod for minimum match length to keep in OTU table.')
    parser.add_argument('-min_pident', default = 0.0, help='Threshod for minimum pident to keep in OTU table.')
    parser.add_argument('-match_length', default = 0.0, help='Threshold for matched length that can be assigned a taxa.')
    parser.add_argument('-pident', default = 97, help='Threshold for pident that can be assigned a taxa.')
    
    args = parser.parse_args(name_space)
    
    input_otu = args.otu
    input_tax = args.tax
    output_otu = args.output
    
    otu_table = ParseOtuTable.parser_otu_table(input_otu)
    
    # Read in BLAST 6 style taxonomy result and store in a dictionary (the out has to be arrange as query+target+ql+pairs+id )
    input_content = []
    with open(input_tax, 'rU') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            input_content.append(line)
            
   
    taxonomy = {}
    for line in input_content:
        taxonomy[line[0]] = line[1:]
    
    # Filther the taxonomy result using the minimal threshold (These records are likely to be in the targeted group, such as fungi).
    taxonomy_keep = {}
    low_confident_count = 0
    for key, value in taxonomy.items():
        match_score = float(int(value[2]))/int(value[1])
        
        pident_score = float(value[3])
        #print([match_score, pident_score, args.min_match, args.min_pident])
        if match_score >= float(args.min_match) and pident_score >= float(args.min_pident):
            taxonomy_keep[key] = value
        else:
            low_confident_count += 1
    
    taxonomy_assign = {}
    no_hit_count = 0
    for key, value in taxonomy_keep.items():
        match_score = float(int(value[2]))/int(value[1])
        pident_score = float(value[3])
        if match_score >= float(args.match_length) and pident_score >= float(args.pident):
            taxonomy_assign[key] = value
        else:
            no_hit_count += 1
    
    print('{0} OTUs can be aligned to the reference database.'.format(len(taxonomy)))
    print('{0} assignments have match length >= {1}, and Pident >= {2}. Other assignments are marked as "below_minimum_threshold".'.format(len(taxonomy_keep), args.min_match, args.min_pident))
    print('{0} assignments have match length >= {1}, and Pident >= {2}. Other assignments are marked as "below_assignment_threshold".'.format(len(taxonomy_assign), args.match_length, args.pident))

    otu_matrix = otu_table.species_matrix
    otu_matrix_assigned = []
    
    low_confident_otu = 0
    no_hit_otu = 0
    hit_otu = 0
    
    for line in otu_matrix:
        #print(taxonomy_keep[line[0]])
        try:
            line += taxonomy_assign[line[0]] + ['above_assignment_threshold']
            otu_matrix_assigned.append(line)
            hit_otu += 1
        except KeyError:
            try:
                line += taxonomy_keep[line[0]] + ['below_assignment_threshold']
                otu_matrix_assigned.append(line)
                no_hit_otu += 1
            except KeyError:
                if args.keep_all:
                    try:
                        line += taxonomy[line[0]] + ['below_minimum_threshold']
                        otu_matrix_assigned.append(line)
                        low_confident_otu += 1
                    except KeyError:
                        line += ['nan','nan','nan','nan','below_minimum_threshold']
                        otu_matrix_assigned.append(line)
                        low_confident_otu += 1
                else:
                    pass

    sample_id = otu_table.sample_id
    sample_id = ['OTU_ID'] + sample_id + ['taxonomy','Query_length','Query_aligned_length','Pident','Threshold']
    
    # Write the new OTU table
    with open(output_otu, 'w') as f:
        f.write('%s\n' % '\t'.join(sample_id))
        for line in otu_matrix_assigned:
            line = [str(x) for x in line]
            f.write('%s\n' % '\t'.join(line))
    
    print('The new OTU table was saved in {0}.'.format(output_otu))
    
if __name__ == '__main__':
    import sys
    main(sys.argv[1:])