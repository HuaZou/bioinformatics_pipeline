#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 17:00:27 2016

This is the main script for the FAST (Fungal Amplicon Sequencing Toolbox) package.

Please feel free to contact me for any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
www.songzewei.org
"""
from __future__ import print_function
from __future__ import division
def main():
    import argparse
    import textwrap
    import sys
    #import importlib
    
    parser = argparse.ArgumentParser(prog='fast.py -function', add_help=False)
    
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('-add_labels', action = "store_true")
    group.add_argument('-add_seqs_size', action = "store_true")
    group.add_argument('-assign_taxonomy', action = 'store_true')
    group.add_argument('-combine_fast_map', action = "store_true")
    group.add_argument('-convert_fastq', action = 'store_true')
    group.add_argument('-correct_fasta', action = 'store_true')
    group.add_argument('-count_seqs', action = 'store_true')
    group.add_argument('-dereplicate', action = "store_true")
    group.add_argument('-document', action='store_true')
    group.add_argument('-filter_database', action = 'store_true')
    group.add_argument('-filter_otu_map', action = "store_true")
    group.add_argument('-filter_seqs', action = "store_true")
    group.add_argument('-filter_taxonomy', action = "store_true")
    group.add_argument('-generate_fast_map', action = "store_true")
    group.add_argument('-generate_mapping', action = "store_true")
    group.add_argument('-make_otu_table', action = 'store_true')
    group.add_argument('-merge_otu_maps', action = 'store_true')
    group.add_argument('-merge_seqs', action = "store_true")
    group.add_argument('-nucl_freq', action = "store_true")
    group.add_argument('-otu_deconstruct', action = 'store_true')
    group.add_argument('-otu_map_info', action = 'store_true')
    group.add_argument('-parse_uc_cluster', action = "store_true")
    group.add_argument('-parse_uparse_cluster', action = "store_true")
    group.add_argument('-pick_seqs', action = 'store_true')
    group.add_argument('-random_dataset', action = 'store_true')
    group.add_argument('-random_pick', action = 'store_true')
    group.add_argument('-rarefy_otu_table', action = 'store_true')
    group.add_argument('-rename_otu_map', action = "store_true")
    group.add_argument('-reorder_otu_table', action = "store_true")
    group.add_argument('-split_taxa', action = 'store_true')
    group.add_argument('-stat_seqs', action = 'store_true')
    group.add_argument('-subset_fast_hybrid', action = 'store_true')
    group.add_argument('-substract_controls', action = 'store_true')
    group.add_argument('-summary_otu_table', action = 'store_true')
    group.add_argument('-truncate_seqs', action = 'store_true')
    group.add_argument('-unmerge_seqs', action = 'store_true')
    
    if len(sys.argv) <= 1:
        print("This is the helping document:")
        sys.exit()
    
    else:
        args = parser.parse_args([sys.argv[1]])
#        function_name = False
#        for option in args.__dict__:
#            if args.__dict__[option]:
#                function_name = option
       
        sub_args = sys.argv[2:]
    
    if args.document:
        print("This is the helping document:")
    
    #    else:
    #        function = importlib.import_module(function_name)
    #        function.main(sub_args)
    
    #    elif args.generate_mapping:
    #        import generate_mapping
    #        generate_mapping.main(sub_args)
        
    #    else:
    #        function = __import__(function_name)
    #        function.main(sub_args)
        
    if args.add_labels:
        import add_labels as function
        function.main(sub_args)
    
    if args.add_seqs_size:
        import add_seqs_size as function
        function.main(sub_args)
    
    if args.assign_taxonomy:
        import assign_taxonomy as function
        function.main(sub_args)
    
    if args.combine_fast_map:
        import combine_fast_map as function
        function.main(sub_args)
    
    if args.convert_fastq:
        import convert_fastq as function
        function.main(sub_args)
    
    if args.correct_fasta:
        import correct_fasta as function
        function.main(sub_args)
    
    if args.count_seqs:
        import count_seqs as function
        function.main(sub_args)
    
    if args.dereplicate:
        import dereplicate as function
        function.main(sub_args)
    
    if args.filter_database:
        import filter_database as function
        function.main(sub_args)
    
    if args.filter_otu_map:
        import filter_otu_map as function
        function.main(sub_args)
     
    if args.filter_seqs:
        import filter_seqs as function
        function.main(sub_args)

    if args.filter_taxonomy:
        import filter_taxonomy as function
        function.main(sub_args)
    
    if args.generate_fast_map:
        import generate_fast_map as function
        function.main(sub_args)
    
    if args.generate_mapping:
        import generate_mapping as function
        function.main(sub_args)
    
    if args.make_otu_table:
        import make_otu_table as function
        function.main(sub_args)
    
    if args.merge_otu_maps:
        import merge_otu_maps as function
        function.main(sub_args)
    
    if args.merge_seqs:
        import merge_seqs as function
        function.main(sub_args)
        
    if args.nucl_freq:
        import nucl_freq as function
        function.main(sub_args)
    
    if args.otu_deconstruct:
        import otu_deconstruct as function
        function.main(sub_args)
    
    if args.otu_map_info:
        import otu_map_info as function
        function.main(sub_args)
    
    if args.parse_uc_cluster:
        import parse_uc_cluster as function
        function.main(sub_args)
    
    if args.parse_uparse_cluster:
        import parse_uparse_cluster as function
        function.main(sub_args)
    
    if args.pick_seqs:
        import pick_seqs as function
        function.main(sub_args)
    
    if args.random_dataset:
        import random_dataset as function
        function.main(sub_args)
    
    if args.random_pick:
        import random_pick as function
        function.main(sub_args)    
    
    if args.rarefy_otu_table:
        import rarefy_otu_table as function
        function.main(sub_args)
    
    if args.rename_otu_map:
        import rename_otu_map as function
        function.main(sub_args)

    if args.reorder_otu_table:
        import reorder_otu_table as function
        function.main(sub_args)
    
    if args.split_taxa:
        import split_taxa as function
        function.main(sub_args)

    if args.stat_seqs:
        import stat_seqs as function
        function.main(sub_args)

    if args.subset_fast_hybrid:
        import subset_fast_hybrid as function
        function.main(sub_args)
    
    if args.substract_controls:
        import substract_controls as function
        function.main(sub_args)

    if args.summary_otu_table:
        import summary_otu_table as function
        function.main(sub_args)
    
    if args.truncate_seqs:
        import truncate_seqs as function
        function.main(sub_args)

    if args.unmerge_seqs:
        import unmerge_seqs as function
        function.main(sub_args)
        
if __name__ == '__main__':
    main()