# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 13:32:33 2015

This script will remove the abundance value in extract control from all samples.
Samples can have different extract control, by providing the names in a file:
Example:

OTU_ID Control_name
Sample_1 Extract_run1
Sample_2 Extract_run1
Sample_3 Extract_run1
Sample_4 Extract_run2
Sample_5 Extract_run2
Sample_6 Extract_run2

Save the name list in a tab-delimited text file.

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
def main(name_space):
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
                                        ------------------------'''), prog = 'fast.py -substract_controls')

    parser.add_argument("-otu", help="Name of the input OTU table (tab delimited file).")
    parser.add_argument("-control", help="Name of the extract control list.")
    parser.add_argument("-o", "--output", help="Name of the output OTU table")
    parser.add_argument("-keep_zero", action="store_true", help="Indicate that all OTU that sum to zero should be kept.")
    args = parser.parse_args(name_space)

    # Parse the extract contrl name list
    control_list = args.control
    name_dict = {} # corresponding extract control for all samples
    with open(control_list, 'rU') as f:
        temp_list = []
        for line in f:
            line = line.strip('\n').split('\t')
            temp_list.append(line)

    for item in temp_list[1:]:
        name_dict[item[0]] = item[1]

    controls = sorted(set(name_dict.values()))
    print("Found %i extract control provided by the user:" % len(controls))
    for item in controls:
        print(item)

    # Read in the OTU table
    print("Reading in OTU table %s ..." % args.otu)
    otu_table = ParseOtuTable.parser_otu_table(args.otu)
    sample_id = otu_table.sample_id
    otu_id = otu_table.species_id
    meta_id = otu_table.meta_id
    print("Current OTU table contains %i samples, %i OTUs, and %i meta columns." % (len(sample_id), len(otu_id), len(meta_id)))
    sample_dict = otu_table.sample_dict()

    # Substract extract control abundances
    print("Substracting extract control from corresponding samples ...")
    sample_dict_new = {}
    for sample in name_dict.keys():
        sample_dict_new[sample] = {}
        for otu in sample_dict[sample].keys():
            try:
                sample_dict[name_dict[sample]][otu]
                control_abundance = sample_dict[name_dict[sample]][otu]
                sample_abundance = sample_dict[sample][otu]
                if sample_abundance >= control_abundance:
                    new_abundance = sample_abundance - control_abundance
                else:
                    new_abundance = 0
            except KeyError:
                new_abundance = sample_dict[sample][otu]
            sample_dict_new[sample][otu] = new_abundance

    #%% Output the new OTU table
    if args.keep_zero:
        remove_zero = False
        print("OTUs with 0 abundance will be kept in the new OTU table.")
    else:
        remove_zero = True
        print("OTUs with 0 abundance will be discarded in the new OTU table.")

    ParseOtuTable.write_sample_dict(sample_dict_new, otu_table.meta_dict(), otu_table.species_id, 'temp.txt')
    new_table = ParseOtuTable.parser_otu_table('temp.txt')
    import os
    os.remove('temp.txt')
    ParseOtuTable.sorted_table(new_table, new_file_path=args.output, remove_zero=remove_zero)

    new_table = ParseOtuTable.parser_otu_table(args.output)
    sample_id = new_table.sample_id
    otu_id = new_table.species_id
    meta_id = new_table.meta_id
    print("New OTU table saved in %s." % args.output)
    print("New OTU table contains %i samples, %i OTUs, and %i meta columns." % (len(sample_id), len(otu_id), len(meta_id)))

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])