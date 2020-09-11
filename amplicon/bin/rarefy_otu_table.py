# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 19:23:08 2015

Randomly sample (rarefying) the OTU table. If -iter is specified larger than 1, each sample
will be repeated, and the sample with a richness in the medium will be kept. Pick the one in 
the center of NMDs cloud could be a better method.

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
    from lib import random_subsample as rs
    from lib import ParseOtuTable
    from lib import File_IO
    import argparse
    import textwrap
    import time
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                    ------------------------
                                    By Zewei Song
                                    University of Minnesota
                                    Dept. Plant Pathology
                                    songzewei@outlook.com
                                    ------------------------'''), prog='fast.py -rarefy_otu_table')
    parser.add_argument('-otu', help='Input OTU table')
    parser.add_argument('-o', '--output', help='Output OTU table')
    parser.add_argument('-d', '--depth', help='Sampling depth for each sample')
    parser.add_argument('-iter', default=1, help='Iteration time for each sample')
    parser.add_argument('-thread', default=1, help='Number of threads')
    parser.add_argument('-keep_all', action='store_true', help='Indicate to keep all samples')
    parser.add_argument('-meta_column', default='taxonomy', help='Name of the first meta data')
    args = parser.parse_args(name_space)
    
    input_otu = args.otu
    
    iter_num = int(args.iter)
    thread = int(args.thread)
    meta_col = args.meta_column
    if args.output:
        output_otu = args.output
    else:
        output_otu = File_IO.name_file(input_otu, '', 'rare')
    

    otu_table = ParseOtuTable.parser_otu_table(input_otu, meta_col=meta_col)
    input_sample = otu_table.sample_matrix
    start = time.time()

    print('Input OTU table: %s' % input_otu)
    if args.depth:
        depth = int(args.depth)
        print('Sampling depth: %i' % depth)
    else:
        depth = min([sum(i[1:]) for i in input_sample])
        print('Sampling depth set to the minimum abundance: %i' % depth)
    print('Iteration time for each OTU: %i' % iter_num)
    print('Threads number: %i' % thread)
    print('Reading in the OTU table ...')

    if args.keep_all:
        count = 0
        for line in input_sample:
            if sum(line[1:]) < depth:
                count += 1
        print('Found %i samples in the OTU table.' % len(input_sample))
        print('%i samples has total abundance less than the sampling depth, but will be kept in the output.' % count)
    else:
        temp = []
        count = 0
        for line in input_sample:
            if sum(line[1:]) >= depth:
                temp.append(line)
            else:
                count += 1
        input_sample = temp
        print('Found %i samples in the OTU table.' % len(input_sample))
        print('%i samples has total abundance less than the sampling depth, and will be excluded.' % (count))

    otu_id = otu_table.species_id
    otu_table_rarefied = [['OTU_ID'] + otu_id + otu_table.meta_id]

    for sample in input_sample:
        print('Rarefying %s ...' % sample[0])
        repeat_sample = rs.repeat_rarefaction_parallel(sample[1:], depth, iter_num, processor=thread)
        repeat_sample.sort(key=lambda x: sum(i > 0 for i in x))
        repeat_sample = [sample[0]] + repeat_sample[int(iter_num / 2)]  # Pick the rarefied sample with the average richness
        otu_table_rarefied.append(repeat_sample[:])

    otu_table_rarefied = [list(i) for i in zip(*otu_table_rarefied)]

    # Add meta data
    meta_data = otu_table.meta_dict()
    for line in otu_table_rarefied[1:]:
        for key in otu_table.meta_id:
            line.append(meta_data[key][line[0]])
    for key in otu_table.meta_id:
        otu_table_rarefied[0].append(key)

    with open(output_otu, 'w') as f:
        for line in otu_table_rarefied:
            line = [str(i) for i in line]
            f.write('%s\n' % '\t'.join(line))
    print('Rarefied OTU table saved in %s.' % output_otu)
    end = time.time()
    used_time = round(float(end - start), 2)
    time_per_sample = round(used_time / len(input_sample), 2)
    print('Total time used: %s seconds (%s seconds per sample)' % (str(used_time), str(time_per_sample)))

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])