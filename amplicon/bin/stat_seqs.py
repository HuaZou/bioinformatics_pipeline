# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 09:43:19 2015

Generate a report file about input FASTA file, include:
1) Total number of sequences.
2) Ambiguous bases count distribution against sequence length.
3) Max length of homopolyer distribution against sequence length (include All bases, and each base separated).
4) Number of sequences distribution against sequence length.

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
    from lib import File_IO
    from lib import Seq_IO
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
                                        ------------------------'''), prog='fast.py -stat_seqs')
    parser.add_argument("-i", "--input", help="Name of the input sequence file")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-o", "--output", default='report.txt', help="Specify a report file for output")
    args = parser.parse_args(name_space)
    seq_file = args.input
    report_file = args.output

    print("Reading in %s ..." % seq_file)
    seq_content = File_IO.read_seqs(seq_file)
    print('Found %d sequences in this curernt file, analyzing ...' % len(seq_content))
    start = time.time()
    #count = 0
    seq_length = {}
    seq_ambiguous = {}
    seq_homop = {'All': {}, 'A': {}, 'T': {}, 'C': {}, 'G': {}}
    seq_total_bases = {'A': 0, 'T': 0, 'C': 0, 'G': 0}

    for record in seq_content:
        temp_seq = record[1]
        temp_length = len(temp_seq)
        try:
            seq_length[temp_length] += 1
        except KeyError:
            seq_length[temp_length] = 1

        temp_ambiguous = Seq_IO.count_ambiguous(temp_seq)
        try:
            seq_ambiguous[temp_ambiguous] += 1
        except KeyError:
            seq_ambiguous[temp_ambiguous] = 1

        temp_homop = Seq_IO.count_homop(temp_seq)
        for base in temp_homop:
            temp_max_length = temp_homop[base]
            try:
                seq_homop[base][temp_max_length] += 1
            except KeyError:
                seq_homop[base][temp_max_length] = 1

        temp_bases_count = Seq_IO.count_bases(temp_seq)
        for key in seq_total_bases:
            seq_total_bases[key] += temp_bases_count[key]

    end = time.time()
    used_time = round(end - start, 2)
    print('Finished analyzing, used %s second, printing report ...' % str(used_time))

    # Make all four homopolyer distribution the same length
    all_bases_homop_len = []
    for base in seq_homop:
        for base_length in seq_homop[base]:
            if base_length not in all_bases_homop_len:
                all_bases_homop_len.append(base_length)
    for base in seq_homop:
        for length in all_bases_homop_len:
            try:
                seq_homop[base][length]
            except KeyError:
                seq_homop[base][length] = 0

    # Get all possible sequence length
    all_length = [i for i in seq_length]
    all_length.sort()

    #%% Write the report
    with open(report_file, 'w') as report:
        report.write('Report:\t%s\n' % report_file)
        report.write('Total number of sequence:\t%d\n\n' % len(seq_content))

        report.write('#' * 100 + '\n')
        report.write('Ambiguous base distribution:\nNumber of N\tNumber of sequences\n')
        for key in sorted(seq_ambiguous.keys()):
            report.write('%d\t%d\n' % (key, seq_ambiguous[key]))

        report.write('#' * 100 + '\n')
        report.write('Max homopolymer distribution:\nMax homopolyer length\tAll bases\tA\tT\tC\tG\n')
        for key in sorted(seq_homop['A'].keys()):
            report.write('%s\t%d\t%d\t%d\t%d\t%d\t\n' % (
                key, seq_homop['All'][key], seq_homop['A'][key], seq_homop['T'][key], seq_homop['C'][key],
                seq_homop['G'][key]))

        report.write('#' * 100 + '\n')
        report.write(
            'Length distribution:\tMaximum length:\t%d\tMinimum length:\t%d\n' % (max(all_length), min(all_length)))
        report.write('Length\tNumber of sequences\n')
        for key in sorted(seq_length.keys(), reverse=True):
            report.write('%d\t%d\n' % (key, seq_length[key]))
    print('Report on %s can be found in %s.' % (seq_file, report_file))

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])