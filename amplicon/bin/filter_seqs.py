# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 15:55:13 2015

Please feel free to contact with any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""
from __future__ import print_function
from __future__ import division
def main(Namespace):
    from lib import File_IO
    from lib import Seq_IO
    import argparse
    import textwrap
    import sys
    import time

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                    ------------------------
                                    By Zewei Song
                                    University of Minnesota
                                    Dept. Plant Pathology
                                    songzewei@outlook.com
                                    ------------------------'''), prog='fast.py -filter_seqs')
    parser.add_argument('-i', '--input', help='Name of the input file, can be FASTA or FASTQ')
    parser.add_argument('-o', '--output', help='Name of the output file')
    parser.add_argument('-maxN', help='Number of maximum ambiguous base')
    parser.add_argument('-maxhomop', help='Maximum length of homopolyer')
    args = parser.parse_args(Namespace)

    start = time.time()
    seq_file = args.input
    filtered_file = args.output

    print('Reading in %s ...' % seq_file)
    seqs = File_IO.read_seqs(seq_file)
    count_total = len(seqs)
    print('Found %d sequences.' % count_total)

    if len(seqs[0]) == 2:
        seqs_type = 'fasta'
    elif len(seqs[0]) == 4:
        seqs_type = 'fastq'
    else:
        print('This is not a corerct FASTA or FASTQ file, please check you file.')
        sys.exit()
    #
    checkN = False
    check_homop = False
    if args.maxN:
        maxN = int(args.maxN)
        checkN = True
        print('Maximum ambiguous base allowed: %d' % maxN)
    if args.maxhomop:
        maxhomop = int(args.maxhomop)
        check_homop = True
        print('Maximum length of homopolyer: %d' % maxhomop)
    else:
        pass
    checker = 0
    if checkN and check_homop:
        checker = 12
    elif checkN:
        checker = 1
    elif check_homop:
        checker = 2

    seqs_filtered = []
    count_pass = 0
    count_total = 0

    if checker == 12:
        for record in seqs:
            count_total += 1
            #sys.stderr.write('Processing %i sequence ...' % count_total + '\b' * 100,)
            current_record = record[1]
            if not Seq_IO.check_ambiguous(current_record, maxN):
                if not Seq_IO.check_homop(current_record, maxhomop + 1):
                    seqs_filtered.append(record)
                    count_pass += 1

    if checker == 1:
        for record in seqs:
            count_total += 1
            #sys.stderr.write('Processing %i sequence ...' % count_total + '\b' * 100,)
            current_record = record[1]
            if not Seq_IO.check_ambiguous(current_record, maxN):
                seqs_filtered.append(record)
                count_pass += 1

    if checker == 2:
        for record in seqs:
            count_total += 1
            #sys.stderr.write('Processing %i sequence ...' % count_total + '\b' * 100,)
            current_record = record[1]
            if not Seq_IO.check_homop(current_record, maxhomop + 1):
                seqs_filtered.append(record)
                count_pass += 1
    end = time.time()
    used_time = round(float(end - start), 2)
    print

    print ('Filtered %d sequences, %d (%s%%) passed. Used %s seconds.' % (
        count_total, count_pass, str(round(float(count_pass) / count_total, 1) * 100), str(used_time)))

    print('Writing to %s ...' % filtered_file)
    count = File_IO.write_seqs(seqs_filtered, filtered_file, checker=False, overwrite=True)
    print('Filtered sequences (%i seqs) store in %s' % (count, filtered_file))

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])