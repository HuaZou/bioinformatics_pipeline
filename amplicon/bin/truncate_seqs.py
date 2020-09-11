#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 08 13:35:44 2015

This script will truncate the given FASTA file to a fixed length.
Any sequences shorter than the provided length will be discard.
Any sequences longer than the provided length will be cut to the length.

Please feel free to contact me for any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""
def main(name_space):
    import argparse
    import textwrap
    from lib import File_IO

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                        ------------------------
                                        By Zewei Song
                                        University of Minnesota
                                        Dept. Plant Pathology
                                        songzewei@outlook.com
                                        ------------------------'''), prog = 'fast.py -truncate_seqs')

    parser.add_argument("-i", "--input", help="Name of the input FASTA file.")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-fixed_length', help='A fixed length to cut on all sequences.')
    group.add_argument('-slice', help='Slice size to cut from head and tail of each sequence in the format of "head,tail".')
    parser.add_argument('-sliced_out', action='store_true', help='Indicate to output sliced sequences.')
    parser.add_argument("-o", "--output", help="Name of the output file.")
    args = parser.parse_args(name_space)

    if args.fixed_length:
        truncate_length = int(args.fixed_length)
        sequences = File_IO.read_seqs(args.input)
        count = len(sequences)
        print("Reading in %s ..." % args.input)
        print("%s contains %i records." % (args.input, count))
        print("Cutting sequences to a fixed length: %i ..." % truncate_length)

        count_fail = 0
        with open(args.output, 'w') as f:
            for record in sequences:
                if len(record[1]) >= truncate_length:
                    if len(record) == 2:
                        f.write('>%s\n' % record[0])
                        f.write('%s\n' % record[1][:truncate_length])
                    elif len(record) == 4:
                        f.write('@%s\n' % record[0])
                        f.write('%s\n' % record[1][:truncate_length])
                        f.write('%s\n' % record[2])
                        f.write('%s\n' % record[3][:truncate_length])
                else:
                    count_fail += 1
        print("%i sequences were cut to %i and save in %s." % (count - count_fail, truncate_length, args.output))

    if args.slice:
        slice_window = args.slice.split(',')
        head = int(slice_window[0])
        tail = int(slice_window[1])
        sequences = File_IO.read_seqs(args.input)
        count = len(sequences)
        print("Reading in %s ..." % args.input)
        print("%s contains %i records." % (args.input, count))
        print("Slicing %i bp from the head and %i bp from the tail ..." % (head, tail))

        count_fail = 0
        with open(args.output, 'w') as f:
            for record in sequences:
                seq_len = len(record[1])
                if seq_len > head + tail:
                    if len(record) == 2:
                        f.write('>%s\n' % record[0])
                        f.write('%s\n' % record[1][head:(seq_len - tail)])
                    elif len(record) == 4:
                        f.write('@%s\n' % record[0])
                        f.write('%s\n' % record[1][head:(seq_len - tail)])
                        f.write('%s\n' % record[2])
                        f.write('%s\n' % record[3][head:(seq_len - tail)])
                else:
                    count_fail += 1
        print("%i sequences were sliced and save in %s." % (count - count_fail, args.output))

        if args.sliced_out:
            if head > 0:
                head_output = 'head.' + args.output
                with open(head_output, 'wb') as f:
                    for record in sequences:
                        if len(record) == 2:
                            f.write('>%s\n' % record[0])
                            f.write('%s\n' % record[1][:head])
                        elif len(record) == 4:
                            f.write('@%s\n' % record[0])
                            f.write('%s\n' % record[1][:head])
                            f.write('%s\n' % record[2])
                            f.write('%s\n' % record[3][:head])
                print('The sliced head sequences wrote to %s.' % (head_output))

            if tail > 0:
                tail_output = 'tail.' + args.output
                with open(tail_output, 'wb') as f:
                    for record in sequences:
                        seq_len = len(record[1])
                        if len(record) == 2:
                            f.write('>%s\n' % record[0])
                            f.write('%s\n' % record[1][(seq_len - tail):])
                        elif len(record) == 4:
                            f.write('@%s\n' % record[0])
                            f.write('%s\n' % record[1][(seq_len - tail):])
                            f.write('%s\n' % record[2])
                            f.write('%s\n' % record[3][(seq_len - tail):])
                print('The sliced tail sequences wrote to %s.' % (tail_output))

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])