#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 20:36:50 2016

Calculte the frequency of nucleotide from the start (and end) of sequences.

Please feel free to contact me with any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""
from __future__ import print_function


def main(Namespace):
    import argparse
    import textwrap
    from lib import File_IO
    from lib import Seq_IO

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                    ------------------------
                                    By Zewei Song
                                    University of Minnesota
                                    Dept. Plant Pathology
                                    songzewei@outlook.com
                                    ------------------------'''), prog='fast.py -nucl_freq')
    parser.add_argument('-i', '--input', help='Name of the input FASTA or FASTQ file.')
    parser.add_argument('-o', '--output', default='nucl_report.txt', help='Name of the reporting file.')
    parser.add_argument('-tail', action='store_true', help='Indicate to also count from the tail of the sequences.')
    args = parser.parse_args(Namespace)

    input_file = args.input
    output_file = args.output
    tail_indicator = args.tail

    print('Reading in file: {0} ...'.format(input_file))
    input_seq = File_IO.read_seqs(input_file)
    print('The file contains {0} sequences.'.format(len(input_seq)))
    if tail_indicator:
        print('Counting nucleotide frequencies from both ends of all sequences...')
    else:
        print('Counting nucleotide frequencies from the head of all sequences...')

    nucl_freq, unidentified_count = Seq_IO.nucl_freq(input_seq, tail = tail_indicator)

    nucl_list = ['A','T','C','G','N']

    # Output the counting result
    header = 'Position\tA\tT\tC\tG\tN\tMost frequent\tFrequency'
    output_content = [header]
    for pos in range(len(nucl_freq)):

        # Get the most frequent nucleotide at this position:
        temp_list = []
        most_freq_nucl = ""
        for nucl in nucl_list:
            temp_list.append([nucl_freq[pos][nucl],nucl])
        temp_list.sort(reverse=True)
        most_freq_nucl = temp_list[0][1]
        sum_nucl_count = sum(nucl_freq[pos].values())
        most_freq_nucl_freq = float(temp_list[0][0]) / sum_nucl_count # calculate the frequency of this nucleotide

        # Get the output for current position
        current_line = []
        current_line = [str(pos + 1)]
        for nucl in nucl_list:
            current_line.append(str(nucl_freq[pos][nucl]))
        current_line.append(most_freq_nucl)
        current_line.append(str(most_freq_nucl_freq))
        current_line = '\t'.join(current_line)
        output_content.append(current_line)

        with open(output_file, 'w') as f:
            for line in output_content:
                f.write("%s\n" %line)
    print('A report has been written to {0}'.format(output_file))
    print('A total of {0} nucleotide has unknown letter (only uppercase was counted).'.format(unidentified_count))

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])