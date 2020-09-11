# -*- coding: utf-8 -*-
"""
Created on Wed Apr 08 13:17:12 2015

@author: Zewei Song
"""
def main(name_space):
    import argparse
    import textwrap
    from lib import File_IO
    import sys
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                        ------------------------
                                        By Zewei Song
                                        University of Minnesota
                                        Dept. Plant Pathology
                                        songzewei@outlook.com
                                        ------------------------'''), prog = 'fast.py -count_seqs')
    parser.add_argument("-i", "--input", help="Name of the input sequence file")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-a", "--fasta", action="store_true", help="Set file type to FASTA")
    group.add_argument("-q", "--fastq", action="store_true", help="Set file type to FASTQ")
    args = parser.parse_args(name_space)
    
    seq_file = args.input
    
    if args.fasta:
        head_symbol = '>'
        seq_type = 'fasta'
    if args.fastq:
        head_symbol = "@"
        seq_type = "fastq"
    else:
        with open(seq_file, 'rU') as f:
            header = f.read(1)
        if header == '>':
            head_symbol = ">"
            seq_type = "fasta"
            print "File type set as FASTA."
        elif header == '@':
            head_symbol = '@'
            seq_type = 'fastq'
            print "File type set as FASTQ."
        else:
            print '%s is not a valid header for FASTA or FASTQ file.' % header
            sys.exit()
    seq_content = File_IO.read_seqs(seq_file, file_type=seq_type)
    seq_count = len(seq_content)
    print "%i records found in %s." % (seq_count, seq_file)

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])