# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 21:59:21 2015

Please feel free to contact me with any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""
def main(name_space):
    import argparse
    from lib import File_IO
    import textwrap
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                        ------------------------
                                        By Zewei Song
                                        University of Minnesota
                                        Dept. Plant Pathology
                                        songzewei@outlook.com
                                        ------------------------'''), prog = 'correct_fasta')
    parser.add_argument("-i", "--input", help="FASTA file need to be fixed")
    parser.add_argument("-o", "--output", help="Name of the output FASTA file")
    parser.add_argument('-head', default='>',help='Specify a head symbol if not >')
    args = parser.parse_args(name_space)
    
    file_origin = args.input
    if args.output:
        file_corrected = args.output
    else:
        file_corrected = 'corrected_'+file_origin
    head = args.head
    
    fasta_corrected = File_IO.read_fasta_multiline(File_IO.read_file(file_origin),head_symbol=head)
    count = File_IO.write_seqs(fasta_corrected,file_corrected,checker=False,overwrite=True)
    
    print 'Checked %d sequences in %s and saved in %s.' % (count, file_origin, file_corrected)

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])