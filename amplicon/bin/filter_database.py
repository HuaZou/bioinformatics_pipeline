# -*- coding: utf-8 -*-
"""
Created on Tue Dec 08 10:51:54 2015

This is a trial script for filtering the UNITE database.
All records with 'unidentified' in its name are discarded, resutling a clean reference database for taxonomic assignment.

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
                                        ------------------------'''), prog = '-filter_database')

    parser.add_argument("-i", "--input", help="Name of the input FASTA file.")
    parser.add_argument("-o", "--output", help="Name of the output FASTA file")
    args = parser.parse_args(name_space)

    database = File_IO.read_seqs(args.input)
    count = len(database)
    print "Reading in %s ..." % args.input
    print "%s contains %i records." % (args.input, count)
    
    count_filter = 0
    database_cleaned = []
    for record in database:
        if record[0].find('unidentified') == -1: # check if current record contain 'unidentified' taxonomic level.
            database_cleaned.append(record)
            count_filter += 1
    
    print "%i records contain 'unidentified' string." % (count - count_filter)
    
    count_write = File_IO.write_seqs(database_cleaned, args.output)
    print "Filtered database is saved in %s with %i records." % (args.output, count_write)

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])