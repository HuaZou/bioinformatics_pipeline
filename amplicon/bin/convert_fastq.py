# -*- coding: utf-8 -*-
"""
Created on Mon Apr 06 16:52:45 2015
Last modified on 4/8/2015
Convert a fastq file to fasta file. Current speed is about 50,000/s or 100,000,000 per 20s for a 4gb FASTQ file.
Smaller file should be faster.
FASTQ file need to be in the write format and identical length in sequence and quality scores.
Right now it will load the entire file into the memory.

python convert_fastq.py -h for help document

Feel free to contact me for any question.
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""
def main(name_space):
    import argparse
    import textwrap
    from lib import File_IO
    import time
    import sys
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                        ------------------------
                                        By Zewei Song
                                        University of Minnesota
                                        Dept. Plant Pathology
                                        songzewei@outlook.com
                                        ------------------------'''))
    parser.add_argument("-i", "--input", help="Convert a FASTQ file.")
    parser.add_argument("-o", "--output", help="Name of the output FASTA file")
    #parser.add_argument("-q", "--qual", action="store_true", help="Output Qual file")
    args = parser.parse_args(name_space)
    
    fasta_file = args.output
    #qual = args.qual
    
    if args.input:
        fastq_file = args.input
        start = time.time()
        print("Loading %s ..." % fastq_file)
        fasta_content = File_IO.read_seqs(fastq_file, file_type='fastq', output='fasta')
        print('Converting to FASTA ...')
        record_num = File_IO.write_seqs(fasta_content, fasta_file, checker=False, overwrite=True)
        print("Converted %d records in %s ..." % (record_num, fastq_file))
    
        end = time.time()
        used_time = round(end - start, 2)
        print("It took %s sec to convert (%s seqs/s).\nFASTA file saved in %s." % (
            str(used_time), str(round(record_num/used_time,0)), fasta_file))
    #    if qual:
    #        print "Quality scores saved in %s." % (File_IO.name_file(fasta_file, '', 'qual'))
    
    else:
        print("Please specify a FASTQ file.")
        sys.exit()

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])