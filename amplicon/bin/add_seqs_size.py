# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 10:38:33 2015

Add a size annotation using a Qiime style OTU map.
The size annotation ";size=XXX;" is required for USEARCH to perform OTU clustering.
The sequences will also be sorted based on the size.

Please feel free to contact me for any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""

def main(Namespace):
    from lib import File_IO
    from lib import ParseOtuMap
    import argparse
    import textwrap
    import sys
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                    ------------------------
                                    By Zewei Song
                                    University of Minnesota
                                    Dept. Plant Pathology
                                    songzewei@outlook.com
                                    ------------------------'''), prog='fast.py -add_seqs_size')
    parser.add_argument('-i', '--input', help='Input FASTA file')
    parser.add_argument('-map', help='Input OTU map file')
    parser.add_argument('-o', '--output', help='Output FASTA file')
    args = parser.parse_args(Namespace)
    
    input_fasta = args.input
    input_map = args.map
    output_fasta = args.output
    
    print 'Reading in OTU map ...'
    otu_map = ParseOtuMap.read_otu_map(input_map)
    print 'Reaidng in sequence file ...'
    fasta = File_IO.read_seqs(input_fasta)
    print 'Found %i OTUs in the map file, found %i sequences in the sequence file.' % (len(otu_map), len(fasta))
    
    count = 0
    print 'Adding size annotation ...'
    for record in fasta:
        try:
            size = len(otu_map[record[0]])
            if record[0][-1] == ';':
                record[0] += ('size=%i;' % size)
            else:
                record[0] += (';size=%i;' % size)
            record.append(size)
            count += 1
            print 'Annotating %i sequence ...' % count + '\b' * 100,
        except KeyError:
            print "Can not find %s in the OTU map file." % record[0]
            sys.exit()
    print
    print 'Sorting the annotated sequences ...'
    fasta.sort(key=lambda x: x[-1], reverse=True)
    
    print 'Writing to a new FASTA file ...'
    with open(output_fasta, 'w') as f:
        for record in fasta:
            f.write('>%s\n' % record[0])
            f.write('%s\n' % record[1])
    print 'Sequences with size annotations saved in %s.' % output_fasta

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])