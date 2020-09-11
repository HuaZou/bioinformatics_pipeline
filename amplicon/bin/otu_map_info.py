# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Create on .

Please feel free to contact me with any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""
def main(name_space):
    import argparse
    import textwrap
    from lib import ParseOtuMap
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                        ------------------------
                                        By Zewei Song
                                        University of Minnesota
                                        Dept. Plant Pathology
                                        songzewei@outlook.com
                                        ------------------------'''), prog = 'fast.py -otu_map_info')
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-qiime_map", help="Input Qiime OTU map")
    group.add_argument('-fast_map', help='Input FAT map.')
    args = parser.parse_args(name_space)
    
    if args.qiime_map != None:
        input_map = args.qiime_map
        method = 'qiime'
    elif args.fast_map != None:
        input_map = args.fast_map
        method = 'fast'
    
    print 'Reading in %s ...' % input_map
    otu_map = ParseOtuMap.read_otu_map(input_map)
    
    if method == 'qiime':
        otu_map_parser = ParseOtuMap.otu_map_parser(otu_map)
        otu_num = otu_map_parser.derep_count
        seq_num = otu_map_parser.seqs_count
        otu_max = otu_map_parser.max_derep
        otu_min = otu_map_parser.min_derep
        otu_ave = otu_map_parser.ave_derep
        
        print 'The map type is: Qiime.'        
        print 'OTU = %i; Sequence = %i' % (otu_num, seq_num)
        print 'Max abundance = %i; Min abundance = %i; Ave abundance = %i' % (otu_max, otu_min, otu_ave)

    elif method == 'fast':
        otu_map = ParseOtuMap.read_fast_output(input_map)
        otu_map_parser = ParseOtuMap.fast_output_parser(otu_map)
        
        seq_list = otu_map_parser.get_seqs()
        otu_num = len(seq_list)
        otu_type = otu_map_parser.fast_type
        seq_num = sum([i[2] for i in seq_list])
        otu_max = max([i[2] for i in seq_list])
        otu_min = min([i[2] for i in seq_list])
        otu_ave = seq_num / otu_num
        
        print 'The map type is: FAST-{0}.'.format(otu_type)
        print 'OTU = %i; Sequence = %i' % (otu_num, seq_num)
        print 'Max abundance = %i; Min abundance = %i; Ave abundance = %i' % (otu_max, otu_min, otu_ave)

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])