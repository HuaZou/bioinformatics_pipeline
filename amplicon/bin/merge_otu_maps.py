# -*- coding: utf-8 -*-
"""
Create on April 20th, 2015.

Merge two Qiime style OTU maps if the smaller one is the clusters of the larger one.

The small map can be a subset of the OTU names of the larger map.

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
    import sys
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                        ------------------------
                                        By Zewei Song
                                        University of Minnesota
                                        Dept. Plant Pathology
                                        songzewei@outlook.com
                                        ------------------------'''), prog = '-merge_otu_maps')
    parser.add_argument('-map_large',help='The large OTU map')
    parser.add_argument('-map_small',help='The small OTU map')
    parser.add_argument('-o','--output', help='Output of the merged map')
    args = parser.parse_args(name_space)
    
    input_mapL = args.map_large
    input_mapS = args.map_small
    output_file = args.output
    
    mapL = ParseOtuMap.read_otu_map(input_mapL)
    mapS = ParseOtuMap.read_otu_map(input_mapS)
    mapL_parser = ParseOtuMap.otu_map_parser(mapL)
    mapS_parser = ParseOtuMap.otu_map_parser(mapS)
    print 'Large map: %i OTUs, %i sequences' % (mapL_parser.derep_count, mapL_parser.seqs_count)
    print 'Small map: %i OTUs, %i sequences' % (mapS_parser.derep_count, mapS_parser.seqs_count)
    
    map_merged = {}
    for key_s, value_s in mapS.items():
        map_merged[key_s] = []
        for name in value_s:
            try:
                map_merged[key_s] += mapL[name]
            except KeyError:
                print 'Cannot find %s in %s. Please check you OTU map files.' % (name, input_mapL)
                sys.exit()
    
    map_merged_parser = ParseOtuMap.otu_map_parser(map_merged)
    print 'Merged OTU map:'
    print 'OTUs=%i, Sequences=%i, Max=%i, Min=%i, Ave=%i' %(map_merged_parser.derep_count, map_merged_parser.seqs_count,\
                                                            map_merged_parser.max_derep, map_merged_parser.min_derep,\
                                                            map_merged_parser.ave_derep)
    ParseOtuMap.write_otu_map(map_merged, output_file)
    print 'Saved to %s.' % output_file

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])