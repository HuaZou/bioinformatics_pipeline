#!usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create on 4/24/2015.

Please feel free to contact me for any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""

from lib import ParseOtuTable
import argparse
import textwrap

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=textwrap.dedent('''\
                                    ------------------------
                                    By Zewei Song
                                    University of Minnesota
                                    Dept. Plant Pathology
                                    songzewei@outlook.com
                                    ------------------------'''))
parser.add_argument("-otu", help="Input OTU folder")
parser.add_argument("-mock_list", help="A list of species name in the mock community, Genus and species should be connected with _")
parser.add_argument("-mock_column", default='mock', help="Name of the mock community sample")
parser.add_argument("-o", "--output", help="Output OTU table")
parser.add_argument('-min_length', default = 0.9, help='Minimum mathced length')
parser.add_argument('-min_pident', default = 90, help='Minimum Pident')
args = parser.parse_args()

input_otu = args.otu
input_mock = args.mock_list
mock_column = args.mock_column
output_report = args.output
min_length = float(args.min_length)
min_pident = float(args.min_pident)

#input_otu = 'otu_table.tax.txt'
#input_mock = 'mock_list.txt'
#mock_column = 'Mock-community'

otu_table = ParseOtuTable.parser_otu_table(input_otu)
otu_sample = otu_table.sample_dict()
mock_sample = otu_sample[mock_column]
otu_meta = otu_table.meta_dict()

#%%
mock_list = []
mock_hit = {}
with open(input_mock, 'rU') as f:
    for line in f:
        mock_list.append(line.strip('\n'))
        mock_hit[line.strip('\n')] = {}
mock_hit['other'] = {}


for otu_name in mock_sample:
    if otu_meta['taxonomy'][otu_name] != 'no_blast_hit':
        coverage = int(otu_meta['Subject_Len'][otu_name])/float(int(otu_meta['Query_Len'][otu_name]))
        if coverage >= min_length and float(otu_meta['Pident'][otu_name]) >= min_pident:
            match_checker = False
            for mock_species in mock_list:
                if otu_meta['taxonomy'][otu_name].find(mock_species) != -1:  # Found mock species in current OTU
                    mock_hit[mock_species][otu_name] = [mock_sample[otu_name], otu_meta['taxonomy'][otu_name], str(coverage), otu_meta['Pident'][otu_name]]
                    match_checker = True
                    break
            if not match_checker:
                mock_hit['other'][otu_name] = [mock_sample[otu_name], otu_meta['taxonomy'][otu_name], str(coverage), otu_meta['Pident'][otu_name]]
        else:
            mock_hit['other'][otu_name] = [mock_sample[otu_name], otu_meta['taxonomy'][otu_name], str(coverage), otu_meta['Pident'][otu_name]]
    else:
        mock_hit['other'][otu_name] = [mock_sample[otu_name], otu_meta['taxonomy'][otu_name], '', otu_meta['Pident'][otu_name]]

#%%
# Print report

mock_list += ['other']
#%%
report_general = [['Mock species','Total sequence','Total OTUs']]
for species in mock_list:
    hit = mock_hit[species]
    seqs_sum = []
    otu_count = len(hit)    
    for k, v in hit.items():
        seqs_sum.append(v[0])
    seqs_sum = sum(seqs_sum)
    report_general.append([species,seqs_sum,otu_count])
#%%
report_detail = []
count = 1
for species in mock_list:
    report_detail.append([str(count) + '. ' + species + ':', 'OTU', 'Abundance', 'Matched length', 'Pident', 'Taxnonomy'])
    for key, value in mock_hit[species].items():
        report_detail.append(['', key, str(value[0]), value[2], value[3], value[1]])
    count += 1
    report_detail.append([''])
#%%
with open(output_report, 'w') as f:
    f.write('%s\n' % 'I. General report on mock community')
    for line in report_general:
        line = [str(i) for i in line]
        f.write('%s\n' % '\t'.join(line))
    f.write('\n%s\n' % 'II. Detail report')
    for line in report_detail:
        f.write('%s\n' % '\t'.join(line))