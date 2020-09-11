# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 22:02:28 2015

Please feel free to contact me with any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""

import argparse, textwrap
import shlex, subprocess

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=textwrap.dedent('''\
                                    ------------------------
                                    By Zewei Song
                                    University of Minnesota
                                    Dept. Plant Pathology
                                    songzewei@outlook.com
                                    ------------------------'''))
parser.add_argument("-otu", help="Input OTU folder")
parser.add_argument("-db", help="Path to BLAST database")
parser.add_argument("-query", help="Path to the representative sequences")
parser.add_argument("-o", "--output", help="Output report")
args = parser.parse_args()

input_otu = args.otu
input_db = args.db
input_query = args.query
output_report = args.output

#command_line = ['blastn','-db','unite_02.03.2015','-query','Corn_Pytobiome_rep_seq.fa','-max_target_seqs','1','-outfmt','"6 qseqid stitle qlen length pident evalue"']
#command_line = ['blastn','-db', input_db,'-query', input_query,'-max_target_seqs','1','-outfmt','"6 qseqid stitle qlen length pident evalue"']

#cmd = 'blastn -db unite_02.03.2015 -query Corn_Pytobiome_rep_seq.fa -max_target_seqs 1 -outfmt "6 qseqid stitle qlen length pident evalue"'
#%%
#import shlex, subprocess
input_otu = 'otu_table.txt'
input_db = 'blast/unite_02.03.2015'
input_query = 'cluster_usearch.fa'
output_report = 'blast_report.txt'
cmd = 'blastn -db ' + input_db + ' -query ' + input_query + ' -max_target_seqs 1 -outfmt "6 qseqid stitle qlen length pident evalue"'

args = shlex.split(cmd)
p = subprocess.Popen(args, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    stdin=subprocess.PIPE)
out, err = p.communicate()
#%%
from lib import ParseOtuTable
otu = ParseOtuTable.parser_otu_table(input_otu)
otu_species = otu.species_matrix
#%%
otu_info = {}
for item in otu_species:
    otu_info[item[0]] = {'abundance':sum(item[1:]),'match_len':0}
#%%
temp = out.split('\n')
blast_result = []
for line in temp:
    blast_result.append(line.split('\t'))
del blast_result[-1]
#%%
for record in blast_result:
    match_len = float(record[3])/int(record[2])
    otu_info[record[0]]['match_len'] = match_len
#%%
threshold = range(1,101)
threshold = [i/100.0 for i in threshold]
#threshold = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
report = {}
for t in threshold:
    report[str(t)] = [0,0]

for t in threshold:
    for key, value in otu_info.items():
        if value['match_len'] >= t:
            report[str(t)][0] += value['abundance']
            report[str(t)][1] += 1
#%%
r = report.items()
r.sort()
#%%
with open('blast_report.txt', 'w') as f:
    for record in r:
        line = record[0] + '\t' + str(record[1][0]) + '\t' + str(record[1][1])
        f.write('%s\n' % line)
#%%
    