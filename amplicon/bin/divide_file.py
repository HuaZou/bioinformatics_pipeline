# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 16:48:19 2015

Divide a large fastq file into smaller files.

Please feel free to contact me with any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""

#from lib import File_IO
import argparse
import sys
#import time

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Name of the input file, can be FASTA or FASTQ')
parser.add_argument('-n', '--number', help='Number of new files')
args = parser.parse_args()

input_file = args.input
number_file = int(args.number)

# Count the total line number of the FASTQ file
line_total = 0
with open(input_file, 'rU') as f:
    for line in f:
        line_total += 1
        sys.stdout.write('Total line number: %i \r' %line_total)
print
if line_total%4 != 0:
    print 'Number of line is not in the set of 4.'
    sys.exit()

#%%
#line_total = 1024
#number_file = 5

seq_total = line_total/4
# Assign lines to each smaller file
subset_seq = [0] * number_file
base_seq = seq_total/number_file
extra_seq = seq_total%number_file
for i in range(number_file):    
    subset_seq[i] = base_seq
    if extra_seq >0 :
        subset_seq[i] += 1
        extra_seq -= 1

# Convert seq number to line number
subset_line = [i*4 for i in subset_seq]
#%%
# Generate the new file list
#input_file = 'test.fastq'
#number_file = 5
name_files = [''] * number_file
n = 0
for i in range(number_file):
    name_files[i] = input_file.replace('.fastq', '_' + str(n) + '.fastq')
    n += 1
#%% Write sequences to the smaller files, Need to work on this
    # Need to save in temp and write at once.
current_file = 0
temp_content = []
with open(input_file, 'rU') as f:
    for line in f:
        # Set for current file
        if subset_line[current_file] > 0:
            subset_line[current_file] -= 1
            temp_content.append(line)
        else:
            output_file = open(name_files[current_file], 'w')
            print ('Writing %s ...' % name_files[current_file])
            for line in temp_content:                
                output_file.write('%s' % line)
            temp_content = []
            current_file += 1

# Write the last sub-file
with open(name_files[-1], 'w') as f:
    print ('Writing %s ...' % name_files[-1])
    for line in temp_content:
        f.write('%s' % line)