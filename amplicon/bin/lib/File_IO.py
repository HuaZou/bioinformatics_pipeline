# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 11:29:35 2015
Some simple manipulation on files and folders
@author: Zewei Song
"""
from __future__ import print_function
from __future__ import division
import gzip

# Create a new folder
def mk_dir(folder):
    # Create a new folder, will check its availability first.
    if check_dir(folder):
        from os import makedirs

        makedirs(folder)


# Check is a folder already exist. exist=Flase
def check_dir(folder):
    # Check if a folder already exists, if not return the value True.
    import os
    import sys

    if os.path.exists(folder):
        if not os.listdir(folder):
            return False
        else:
            print('Folder: <%s> already exists. Program aborted.' % folder)
            sys.exit()
    else:
        return True


# Return a list of file names in a folder
def file_list(folder):
    # Get the file list in a folder, but not any sub-folders.
    from os import walk

    f = []
    for (dirpath, dirnames, filenames) in walk(folder):
        f.extend(filenames)
        break
    return f


# Check the availability of a file, no file = True, exist=False
def check_file(input_filename):
    # Check if a file exists, if not return the value True.
    from os import path

    if path.isfile(input_filename):
        return False
    else:
        return True


# Read the entire content of the file by lines and store in list
def read_file(input_filename):
    with open(input_filename, 'rU') as f:
        content = f.readlines()  # f.read().split('\n') take more memory than f.readlines()
    return content  # Return a list with all lines in the file. Lines contain \n at the end


# Read in a sequence file (FASTA or FASTQ) and store records in a list
def read_seqs(input_filename, file_type='default', output='default'):
    import sys
    # Check for record header
    if input_filename[-3:] == '.gz': # Check if the file extension is .gz, beware that this is a loose condition for gz file.
        with gzip.open(input_filename, 'r') as f:
            content = [i.decode(encoding='utf-8') for i in f.readlines()]
            head_symbol = content[0][0]
    else:
        with open(input_filename, 'rU') as f:
            content = f.readlines()
            head_symbol = content[0][0]

    if file_type == 'default':
        # Try to guess file format if not provided
        if head_symbol == '>':
            file_type = 'fasta'
        elif head_symbol == '@':
            file_type = 'fastq'
        else:
            print('%s is not a correct header for FASTA or FASTQ, please check you file.' % head_symbol)
            sys.exit()
    if file_type in ['fasta', 'Fasta', 'FASTA']:
        if content[0][0] == '@':
            print('%s seems to be a FASTQ file. Please use the correct format.' % input_filename)
            sys.exit()
        line_num = 2
        seq_line_num = 2
    elif file_type in ['fastq', 'Fastq', 'FASTQ']:
        if content[0][0] == '>':
            print('%s seems to be a FASTA file. Please use the correct format.' % input_filename)
            sys.exit()
        line_num = 4
        seq_line_num = 4
    else:
        print('Please specify a right format [fasta | fastq].')
        sys.exit()

    # set line number to 2 if want to export FASTA for a FASTQ file.
    if output in ['fasta', 'Fasta', 'FASTA']:
        seq_line_num = 2

    # Store record in list
    seq_num = int(len(content) / line_num)  # Total number of records
    output_content = []
    for i in range(seq_num):
        temp = []  # Store current record
        for line in range(seq_line_num):
            # Loop through each line in the current record (2 for fasta, 4 for fastq)
            temp.append(content[i * line_num + line][:-1])  # [:-1] remove '\n' at the end of each line
            content[i * line_num + line] = ''  # Remove processed content make it faster
        temp[0] = temp[0][1:]  # Remove header from sequence label
        output_content.append(temp)
    return output_content  # Return a list with record in the file


# Read in a multiple line fasta file, much slower than read_seqs().
def read_fasta_multiline(input_content, head_symbol='>'):
    corrected_content = []
    for line in input_content:
        if line[0] == head_symbol:
            corrected_content.append([line.strip('\n')[1:], ''])
        else:
            corrected_content[-1][1] += line.strip('\n')
    return corrected_content


# Write a sequence list to a new file.
def write_seqs(input_content, output_file, checker=True, overwrite=False):
    global head_symbol, head_symbol
    import sys
    seqs_line_num = len(input_content[0])
    if seqs_line_num == 2:
        head_symbol = '>'
        line_num = 2
    elif seqs_line_num == 4:
        head_symbol = '@'
        line_num = 4

    if checker:
        if not check_file(output_file):  # File already exist
            print('%s already exist, please check your working folder to avoid any potential loss.' % output_file)
            sys.exit()

    if not overwrite:  # Append to current file
        with open(output_file, 'a') as f:
            for record in input_content:
                record[0] = head_symbol + record[0]
                for i in range(line_num):
                    f.write('%s\n' % record[i])
    elif overwrite:  # Overwrite current file
        with open(output_file, 'w') as f:
            for record in input_content:
                record[0] = head_symbol + record[0]
                for i in range(line_num):
                    f.write('%s\n' % record[i])
    return len(input_content)  # Return total sequence number


# Name a file using prefix and suffix
def name_file(input_filename, prefix, suffix):
    # Add prefix and "_" before each file
    # Replace the last .extension with .suffix
    # Add .suffix is no extension in the file
    dot_position = [i for i in range(len(input_filename)) if input_filename.startswith('.', i)]
    if not dot_position:  # input filename does not has "."
        if prefix == "":
            input_filename = input_filename + "." + suffix
        else:
            input_filename = prefix + "_" + input_filename + "." + suffix
    else:
        if prefix == "":
            input_filename = input_filename[:dot_position[-1]] + "." + suffix + input_filename[dot_position[-1]:]
        else:
            input_filename = prefix + "_" + input_filename[:dot_position[-1]] + "." + suffix + input_filename[dot_position[-1]:]
    return input_filename