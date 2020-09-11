# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 15:57:35 2015

Please feel free to contact with any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""


def check_ambiguous(seq, length):
    if seq.count('N') > length:
        return True
    else:
        return False
    

def check_homop(seq, length):
    seq_homop = [i*length for i in ['A','T','C','G']]
    for i in seq_homop:
        if i in seq:
            return True
    return False
    

def count_bases(seq):
    bases_number = {'A':0,'T':0,'C':0,'G':0}
    for base in bases_number:
        bases_number[base] = seq.count(base)
    return bases_number

        
def count_ambiguous(seq):
    return seq.count('N')


def count_homop(seq):
    count_homop = {'All':0,'A':0,'T':0,'C':0,'G':0}
    base_list = [i for i in count_homop]
    seq_length = len(seq)
    for base in base_list:
        for i in range(seq_length):
            homop = base*(i+1) # Generate a homopolymer with certain length
            if homop in seq:
                count_homop[base] = i+1
            else:
            # Do not found homopolyer at the given length, break the loop
                break
    max_length = max([count_homop[i] for i in count_homop])
    count_homop['All'] = max_length
    return count_homop


def pick_seqs(fasta,name_list):
    # Pick sequences based on the name list input:
    # This method is slow compared to using a dictionary
    fasta_picked = []
    for record in fasta:
        if record[0] in name_list:
            fasta_picked.append(record)
            del name_list[name_list.index(record[0])]
        if name_list == []:
            break
    return fasta_picked


def make_dict(seqs):
    # Convert a seqs list into dictionary using sequence name as keys
    seqs_dict = {}
    for record in seqs:
        seqs_dict[record[0]] = record[1:]
    return seqs_dict
    
def nucl_freq(input_seq, tail = False):
# Count the frequency of nucleotide. Return a dictionary.
    # Reverse the input sequences if tail = True
    if tail == True:
        for index, record in enumerate(input_seq):
            input_seq[index][1] = record[1][::-1]
    else: 
        pass

    # Get the maximum sequence length
    seq_len = []
    for record in input_seq:
        seq_len.append(len(record[1]))
    max_seq_len = max(seq_len)
    
    # Create dictionary for nucleotide counting    
    nucl_dict = {}
    nucl_list = ['A','T','C','G','N']
    for i in range(max_seq_len):
        nucl_dict[i] = {}
        for nucl in nucl_list:
            nucl_dict[i][nucl] = 0
    
    # Count the frequency of nucleotide  
    unidentified_count = 0
    for record in input_seq:
        for index, nucl in enumerate(record[1]):
            try:
                nucl_dict[index][nucl] += 1
            except KeyError:
                unidentified_count += 1
                #print "Unidentified nucleotide %s in %s" % (nucl, record[0])
    return nucl_dict, unidentified_count