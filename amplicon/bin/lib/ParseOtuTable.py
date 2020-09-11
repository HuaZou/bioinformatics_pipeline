# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Create on 4/21/2015.

This library will parse a tab-delimited OTU table into dictionaries.

Please feel free to contact me with any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""
from __future__ import print_function
from __future__ import division

class parser_otu_table(object):
    def __init__(self, file_path, meta_col='taxonomy'):
        with open(file_path, 'r') as f:
            temp = f.readlines()
        table = []
        for line in temp:
            table.append(line.strip('\n').split('\t'))
        size = len(table[0])
        for line in table:
            if len(line) < size:
                line += [''] * (size - len(line))

        # Get id for sample, species, amd meta data
        try:
            meta_position = table[0].index(meta_col)  # Begining position of meta data
        except ValueError:
            meta_position = len(table[0]) + 1
        self.sample_id = table[0][1:meta_position]  # Second column till the first meta column
        self.meta_id = table[0][meta_position:]  # Start from the first meta column to the end
        self.species_id = [i[0] for i in table[1:]]  # First column starting from the second row
        self.header = ['OTU_ID'] + self.sample_id + self.meta_id

        # Convert all abundance to intger
        for line in table[1:]:
            try:
                line[1:meta_position] = map(int, line[1:meta_position])
            except ValueError:
                print("There are non-number value in your OTU table.")
                import sys
                sys.exit()

        # Get sample, species, and meta data matrix with head names
        self.species_matrix = [i[:meta_position] for i in table[1:]]
        temp = [i[1:meta_position] for i in table]
        self.sample_matrix = [list(i) for i in zip(*temp)]
        temp = [i[meta_position:] for i in table]
        self.meta_matrix = [list(i) for i in zip(*temp)]


    # Generate a dictionary using sample name, OTU name
    def sample_dict(self):
        sample = {}
        for s in self.sample_id:
            sample[s] = {}
        for line in self.sample_matrix:
            abundance = line[1:]
            for i in range(len(self.species_id)):
                if int(abundance[i]) > 0:
                    sample[line[0]][self.species_id[i]] = int(abundance[i])
        return sample

    # Generate a dictionary using OTU name, sample name
    def species_dict(self):
        species = {}
        for s in self.species_id:
            species[s] = {}
        for line in self.species_matrix:
            for i in range(len(self.sample_id)):
                if int(line[1:][i]) > 0:
                    species[line[0]][self.sample_id[i]] = int(line[1:][i])
        return species


    # Generate a dictionary using mata name, OTU name
    def meta_dict(self):
        meta = {}
        for m in self.meta_id:
            meta[m] = {}
        for line in self.meta_matrix:
            for i in range(len(self.species_id)):
                meta[line[0]][self.species_id[i]] = line[1:][i]
        return meta


# Generate a new tab delimited OTU table sorted by sum of OTU abundance
def sorted_table(otu_table_parser, new_file_path='otu_table_sorted.txt', remove_zero = False):
    otu_id = otu_table_parser.species_id
    otu_dict = otu_table_parser.species_dict()
    meta_id = otu_table_parser.meta_id
    meta_dict = otu_table_parser.meta_dict()

    otu_abundance_list = []
    for otu in otu_id:
        sum_abundance = sum(otu_dict[otu].values())
        otu_abundance_list.append([sum_abundance, otu])
    otu_abundance_list = sorted(otu_abundance_list, reverse = True)

    if remove_zero:
        new_otu_list = []
        for item in otu_abundance_list:
            if item[0] > 0:
                new_otu_list.append(item)
            else:
                pass
        otu_abundance_list = new_otu_list

    output_content = [otu_table_parser.header]
    for item in otu_abundance_list:
        otu_name = item[1]
        current_line = [otu_name]
        for sample in otu_table_parser.sample_id:
            try:
                current_line.append(otu_dict[otu_name][sample])
            except KeyError:
                current_line.append(0)

        if len(meta_id) >0:
            for meta_column in meta_id:
                current_meta = meta_dict[meta_column][otu_name]
                current_line.append(current_meta)
        else:
            pass

        output_content.append(current_line)
    write_content(output_content, new_file_path)

    return


# Output a tab delimited OTU table using a sample dict and meta_dict
def write_sample_dict(sample_dict, meta_dict, otu_id, output_file_path):
    sample_id = list(sample_dict.keys())
    sample_id.sort()
    meta_id = list(meta_dict.keys())

    header = ['OTU_ID'] + sample_id + meta_id
    content = [header]
    for otu in otu_id:
        current_line = [otu]
        for sample in sample_id:
            try:
                current_abundance = sample_dict[sample][otu]
                current_line.append(str(current_abundance))
            except KeyError:
                current_line.append('0')
        if len(meta_id) > 0:
            for meta_column in meta_id:
                current_meta = meta_dict[meta_column][otu]
                current_line.append(current_meta)
        else:
            pass
        content.append(current_line)

    write_content(content, output_file_path)
    return

# Output a tab delimited OTU table using a sample dict and meta_dict, and re-organized the sample order using a list
# Need to add a checking for sample name later
def write_sample_dict_newlist(sample_dict, meta_dict, otu_id, output_file_path, new_order):
    #sample_id = sample_dict.keys()
    #sample_id.sort()
    meta_id = meta_dict.keys()

    header = ['OTU_ID'] + list(new_order) + list(meta_id)
    content = [header]
    for otu in otu_id:
        current_line = [otu]
        for sample in new_order:
            try:
                current_abundance = sample_dict[sample][otu]
                current_line.append(str(current_abundance))
            except KeyError:
                current_line.append('0')
        if len(meta_id) > 0:
            for meta_column in meta_id:
                current_meta = meta_dict[meta_column][otu]
                current_line.append(current_meta)
        else:
            pass
        content.append(current_line)

    write_content(content, output_file_path)
    return

def write_content(input_content, output_file_path):
    with open(output_file_path, 'w') as f:
        for line in input_content:
            line = [str(i) for i in line]
            line = '\t'.join(line)
            f.write('%s\n' % line)