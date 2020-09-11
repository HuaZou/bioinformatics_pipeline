# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 13:43:00 2015

This library will parse the taxonomy in an OTU table (tab delimited) into a given taxonomic level.

Please feel free to contact me with any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""

# Parse a taxonomic string, break it down into taxonomic levels. Now only support the UNITE format.
# The broken down string is save as a dictionary
def parse_taxonomy(tax_string):
    parser = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    if tax_string == 'no_blast_hit' or tax_string == "no_hit":
        string = ['k__Fungi']
        for level in parser[1:]:
            string.append(level + 'unidentified')
    else:
        string = tax_string[tax_string.index('k__'):]
        string = string.split(';')
    
    string_dict = {}
    for index, parser in enumerate(parser):
        string[index] = string[index].strip(parser)
        string_dict[parser[0]] = string[index]
    return string_dict


# Parse the taxonomy string of all OTUs
def taxonomy_info(input_otu_table, tax_level):
    from lib import ParseOtuTable
    otu_table = ParseOtuTable.parser_otu_table(input_otu_table)
    tax_dict = otu_table.meta_dict()['taxonomy']
    otu_dict = otu_table.species_dict()
    sample_id = otu_table.sample_id
    
    abundance_dict = {}
    richness_dict = {}
    
    for otu in otu_dict.keys():
        current_tax = parse_taxonomy(tax_dict[otu])
        try:
            abundance_dict[current_tax[tax_level]]
            for sample in otu_dict[otu].keys():
                abundance_dict[current_tax[tax_level]][sample] += otu_dict[otu][sample]
                richness_dict[current_tax[tax_level]][sample] += 1
        except KeyError:
            abundance_dict[current_tax[tax_level]] = {}
            richness_dict[current_tax[tax_level]] = {}
            for sample in otu_table.sample_id:
                abundance_dict[current_tax[tax_level]][sample] = 0
                richness_dict[current_tax[tax_level]][sample] = 0
            for sample in otu_dict[otu].keys():
                abundance_dict[current_tax[tax_level]][sample] = otu_dict[otu][sample]
                richness_dict[current_tax[tax_level]][sample] += 1
    
    report_dict = {'abundance': abundance_dict, 'richness': richness_dict, 'taxonomy_level': tax_level, 'sample_id': sample_id}
    
    return report_dict

# Change the output of taxonomy_info() into writable form (a list):
def taxonomy_output(input_report_dict):
    header = ['tax_level'] + input_report_dict['sample_id']    
    
    output_dict = {'abundance':[header], 'richness':[header]}    
    for tax_level in input_report_dict['abundance'].keys():
        current_abundance = [tax_level]
        current_richness = [tax_level]
        for sample in header[1:]:
            current_abundance.append(str(input_report_dict['abundance'][tax_level][sample]))
            current_richness.append(str(input_report_dict['richness'][tax_level][sample]))
        output_dict['abundance'].append(current_abundance)
        output_dict['richness'].append(current_richness)
    return output_dict
        