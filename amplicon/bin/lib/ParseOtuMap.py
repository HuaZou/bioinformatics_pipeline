# -*- coding: utf-8 -*-
"""
Created on Thu Apr 02 22:01:14 2015

Please feel free to contact me for any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""
from __future__ import print_function
from __future__ import division
#%%
def read_otu_map(filename):
    from lib import File_IO
    OtuMap = File_IO.read_file(filename)
    MapDict = {}
    for line in OtuMap:
        names = line.strip('\n').split('\t')
        MapDict[names[0]]=names[1:]
    return MapDict
#%%
class otu_map_parser(object):
    # Get basic information of an OTU map
    def __init__(self, Map):
        self.derep_count = len(Map)
        derep_size = [len(Map[i]) for i in Map]
        self.seqs_count = sum(derep_size)
        self.max_derep = max(derep_size)
        self.min_derep = min(derep_size)
        self.ave_derep = self.seqs_count/float(self.derep_count)


def get_rep_list(Map):
    rep_list = [i for i in Map]
    return rep_list

def filter_by_name(MapDict,otu_list,method='discard'):
# Pick or discard OTUs in the otu list from a Qiime OTU map.
    if method == 'discard':
        for otu in otu_list:
            del MapDict[otu]
        return MapDict
    elif method == 'keep':
        MapDict_keep = {}
        for otu in otu_list:
            MapDict_keep[otu] = MapDict[otu]
        return MapDict_keep
    else:
        print('Wrong method: %s is not a valid method.'% method)
        print('use discard or keep')

def filter_by_size(MapDict,min_size=2):
#Filter OTU map by minimum size.
#By default it will remove all singletons.
    MapDict_Size = MapDict.copy()
    for item in MapDict:
        cluster_size = len(MapDict[item])
        if cluster_size < min_size:
            del MapDict_Size[item]
        else:
            pass
    return MapDict_Size


#%%
def write_otu_map(MapDict,output_file='new_map.txt'):
#Not sure for now if I should sort the otu names by number
#Probably not, it can be converted to otu table in Qiime.
    with open(output_file, 'w') as f:
        for key, value in MapDict.items():
            cluster_seqs = '\t'.join(value)
            f.write('%s\t%s\n'%(key,cluster_seqs))


def extract_all_seqs(Map, otu_list):
# Extract all seqeunce names associated with the OTUs in the OTU list provided.
    import sys
    extracted = []
    for name in otu_list:
        try:
            for seq_name in Map[name]:
                extracted.append(seq_name)
        except ValueError:
            print('Cannot find %s in the provided Qiime map.')
            sys.exit()
    return extracted



#%%
def generate_fast_output(input_otu_map, input_centroid, real_sample = False, separator = ';'):
# Create a FAST style file that contains information of sequences and sample counts.
    fast_dict = {}
    centroid_dict = {}
    centroid_list = []
    
    # Remove USEARCH size label from sequence labels
    for record in input_centroid:
        if record[0].find(separator) != -1:
            current_label = record[0][:record[0].find(separator)]
        else:
            current_label = record[0]
        centroid_dict[current_label] = record[1]
        centroid_list.append(current_label)

    new_otu_map = {}
    for key, value in input_otu_map.items():
        if key.find(separator) != -1:
            new_key = key[:key.find(separator)]
        else:
            new_key = key

        new_sample = []
        for item in value:
            if real_sample:
                new_sample.append(item[:item.find('_')]) # if it is real sample, get the sample name by locate the first "_"
            else:
                if item.find(separator) != -1:
                    new_sample.append(item[:item.find(separator)]) # if it is derep unit, get the derep name by locagte the first ";"
                else:
                    new_sample.append(item)
        
        new_otu_map[new_key] = new_sample

    for element in centroid_list:
        fast_dict[element] = {}
        fast_dict[element]['seq'] = centroid_dict[element]
        
        fast_dict[element]['sample'] = {}
        for sample in new_otu_map[element]:
            try:
                fast_dict[element]['sample'][sample] += 1
            except KeyError:
                fast_dict[element]['sample'][sample] = 1
    
    return fast_dict
#%%    

def merge_fast_output(input_fast_derep, input_fast_otu):
# Merge the derep output with the OTU output
    hybrid_fast_dict = {}
    for key, value in input_fast_otu.items():
        
        hybrid_fast_dict[key] = {}
        hybrid_fast_dict[key]['seq'] = value['seq']
        hybrid_fast_dict[key]['sample'] = {}
        
        for derep_unit in value['sample'].keys():

            hybrid_fast_dict[key]['sample'][derep_unit] = {}
            hybrid_fast_dict[key]['sample'][derep_unit]['seq'] = input_fast_derep[derep_unit]['seq']
            hybrid_fast_dict[key]['sample'][derep_unit]['sample'] = input_fast_derep[derep_unit]['sample']

    
    return hybrid_fast_dict

#%%
def merge_derep_fast_with_Qiime_otu():
    pass

def write_fast_output(input_fast_dict, output_file):
# Save the FAST output in JSON format
    import json
    json.dump(input_fast_dict, open(output_file, 'w'))

def read_fast_output(input_fast_file):
# Load the FAST output
    import json
    return json.load(open(input_fast_file))

#%%
class fast_output_parser(object):
    def __init__(self, input_fast):
        self.fast = input_fast
        temp_check = input_fast[list(input_fast.keys())[0]]['sample']
        temp_value = type(temp_check[list(temp_check.keys())[0]])
        #print temp_value is dict
        if temp_value is int:
            self.fast_type = 'individual'
        elif temp_value is dict:
            self.fast_type = 'hybrid'
        #print self.fast_type
        self.unit_count = len(input_fast)
    
    def get_seqs(self, sort_by_size = True, sizeout = False):
        if sort_by_size:
            seq_list = []
            for key, value in self.fast.items():
                current_record = []
                current_record.append(key)
                current_record.append(value['seq'])
                
                current_size = 0
                
                if self.fast_type == 'individual':
                    current_size = value['sample'].values()
                    current_size = sum([int(i) for i in current_size])
                    current_record.append(current_size)                            
                
                if self.fast_type == 'hybrid':
                    current_size = 0
                    for unit, sub_value in value['sample'].items():
                        for sample, size in sub_value['sample'].items():
                            current_size += int(size)
                    current_record.append(current_size)
                
                
                seq_list.append(current_record)
            
            seq_list = sorted(seq_list, key=lambda x:x[2], reverse=True)
        
        else:
            seq_list = []
            for key, value in self.fast.items():
                current_record = []
                current_record.append(key)
                current_record .append(value)['seq']
                seq_list.append(current_record)
        
        if sizeout:
            for index, record in enumerate(seq_list):
                seq_list[index][0] = record[0] + ';size=' + str(record[2])
            return seq_list
        else:
            return seq_list
    
    def get_samples(self, sort_by_size = True):
        pass
    
    def detail_sample_unit(self, target_unit):
    # Get an output for a single sample unit        
        if self.fast_type == 'individual':
            print('This is not a hybrid FAST map.')
            return ''
        
        elif self.fast_type == 'hybrid':
            unit_dict = self.fast[target_unit]
            unit_list = []
            
            for derep_unit, value in unit_dict['sample'].items():
                unit_list += value['sample'].keys()
            unit_list = set(unit_list)
            unit_list = sorted(list(unit_list))
            
            header = [''] + list(unit_list) + ['Sequence']
            output_content = []                      
            
            
            for derep_unit, value in unit_dict['sample'].items():
                current_line = []
                current_line.append(derep_unit)
                for sample in unit_list:
                    try:
                        current_line.append(int(value['sample'][sample]))
                    except KeyError:
                        current_line.append(0)
                current_line.append(value['seq'])
                output_content.append(current_line)
            output_content = sorted(output_content, key=lambda x:sum(x[1:-1]), reverse=True)
            
            output_content = [header] + output_content
            return output_content
    
    def parse_otu_table(self):
        if self.fast_type == 'individual':
            print('This is not a hybrid FAST map.')
            return ''
        elif self.fast_type == 'hybrid':
            sample_list = []
            otu_table_dict = {}
            for otu, value in self.fast.items():
                otu_table_dict[otu] = {}
                for derep_unit, derep_value in value['sample'].items():
                    for sample, abundance in derep_value['sample'].items():
                        treatment = sample
                        if treatment not in sample_list:
                            sample_list.append(treatment)
                        try:
                            otu_table_dict[otu][treatment] += abundance
                        except KeyError:
                            otu_table_dict[otu][treatment] = abundance
            return sample_list, otu_table_dict
#%%