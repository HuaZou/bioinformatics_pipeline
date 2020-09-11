# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 10:02:35 2015
Parse a tab delimited OTU table into its basic elements.
I used list instead of numpy array for simplicity.
It seems enough for current size of otu table.
@author: Zewei Song
"""
class parse_OTU(object):
    def __init__(self,filepath, meta_column='taxonomy'):
        #Use to store the entire table
        #metaSize specify the number of column of the meta Data (usually 1 for taxonomy)
        import sys
        self.otu_table = []      
        #Read in the entire OTU table        
        input = open(filepath, 'r')
        for line in input:
            self.otu_table.append(line.strip('\n').split('\t'))
        input.close()        
        
        self.header = self.otu_table[0]
        if meta_column in self.header:
            self.meta_start = self.header.index(meta_column)
        else:
            self.meta_start = len(self.header) + 1
        self.sample_id = self.header[1:self.meta_start]
        self.sample_number = len(self.sample_id)
        self.species_number = len(self.otu_table) -1
        
        self.otu = []
        self.meta = []
        self.otu_id = []   #need to change otu to int   
        for line in self.otu_table[1:]:
            self.meta.append(line[self.meta_start:]) #per species
            self.otu_id.append(line[0]) 
            try:
                self.otu.append(map(int,line[1:self.meta_start])) #per species
            except ValueError:
                print 'There is non-number value in the sample column of the otu table.'
                sys.exit()
    
            
    def get_sample(self): #Get a transposed otu abundance matrix, by sample. This will be the input of random_subsample.
                          #The returned sample will be a dict, use sample['abundance'] to get its abundance
                          #And sample['name'] for its name.
                          #The input for random_subsample is sample['abundance']
        abundance = [list(i) for i in zip(*self.otu)]
        sample = []
        n = 0
        for line in abundance:
            sample.append({'abundance':line,'name':self.sample_id[n]})
            n += 1
        return sample
        
    def get_species(self): #Get a otu abundance matrix, by species
        abundance = self.otu
        species = []
        n = 0
        for line in abundance:
            species.append({'abundance':line,'name':self.otu_id[n]})
        return species

'''Usage:''' 
if __name__ == '__main__':
    fomes = parse_OTU('fomes_otu.txt')
    sample = fomes.get_sample()
    species = fomes.get_species()