# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 13:07:49 2015
Random_subsample.py will randomly sample the giving abundance data to a certain
depth.

rarefaction() will rarefy a single sample to certain depth, returning a list of interger,
this corresponds to the OTU in the orginal table.

Samples with total sequence less than rarefy depth will be kept untouched.

repeat_rarefaction() will repeatly rarefy a single sample for a given time, and return
a list by species. It does not return the header.

repeat_rarefaction_parallel() will do the same job as repeat_rarefaction(), but using
multiple processor.

The number of processor can be checked by:
    import psutil
    psutil.cpu_count()

Input variable can be obtained by using parse_OTU() from a tab delimited OTU table.
otuid is a list of the name of OTUs
sample is a list wiht the abundance of each OTU
depth is the number of times for the random sampling
rep is the repeated time for a single sample

@author: Zewei Song
songzewei@outlook.com
"""
from __future__ import print_function
from __future__ import division
#%%##############################################################################
def rarefaction(sample,depth):
    sample_single = generate_repeat_otu(sample)
    sample_single_rare = rarefy_repeat_otu(sample_single,depth)
    sample_rarefied_corrected = degenerate_rarefied_sample(sample_single_rare)
    return sample_rarefied_corrected
#%%#############################################################################
def repeat_rarefaction(sample,depth,rep): #has problem
    sample_rarefied_corrected = []
    sample_single = generate_repeat_otu(sample)
    for i in range(rep):
        sample_single_rare = rarefy_repeat_otu(sample_single,depth)
        sample_rarefied_corrected.append(degenerate_rarefied_sample(sample_single_rare))
    return sample_rarefied_corrected
#%%#############################################################################
def repeat_rarefaction_parallel(sample,depth,rep,processor=4):
    if sum(sample) <= depth:
        return [sample]*rep
    else:
        from multiprocessing import Pool
        p = Pool(processor)
        worker_input = allocate_processor(sample,depth,rep,processor=processor)
        sample_rarefied = (p.map(worker_parallel,worker_input))
        p.close()
        sample_rarefied_cat = []
        for item in sample_rarefied:
            for element in item:
                sample_rarefied_cat.append(element)
        return sample_rarefied_cat
#%%#############################################################################
# Generate a new list that contains the repeat for each OTU.
# The variable sample should be a list of OTU abundance (include 0), without sample name
def generate_repeat_otu(sample):
    sample_single = []
    n = 0
    for item in sample:
        temp = []
        if int(item) > 0:
            for i in range(int(item)):
                temp.append(n)
            sample_single += temp
        n += 1
    sample_single = {'single':sample_single,'original_length':len(sample)} # Save original length of degenerate to keep the order of input OTU
    return sample_single
#%%###############################################################################
#Random sample the repeat list for a given time (depth)
def rarefy_repeat_otu(sample_single,depth):
    from random import randrange
    single = sample_single['single'][:] #if use sample_single as variable, w/o [:] will change it.
    original_length = sample_single['original_length']
    if  len(single) <= depth: # Should be len instead of sum here
        sample_single_rare = single #Total abundance <= depth, ignore rarefaction
    else:
        sample_single_rare = []
        for i in range(depth):
            pick_position = randrange(len(single))
            sample_single_rare.append(single[pick_position])
            single[pick_position] = single[-1]#switch position between picked element and the last one
            del single[-1]                           #to reduce calculation time
    sample_single_rare = {'single':sample_single_rare,'original_length':original_length}
    return sample_single_rare
#%%###############################################################################
def degenerate_rarefied_sample(sample_single_rare):
    from collections import Counter
    count = Counter(sample_single_rare['single'])
    sample_rarefied_corrected = [0] * sample_single_rare['original_length'] # Create a list with the length of original OTUs
    for item in count:
        sample_rarefied_corrected[int(item)] = count[item]
    return sample_rarefied_corrected
#%%##############################################################################
def worker_parallel(worker_input):
    sample = worker_input[0]
    depth = worker_input[1]
    rep = worker_input[2]
    return repeat_rarefaction(sample,depth,rep)
#%%##############################################################################
def allocate_processor(sample,depth,rep,processor=4):
#Generate input parameter for worker_parallel().
    rep_parallel = []
    worker_input = []
    for i in range(processor):
        rep_parallel.append(rep//processor)
    n = 0
    for i in range(rep%processor):
        rep_parallel[n] += 1
        n += 1
    for i in range(processor):
        worker_input.append([sample,depth,rep_parallel[i]])
    return worker_input
#%%##############################################################################
def generate_random_index(sample_size, random_size):
    import numpy as np
    index_list = range(sample_size)
    random_index_list = np.random.choice(index_list, random_size, replace=False)
    return random_index_list
#%%##############################################################################
if __name__ == "__main__":
    print("################################################################################")
    print("random_subsample.py")
    print("Zewei Song")
    print("University of Minnesota")
    print("songzewei@outlook.com\n")
    import psutil
    print("The number of processor in this computer is: "+str(psutil.cpu_count()))
    print("\n")
    otuid = []
    sample = []
    for i in range(20):
        otuid.append('OTU_%d'%(i+1))
        sample.append(100)
    print("Example otuid, length = 20, otu_id is not needed:")
    print(otuid)
    print()
    print("Example sample, length = 20:")
    print(sample)
    print()
    print("Rarefy sample to the level of 500 from a total of 2000:")
    print("rarefaction(sample,500)")
    print(rarefaction(sample,500))
    print()
    print("Rarefy sample to the level of 500 for 3 times:")
    print("repeat_rarefaction(sample,500,3)")
    repeat_sample = repeat_rarefaction(sample,500,3)
    for sample in repeat_sample:
        print(sample)
    print()
    print("Rarefy sample to the level of 500 for 3 times and using 4 processors:")
    print("repeat_rarefaction_parallel(sample,500,3,processor=4)")
    repeat_sample = repeat_rarefaction_parallel(sample,500,3,processor=4)
    for sample in repeat_sample:
        print(sample)
    print("################################################################################")