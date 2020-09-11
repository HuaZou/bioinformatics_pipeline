# -*- coding: utf-8 -*-
"""
Created on Fri May 01 14:27:45 2015

Please feel free to contact me for any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""

import argparse
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Input FASTA file to be dereplicated.')
parser.add_argument('-o', '--output', help='Name for output OTU map and FASTA file.')
parser.add_argument('-t', '--thread', default = 1, help='Number of threads to be used.')

args = parser.parse_args()

input_file = args.input
output_name = args.output
output_map = output_name + '.txt'
output_fasta = output_name + '.fasta'
thread = int(args.thread)

def dereplicate_worker(input_seqs, n, count):
    # Dereplicate a chunk of input sequences
    # n is the iterate number of the worker
    # count is a shared list on number of processed sequences by each worker
    # total is a shared value for total number of sequences
    derep = {}
    for record in input_seqs:
        try:
            derep[record[1]].append(record[0])
        except KeyError:
            derep[record[1]] = [record[0]]
        count[n] += 1
    dump_file = 'file_'+str(n)+'.txt'
    with open(dump_file, 'w') as f:
        pickle.dump(derep, f)
    
def divide_seqs(total, thread_num):
    # Set break point for input sequences
    seqs_divide = []
    size = total / thread_num
    for i in range(thread_num):
        seqs_divide.append(size)
    seqs_divide[0] += total % thread_num
    return seqs_divide

if __name__ == '__main__':
    import time
    from lib import File_IO
    from multiprocessing import Process, Manager
    import os, sys

    print 'Using %i threads ...' % thread

    input_file = input_file
    print 'Loading %s ...' % input_file
    seqs = File_IO.read_seqs(input_file)
    seqs_num = len(seqs)
    print 'Read in %i sequences.' % seqs_num

    # Separated seqs into pools
    print 'Separating raw sequences into %d jobs ...' % thread
    d = divide_seqs(seqs_num, thread)

    start = time.time()
    # Create shared list for store dereplicated dict and progress counter
    manager = Manager()
    count = manager.list([0] * thread)

    print 'Starting dereplicating ...'
    workers = []
    for i in range(thread):
        current_size = d[i]
        workers.append(Process(target=dereplicate_worker,
                               args=(seqs[0:current_size], i, count)))
        seqs[0:current_size] = []
    
    print 'Starting %i jobs ...' % thread  
    count_worker = 1
    for job in workers:
        job.start()
        print 'Starting thread No. %i ...' % count_worker
        count_worker += 1
    
    job_alive = True
    while job_alive:
        time.sleep(0.01)
        job_alive = False
        for job in workers:            
            if job.is_alive():
                job_alive = True
        progress = "Dereplicating: " + str(round(sum(count)/float(seqs_num)*100,2)) + "%" + "\r"
        sys.stderr.write(progress)

    for derep_worker in workers:
        derep_worker.join()
    print "Dereplicating: " + "100.00% \n"
    print 'Finished dereplication.'

    # Merged dereplicated dictionaries into a single dict
    file_list = []
    for i in range(thread):
        file_list.append('file_'+str(i)+'.txt')
    
    merged_dict = {}
    count_merge = 0
    for pickle_file in file_list:
        with open(pickle_file, 'rU') as f:
            temp = pickle.loads(f.read())
        os.remove(pickle_file)
        for key, value in temp.items():
                try:
                    merged_dict[key] += value
                except KeyError:
                    merged_dict[key] = value
                count_merge += 1
                merge_progress = 'Merging %i record ... \r' %count_merge
                sys.stderr.write(merge_progress)
                
    print
    print "Sequences dereplicated, clapsed from %i into %i sequences." % (seqs_num, len(merged_dict))
    s = [len(merged_dict[i]) for i in merged_dict]
    print 'Dereplicated OTU size: Max=%i, Min=%i, Average=%i.' % (max(s), min(s), round(float(sum(s) / len(s)), 2))
    end = time.time()
    used_time = round(end-start,2)
    print "Used time: " + str(used_time) + ' seconds.'
    print

    # Name the dereplicated group
    count = 0
    for key, value in merged_dict.items():  # Add group name to the end of the the name list of each group
        derep_name = 'derep_' + str(count)
        value.append(derep_name)
        count += 1

    # Output dereplicated FASTA file
    print 'Writing dereplicated sequence and OTU map ...'
    output_seq_file = output_fasta
    with open(output_seq_file, 'w') as f:
        for key, value in merged_dict.items():
            f.write('>%s\n' % value[-1])
            f.write('%s\n' % key)
    print '%s contains dereplicated sequences.' % output_fasta

    # Output Qiime style map
    with open(output_map, 'w') as f:
        for key, value in merged_dict.items():
            f.write('%s\t%s\n' % (value[-1], '\t'.join(value[:-1])))  # Use the last element as group name
    print '%s contains an OTU map for dereplicated sequences.' % output_map