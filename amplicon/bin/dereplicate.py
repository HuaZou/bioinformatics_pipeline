# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 18:18:32 2015

Dereplicate a give FASTA file. All record with identical sequence will be grouped into one OTU
in the output OTU map.

Technically, dereplication is OTU clustering with similarity = 1.0

This script is slower than USEARCH 8.0, but does not have limitation on memory usage.

Under windows, in order to use multiple thread, you need to run this program stand alone,
i.e. you can not summon this function using fast.py, but need to run dereplicate.py instead.

Please feel free to contact me for any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""
from __future__ import print_function
from __future__ import division

def dereplicate_worker(input_seqs, output_derep, n, count):
    # Dereplicate a chunk of input sequences
    # n is the iterate number of the worker
    # count is a shared list on number of processed sequences by each worker
    # total is a shared value for total number of sequences
    derep = output_derep[n]
    for record in input_seqs:
        try:
            derep[record[1]].append(record[0])
        except KeyError:
            derep[record[1]] = [record[0]]
        count[n] += 1
    output_derep[n] = derep

def dereplicate_single_thread(input_seqs):
    #import sys
    derep = {}
    count = 0
    for record in input_seqs:
        try:
            derep[record[1]].append(record[0])
        except KeyError:
            derep[record[1]] = [record[0]]
        count += 1
        #sys.stderr.write('Dereplicating %i seq ... \r' % count)
    return derep

def get_treatment(input_string):
    #import sys
    treatment = input_string[:input_string.find('_')]
    return treatment


def divide_seqs(total, thread_num):
    # Set break point for input sequences
    seqs_divide = []
    start = 0
    size = int(total / thread_num)
    for i in range(thread_num):
        seqs_divide.append([start, start + size])
        start += size
    seqs_divide[-1][-1] += total % thread_num
    return seqs_divide

def main(Namespace):
    import argparse
    import textwrap

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
                                    ------------------------
                                    By Zewei Song
                                    University of Minnesota
                                    Dept. Plant Pathology
                                    songzewei@outlook.com
                                    ------------------------'''), prog='fast.py -dereplicate')
    parser.add_argument('-i', '--input', help='Input FASTA file to be dereplicated.')
    parser.add_argument('-o', '--output', help='Name for output OTU map and FASTA file.')
    parser.add_argument('-t', '--thread', default = 1, help='Number of threads to be used.')
    parser.add_argument('-fast', default = "", help="Name of FAST style output file.")
    parser.add_argument('-sizeout', action = 'store_true', help='Specify to add a USEARCH style size label: ";szie=XXX"')

    args = parser.parse_args(Namespace)

    input_file = args.input
    output_name = args.output
    output_map = output_name + '.txt'
    output_fasta = output_name + '.fasta'

    thread = int(args.thread)

    import time
    from lib import File_IO
    from multiprocessing import Process, Manager
    import sys

    print('Using %i threads ...' % thread)
    start = time.time()

    input_file = input_file
    print('Loading %s ...' % input_file)
    seqs = File_IO.read_seqs(input_file)
    seqs_num = len(seqs)
    print('Read in %i sequences.' % seqs_num)

    # Disable multiprocess if using single thread
    if thread == 1:
        derep_dict = dereplicate_single_thread(seqs)
    else:
        # Separated seqs into pools
        print('Separating raw sequences into %d jobs ...' % thread)
        d = divide_seqs(seqs_num, thread)


        # Create shared list for store dereplicated dict and progress counter
        manager = Manager()
        derep_dict = manager.list([{}] * thread)
        count = manager.list([0] * thread)

        print('Starting dereplicating ...')
        workers = []
        for i in range(thread):
            current_range = d[i]
            workers.append(Process(target=dereplicate_worker,
                                   args=(seqs[current_range[0]:current_range[1]], derep_dict, i, count)))
        del seqs

        print('Starting %i jobs ...' % thread)
        count_worker = 1
        for job in workers:
            job.start()
            print('Starting thread No. %i ...' % count_worker)
            count_worker += 1

        job_alive = True
        while job_alive:
            time.sleep(0.01)
            job_alive = False
            for job in workers:
                if job.is_alive():
                    job_alive = True
            #progress = "Dereplicating: " + str(round(sum(count)/float(seqs_num)*100,2)) + "%" + "\r"
            #sys.stderr.write(progress)

        for derep_worker in workers:
            derep_worker.join()
        print('Finished dereplicating.')
        seqs = []  # Empty sequences list to free memory.

    # Merged dereplicated dictionaries into a single dict
    sys.stderr.write('\n')

    if thread > 1:
        sys.stderr.write('Merging %i dictionaries into one ...' % len(derep_dict))
        merged_dict = {}
        count = 0
        for d in derep_dict:
            for key, value in d.items():
                count += 1
                try:
                    merged_dict[key] += value
                except KeyError:
                    merged_dict[key] = value
                #sys.stderr.write('Merging %i sequence ...' % count + '\b' * 50,)
            derep_dict[0] = ''  # Empty finished dictionary to free memory.
    else:
        merged_dict = derep_dict
    print
    print("Sequences dereplicated, clasped from %i into %i sequences." % (seqs_num, len(merged_dict)))
    s = [len(merged_dict[i]) for i in merged_dict]
    print('Dereplicated OTU size: Max=%i, Min=%i, Average=%i.' % (max(s), min(s), round(float(sum(s) / len(s)), 2)))
    end = time.time()
    print("Used time: " + str(end - start) + ' seconds.')
    print


    # Name the dereplicated group
    size_list = sorted([[len(merged_dict[i]), i] for i in merged_dict], reverse=True)
    count = 0
    for element in size_list:
        derep_name = 'derep_' + str(count)
        element.append(derep_name)
        count += 1


    # Output dereplicated FASTA file
    print('Writing dereplicated sequence and OTU map ...')
    output_seq_file = output_fasta
    with open(output_seq_file, 'w') as f:
        if args.sizeout:
            for element in size_list:
                output_label = element[2] + ";size=" + str(element[0])
                f.write('>%s\n' % output_label)
                f.write('%s\n' % element[1])

        else:
            for element in size_list:
                output_label = element[2]
                f.write('>%s\n' % output_label)
                f.write('%s\n' % element[1])

    print('%s contains dereplicated sequences.' % output_fasta)

    # Output Qiime style map
    with open(output_map, 'w') as f:
        for element in size_list:
            name_list = merged_dict[element[1]]
            f.write('%s\t%s\n' % (element[2], '\t'.join(name_list)))  # Use the last element as group name
    print('%s contains an OTU map for dereplicated sequences.' % output_map)

    # Generate FAST style derep output file (a single file with sample names, counts, and dereplicated sequences)
    if args.fast != "":
        fast_file = args.fast

        fast_dict = {}

        for element in size_list:

            fast_dict[element[2]] = {} # Crearte a new dict for current derep unit
            fast_dict[element[2]]['seq'] = element[1] # Save dereplicated sequence

            sample_dict = {} # Create a dict for sample sequence count
            name_list = merged_dict[element[1]]
            for sample in name_list:
                current_sample = get_treatment(sample)
                try:
                    sample_dict[current_sample] += 1
                except KeyError:
                    sample_dict[current_sample] = 1

            fast_dict[element[2]]['sample'] = sample_dict

        import json
        json.dump(fast_dict, open(fast_file, "wb"))

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])