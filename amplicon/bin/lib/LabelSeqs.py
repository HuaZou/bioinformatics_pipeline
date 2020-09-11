# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 15:22:43 2015

Rename sequences in the fastq or fasta file to real sample name.

It will add a Qiime style label (SampleName_132) at first, following a customized label and a USEARCH style label.
The sequence head should looks like:
    @SampleA_1;SampleName=SampleA-1-R1;barcodelabel=SampleA;
For default it will add both Qiime and USEARCH style label. As far as I tested it, Qiime parse the OTU map by seaching
the first "_", so there should be no conflict with other labels.

The program will stop if the folder labeled already exist.

Please feel free to contact me for any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""
from __future__ import print_function
from __future__ import division
# %%+++++Single thread relabeling++++++++++++++++++++++++++++++++++++++++++++++++++++++
def ReLabelFastQ(file_name, label, read_type, input_folder, output_folder='labeled', file_type='fastq',
                 label_type='qiime'):
    #%% Read in sequence file and change the header
    from lib import File_IO

    file_content = File_IO.read_seqs(input_folder + '/' + file_name, file_type=file_type)
    head_symbol = '@'
    if len(file_content[0]) == 2:
        head_symbol == '>'

    count = 0
    for record in file_content:  #Loop through header of the records
        record[0] = ChangeName(label, count, read_type, label_type=label_type)
        count += 1

    file_labeled = output_folder + '/labeled_' + file_name
    with open(file_labeled, 'w') as f:
        for record in file_content:
            record[0] = head_symbol + record[0]  # Add head symbol to sequence name
            for line in record:
                f.write('%s\n' % line)
    return count


def ParseMapping(mapping_file, input_folder):
    # Return a dictionary containing filename, label and read type
    with open(mapping_file, 'rU') as f:
        mapping_content = f.readlines()
        header = mapping_content[0].strip('\n').split('\t')
        pos_SampleID = header.index('#SampleID')
        pos_InputFileName = header.index('InputFileName')
        pos_ReadType = header.index('ReadType')
        temp = mapping_content[1:]
        mapping = []
        for record in temp:
            sample = record.strip('\n').split('\t')
            mapping.append(
                {'file': sample[pos_InputFileName], 'input_folder': input_folder + '/', 'label': sample[pos_SampleID], 'read_type': sample[pos_ReadType]})
        return tuple(mapping)


def ChangeName(barcodelabel, count_seq, read_type='R1', label_type='qiime'):
    # Change the sequence name
    count_seq = str(count_seq)
    usearch_label = 'barcodelabel=' + barcodelabel
    qiime_label = barcodelabel + '_' + count_seq
    if label_type == 'usearch':
        new_header = usearch_label + ';' + 'count=' + count_seq + 'read_type=' + read_type + ';'
        return new_header
    elif label_type == 'qiime':
        new_header = qiime_label
        return new_header
    elif label_type == 'both':
        new_header = qiime_label + ';' + usearch_label + ';' + 'read_type=' + read_type + ';'
        return new_header


#####################################################################################

#%%++++Multiple threads relabeling+++++++++++++++++++++++++++++++++++++++++++++++++++
def LabelFiles(mapping):
    # Received splitted mapping file and other parameters.
    for item in mapping:
        count = ReLabelFastQ(item['file'], item['label'], item['read_type'], item['input_folder'], \
                             output_folder=item['output_folder'], file_type=item['file_type'],
                             label_type=item['label_type'])
        print("%s sequences in %s relabeled to %s as %s file.\n" % (
            count, item['file'], item['label'], item['read_type']))


def CreateWorker(mapping_multithreads, threads=4):
    # Assigned multithreads mapping and parameters to the worker
    from multiprocessing import Process
    worker = []
    for i in range(threads):
        worker.append(Process(target=LabelFiles, args=(mapping_multithreads[i],)))
    return worker


def SplitMapping(mapping_file, input_folder, output_folder='labeled', file_type='fastq', label_type='both',
                 processor=4):
    # Allocate the mapping variable evenly to all processors
    mapping = ParseMapping(mapping_file, input_folder)
    for item in mapping:
        item['output_folder'] = output_folder
        item['file_type'] = file_type
        item['label_type'] = label_type
    file_number = len(mapping)
    threads = []
    for i in range(processor): threads.append(file_number // processor)
    for i in range(file_number % processor): threads[i] += 1

    temp = list(mapping[:])
    n = 0
    mapping_multithreads = []
    for i in range(processor): mapping_multithreads.append([])
    for task in threads:
        for i in range(task):
            mapping_multithreads[n].append(temp.pop(0))
        n += 1
    return tuple(mapping_multithreads)


#####################################################################################

#%%++++++++Main Function+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def MainLabelFiles(mapping_file, input_folder, threads=1, output_folder='labeled', file_type='fastq',
                   label_type='both'):
    #Create a new folder for relabeled files    
    from lib import File_IO

    File_IO.mk_dir(output_folder)
    if threads == 1:
        print("Relabeling files using %d thread ..." % threads)
        mapping = ParseMapping(mapping_file, input_folder)
        file_num = len(mapping)
        for item in mapping:
            count = ReLabelFastQ(item['file'], item['label'], item['read_type'], item['input_folder'], \
                                 output_folder=output_folder, file_type=file_type, label_type=label_type)
            print("%s sequences in %s relabeled to %s as %s file.\n" % (
                count, item['file'], item['label'], item['read_type']))

    elif threads > 1:
        print("Relabeling files using %d threads ..." % threads)
        mapping_multithreads = SplitMapping(mapping_file, input_folder, output_folder=output_folder,
                                            file_type=file_type, label_type=label_type, processor=threads)
        file_num = sum([len(i) for i in mapping_multithreads])
        worker = CreateWorker(mapping_multithreads, threads=threads)

        for item in worker:
            # Start workers
            item.start()
        for item in worker:
            # Wait until all workers finishes
            item.join()

    else:
        print("The number of threads cannot be negative.")
        import sys

        sys.exit()
    return file_num

#%%