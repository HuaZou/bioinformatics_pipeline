# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 13:12:40 2015

Create a class object that iterate through the given 'fasta' or 'fastq' file.
There is an example on how to use it at the bottom. You can create a mock test.fasta file to try it.
It will check if the head symbol is matched to the format (fasta = ">"; fastq = "@").
For now, it will NOT check other features in fastq file, such as
    quality line header "+"
    same length for sequence and quality line
This is mainly for the reason of speed, and in my observations that miss matched in fastq is very rare.
Please feel free to contact me for any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
"""

class ParseSeq(object):
    def __init__(self,filename,filetype):
        self.file = filename
        if filetype == "fasta" or filetype == "Fasta":
            self.identifier = ">"
            self.range = 2
        elif filetype == "fastq" or filetype == "Fastq":
            self.identifier = "@"
            self.range = 4
        else:
            import sys
            print "Please specified the right format: fasta or fastq"
            sys.exit()
        self.file = open(filename,'rU')
        self.count = 0
    def __iter__(self):
        return self
    def next(self):
        currentRecord = []
        self.count+=1
        for i in range(self.range):
            line = self.file.readline()
            if line:
                currentRecord.append(line)
            else:
            #Stop iteration at the EOF.
                raise StopIteration
        #Raise error when the header is not matched to the right format.
        assert currentRecord[0][0] == self.identifier,\
            "The %dth record in %s does not start with %s." %(self.count,self.filename,self.identifier)
        return currentRecord

#%%
class ParseFasta(object):
    def __init__(self, filepath):
        self.file = open(filepath,'rU')
        self.count = 0
    def __iter__(self):
        return self
    def next(self):
        currentRecord = []
        self.count+=1
        for i in range(2):
            line = self.file.readline()
            if line:
                currentRecord.append(line[:-1])
            else:
                raise StopIteration
        currentRecord[0] = currentRecord[0][1:]
        return currentRecord
#%%            

if __name__ == "__main__":
    parser = ParseSeq('test.fasta','fasta')
    print "The identifier is: %s" %(parser.identifier)
    count = 0
    for record in parser:
        count += 1
        print "Record No. %d: "%count
        print "Sequence Name: %s" %record[0].strip('\n')
        print "Sequence: %s" %record[1].strip('\n')