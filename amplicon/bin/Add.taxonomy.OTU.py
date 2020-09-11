#!/usr/bin/env python
import sys

def Usage(command_args):
    print ("Usage:")
    print (f"python {command_args[0]} [OTU] [taxonomy] [out]")
    sys.exit()

def OTU_process(file_name,mytaxon,out_name):
    with open(file_name,'r') as fh_otu, open(out_name,'w') as fh_out:
        header = fh_otu.readline().rstrip()
        fh_out.write(f"{header}\tTaxonomy\n")
        for line in fh_otu:
            line = line.rstrip()
            tmp = line.split("\t")
            if tmp[0] in mytaxon:
                fh_out.write(f"{line}\t{mytaxon[tmp[0]]}\n")
            else:
                print(f"{tmp[0]} not in taxonomy files")
                sys.exit()

def Read_taxon(file_name):
    res = {}
    with open(file_name,'r') as fh:
        for line in fh:
            line = line.rstrip().split("\t")
            if line[0] not in res:
                res[line[0]] = line[1]
            else:
                print(f"{line[0]} replicated in taxonomy file")
                sys.exit()
    return res

if __name__ == "__main__":
    if len(sys.argv) != 4:
        Usage(sys.argv)
    # args parser
    otu_file_name,taxon_name,out=sys.argv[1:]
    # read files
    taxon_dat = Read_taxon(taxon_name)
    OTU_process(otu_file_name,taxon_dat,out)