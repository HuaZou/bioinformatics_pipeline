#!/usr/bin/env python 

import sys, os, re
import argparse as ap 

cwd = os.getcwd()

def parse_argument(args):
	parser = ap.ArgumentParser(description='cut adapter by cutadapt')
    parser.add_argument(
		'-i', '--id', metavar='sampleid', type=str,
		help='id', required=True)
	parser.add_argument(
		'-r1', '--read1', metavar='fastq1', type=str,
		help='fastq1', required=True)
	parser.add_argument(
		'-r2', '--read2', metavar='fastq2', type=str, default=None,
        help='fastq2')
	parser.add_argument(
		'-f', '--fwd', metavar='foward', type=str,
        help='foward sequence of adapter')  
	parser.add_argument(
		'-r', '--rev', metavar='reverse', type=str,
        help='foward sequence of adapter')
	parser.add_argument(
		'-b', '--base', metavar='match_bases', type=int,
        help='match_bases')  
	parser.add_argument(
		'-e', '--err', metavar='error_rate', type=int,
        help='error_rate') 
	parser.add_argument(
		'-m', '--min', metavar='min_length', type=int,
        help='min_length')
	parser.add_argument(
		'-o', '--outdir', metavar='directory', type=str,
		 help='results of output directory', required=True)

	return parser.parse_args()


def trim_adapter(id, read1, read2, fwd_seq, rev_seq, bases, error, min, dir):

    trim_FWD = dir + "/FWD"
    if not os.path.exists(trim_FWD):
        os.makedirs(trim_FWD)
    
    trim_REV = dir + "/REV"
    if not os.path.exists(trim_REV):
        os.makedirs(trim_REV)

    trim_file1 = trim_FWD + "/" + id + "_trimmed.1.fastq.gz"
    trim_file2 = trim_REV + "/" + id + "_trimmed.2.fastq.gz"      
    info_file   = "--info-file=" + outdir + "/" + id + "_reads.adapter.txt"

    shell_list = ["cutadapt", "-g", fwd_seq, "-G", rev_seq, "-O", bases, "-e", error, \
                        "-o", trim_file1, "-p", trim_file2, read1, read2,\
                        "--minimum-length", min, "--discard-untrimmed", info_file]
   
    shell = " ".join(str(x) for x in shell_list)
    os.system(shell)



def main():
    args = parse_argument(sys.argv)
    sampleid    = args.id
    fq1         = args.read1
    fq2         = args.read2
    fwd         = args.fwd
    rev         = args.rev
    match_bases = args.base 
    error_rate  = args.err 
    min_length  = args.min 
	out_dire    = args.outdir
    
    trim_adapter(sampleid, fq1, fq2, fwd, rev, match_bases,\
        error_rate, min_length, out_dire)

if __name__ == "__main__":
    main()
