#!/usr/bin/env python

import sys, re, os 
import argparse as ap 


def parse_argument(args):
	parser = ap.ArgumentParser(description='fq scan by fastqc')
	parser.add_argument(
		'-r1', '--read1', metavar='fastq1', type=str,
		help='fastq1', required=True)
	parser.add_argument(
		'-r2', '--read2', metavar='fastq2', type=str, default=None,
        help='fastq2')
	parser.add_argument(
		'-t', '--threads', metavar='threads', type=int, default=2,
        help='the number of threads')        
	parser.add_argument(
		'-o', '--out', metavar='directory', type=str,
		 help='results of output directory', required=True)

	return parser.parse_args()


def fastqc(read1, read2, threads, directory):
    if not read2 :
        shell_lst = ["fastqc", "-o", directory, "-t", threads, read1, read2]
    else:
        shell_lst = ["fastqc", "-o", directory, "-t", threads, read1]
    shell_fqc = " ".join(str(x) for x in shell_lst)

    os.system(shell_fqc)


def main():
    args = parse_argument(sys.argv)
    fastqc(args.read1, args.read2, args.threads, args.out)


if __name__ == "__main__":
    main()
