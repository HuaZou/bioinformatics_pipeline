#!/usr/bin/python

__author__ = ('ZouHua (zouhua@genomics.cn)')
__version__ = '0.1'
__date__ = '20190806'

import sys
import re
import os
import argparse as ap
cwd = os.getcwd()
result = "/".join([cwd, "/result/01.trimmed/"])

def parse_arguments(args):
    parser = ap.ArgumentParser(
        description="DESCRIPTION\n"
        "generate_qiime2_config version "+__version__+" ("+__date__+"): \n"
        "generate_qiime2_config.py \n"
        "AUTHORS: "+__author__+"\n\n"
        "",
        formatter_class=ap.RawTextHelpFormatter,
        prog="generate_qiime2_config.py")
    parser.add_argument(
        '-p', '--path', metavar='<filepath>', type=str,
        help="trimmed fastq directory path\n",
        required=True)
    parser.add_argument(
        '-o', '--out', metavar='<out>', type=str,
        help="table of qiime2 configure file\n",
        required=True)     
    parser.add_argument(
        '-v', '--version', action='version',
        version="generate_qiime2_config version {} ({})".format(__version__, __date__),
        help="generate_qiime2_config.py version and exit")
    return parser.parse_args()


def make_table(path, out):
    outf = open(out, "w")
    head = "\t".join(["sampleid", "forward-absolute-filepath", "reverse-absolute-filepath"])
    outf.write(head+"\n")
    list_dir = os.listdir(path)
    for line in list_dir:
        if re.findall(r'(\S+)\_1\.trimmed.fq.gz', line):
            sampleid = re.findall(r'(\S+)\_1\.trimmed.fq.gz', line)[0]
            #file_path = os.path.abspath(line)
            file_path = os.path.join(result, line)
            print(file_path)
            file_path2 = file_path.replace("1.trimmed.fq.gz", "2.trimmed.fq.gz")
            res = "\t".join([sampleid, file_path, file_path2])
            outf.write(res+"\n")
    outf.close()


def main():
    args = parse_arguments(sys.argv)
    make_table(args.path, args.out)
    print("Perfectly Running")


main()
