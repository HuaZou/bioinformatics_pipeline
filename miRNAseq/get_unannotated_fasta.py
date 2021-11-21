import sys
import argparse as ap 


def parse_argument(args):
        parser = ap.ArgumentParser(description='make blast index of ncRNA in Mus musculus')
        parser.add_argument(
                '-f', '--fasta', metavar='<fasta>', type=str,
                help='fasta file', required=True)
        parser.add_argument(
                '-d', '--id', metavar='<id>', type=str,
                 help='idlist file', required=True)
        parser.add_argument(
                '-o', '--out', metavar='<output>', type=str,
                help='output prefix', required=True)

        return parser.parse_args()


def main():
    args = parse_argument(sys.argv)
    outf = open(args.out, 'w')
    dict = {}
    with open(args.fasta, 'r') as fastaf:
        for line in fastaf:
            if line.startswith('>'):
                name = line.strip().split()[0][1:]
                dict[name] = ''
            else:
                dict[name] += line.replace('\n','')
     
    with open(args.id,'r') as listf:
        for row in listf:
            row = row.strip()
            for key in dict.keys():
                if key == row:
                    outf.write('>' + key+ '\n')
                    outf.write(dict[key] + '\n')
    outf.close()


if __name__ == "__main__":
        main()
