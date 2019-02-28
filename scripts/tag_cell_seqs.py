import argparse
import re
from Bio import SeqIO


def tag_cell_seqs(in_file, out_file):
    handle = open(out_file, 'w')
    for record in SeqIO.parse(in_file, "fasta"):
        if 'cell' in record.id:
            new_id = re.sub(r'(?=\|(flu)\|)', r'-cell', record.id)
            record.id = new_id
        handle.write(">"+record.id+'\n')
        handle.write(str(record.seq) + "\n")
    handle.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Append -cell to strain name in fasta file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--in_file', help='path to fasta file to edit')
    parser.add_argument('--out_file', help='path to edited fasta file')
    args = parser.parse_args()

    tag_cell_seqs(args.in_file, args.out_file)
