import argparse
from Bio import SeqIO


def find_egg_seqs(fasta_file, references_file, force_include_file):

    handle = open(force_include_file, 'w+')

    with open(references_file, 'r') as ref_handle:
        ref_lines = ref_handle.readlines()
        for line in ref_lines:
            handle.write(str(line)+'\n')

    egg_seqs = []

    #Find all egg-passaged sequences
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        virus_name = str.split(seq_record.id, '|')[0]

        if 'egg' in virus_name:
            egg_seqs.append(virus_name)

    #Find all pairs for egg-passaged sequences and add to list of viruses to include during augur sampling
    for egg_seq in egg_seqs:
        seq_name = str.split(egg_seq,'-egg')[0]

        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            if seq_name in seq_record.id:
                handle.write(str.split(seq_record.id, '|')[0]+'\n')
    handle.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Find all egg-passaged strains and their cell-passaged and unpassaged matched strains",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--fasta_file', help='path to fasta file')
    parser.add_argument('--references_file', help='path to references file')
    parser.add_argument('--force_include_file',
                        help='file to save names of references and all egg-passaged sequences to')
    args = parser.parse_args()

    find_egg_seqs(args.fasta_file, args.references_file, args.force_include_file)
