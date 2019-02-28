import argparse
import pandas as pd

def concatenate_titer_files(titer_tsvs, output):
    concat = pd.DataFrame()
    for titer_tsv in titer_tsvs:
        df = pd.read_csv(titer_tsv, sep='\t', header=None).fillna('')
        concat = pd.concat([concat, df], ignore_index=True)

    concat.to_csv(output, sep='\t', header=None, index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Concatenate egg and cell titers tsv files for each assay type",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--titer_tsvs', nargs='+', help='tsv files to concatenate')
    parser.add_argument('--output', help="name of the file to write concatenated tsv to")
    args = parser.parse_args()

    concatenate_titer_files(args.titer_tsvs, args.output)
