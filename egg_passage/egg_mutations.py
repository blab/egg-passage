"""
Find egg-passaging-specific mutations
"""

import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO

def find_mutation_prevalence(prefix):

    df = pd.read_csv('data/'+prefix+'_df.csv')

    mut_sites = [col for col in df.columns if col[0:3]=='mut']

    mut_prev = {}

    #Find prevalence of each of the top mutations
    for pas_type in df['passage'].unique():
        mut_prev_pas = {}
        for mut_site in mut_sites:
            mut_count = float(len(df[df.passage == pas_type & df[mut_site]==1]))/float(len(df[df.passage == pas_type]))
            mut_prev_pas[mut_site] = mut_count
        mut_prev[pas_type] = mut_prev_pas

    mut_prev_df = pd.DataFrame(mut_prev)

    print(mut_prev_df)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Determines egg-specific mutations")
    parser.add_argument('--prefix', default= 'h3n2', help= "specify prefix for naming data files")
    args = parser.parse_args()

    find_mutation_prevalence(prefix = args.prefix)
