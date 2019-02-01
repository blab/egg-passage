"""
Analysis of strains where sequences are available for the virus before and after egg-passaging (or after cell-passaging and after egg-passaging)
"""

import argparse
import ast
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def find_paired_mutations(prefix, positions):
    """
    Directly compare sequences of viruses before and after egg-passaging
    """

    df = pd.read_csv('data/'+prefix+'_df.csv')
    #filter data for only paired sequences
    df = df[df['pair_id'].notnull()]

    positions = ast.literal_eval(positions)

    #Re-organize DF to one row per pair
    sub_egg = df[df['passage']=='egg'][['source'] + [str(pos) for pos in positions]].rename(columns = dict((str(pos), (str(pos)+'_egg')) for pos in positions))
    sub_u = df[df['passage']=='0'][['source'] + [str(pos) for pos in positions]].rename(columns = dict((str(pos), (str(pos)+'_u')) for pos in positions))
    sub_cell = df[df['passage']=='cell'][['source'] + [str(pos) for pos in positions]].rename(columns = dict((str(pos), (str(pos)+'_cell')) for pos in positions))

    pairs_u_df = sub_egg.merge(sub_u)
    pairs_cell_df = sub_egg.merge(sub_cell)
    pairs_cell_u_df = sub_u.merge(sub_cell)
    pairs_df = pairs_u_df.merge(pairs_cell_df, how='outer')

    #Find proportion of paired strains where egg is diff than pair
    egg_u_diff = {}
    egg_cell_diff = {}
    cell_u_diff = {}

    for pos in positions:
        pairs_u_df[str(pos)+'_diff'] = np.where((pairs_u_df[str(pos)+'_egg'] != pairs_u_df[str(pos)+'_u']), 1, 0)
        pairs_cell_df[str(pos)+'_diff'] = np.where((pairs_cell_df[str(pos)+'_egg'] != pairs_cell_df[str(pos)+'_cell']), 1, 0)
        pairs_cell_u_df[str(pos)+'_diff'] = np.where((pairs_cell_u_df[str(pos)+'_u'] != pairs_cell_u_df[str(pos)+'_cell']), 1, 0)

    for pos in positions:
        egg_u_diff[str(pos)] = float(len(pairs_u_df[pairs_u_df[str(pos)+'_diff'] == 1]))/float(len(pairs_u_df))
        egg_cell_diff[str(pos)] = float(len(pairs_cell_df[pairs_cell_df[str(pos)+'_diff'] == 1]))/float(len(pairs_cell_df))
        cell_u_diff[str(pos)] = float(len(pairs_cell_u_df[pairs_cell_u_df[str(pos)+'_diff'] == 1]))/float(len(pairs_cell_u_df))
    egg_u_diff_df = pd.DataFrame(egg_u_diff, index=['egg_unpassaged'])
    egg_cell_diff_df = pd.DataFrame(egg_cell_diff, index=['egg_cell'])
    cell_u_diff_df = pd.DataFrame(cell_u_diff, index=['cell_unpassaged'])

    plot_df = pd.concat([egg_u_diff_df, egg_cell_diff_df, cell_u_diff_df])
    plot_df = plot_df.unstack().reset_index().rename(columns={'level_0':'site', 'level_1':'comparison', 0:'prevalence'})

    sns.set(style="white")
    fig, ax = plt.subplots()
    sns.barplot(x='site', y='prevalence', hue='comparison', data=plot_df)
    ax.set(xlabel='HA position', ylabel='proportion of pairs with different genotypes')
    fig.suptitle('Direct comparison of single strains sequenced under multiple conditions')
    fig.savefig('plots/paired_comparisons_'+str(prefix)+'.pdf', bbox_inches='tight')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Analyze strains that were sequenced before and after egg-passaging")
    parser.add_argument('--prefix', default= 'h3n2_6y_hi', help= "specify prefix for naming data files")
    parser.add_argument('-pos', '--positions', default = '[160, 194, 186, 225, 219, 203, 156, 138]', help="specify a list of HA1 positions to analyze")
    args = parser.parse_args()

    find_paired_mutations(prefix = args.prefix, positions = args.positions)
