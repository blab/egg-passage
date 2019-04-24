"""
Compare cTiterSub between pairs of egg-passaged/non-egg-passaged strains
"""

import argparse, ast
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_ctitersub(input_df):
    df_name = str.split(input_df, 'dataframes/')[1]
    prefix = str.split(df_name, '.csv')[0]
    pre_prefix = str.split(prefix, '_hi')[0]

    muts_to_consider = ['T160K', 'G186V', 'L194P']
    titer_diffs_singlemut = []
    titer_diffs = []
    for assay in ['hi', 'fra']:
        all_df = pd.read_csv('dataframes/'+pre_prefix+"_"+assay+'.csv')
        df = all_df[all_df['pair_id']!=0]

        #Re-organize DF to one row per pair
        sub_egg = df[df['passage']=='egg'][['source', 'egg_muts', 'cTiterSub']].rename(columns = {'cTiterSub':'egg_titer'})
        sub_u = df[df['passage']=='unpassaged'][['source', 'cTiterSub']].rename(columns = {'cTiterSub':'pair_titer'})
        sub_cell = df[df['passage']=='cell'][['source', 'cTiterSub']].rename(columns = {'cTiterSub':'pair_titer'})

        pairs_u_df = sub_egg.merge(sub_u)
        pairs_cell_df = sub_egg.merge(sub_cell)
        pairs_cell_u_df = sub_u.merge(sub_cell)
        pairs_df = pairs_u_df.merge(pairs_cell_df, how='outer')

        for k,v in pairs_df.iterrows():
            if len(ast.literal_eval(v['egg_muts']))==1:
                diff = v['egg_titer'] - v['pair_titer']
                titer_diffs_singlemut.append({'mutation':'overall', 'titer_diff': diff, 'assay': assay})
                for egg_mut in muts_to_consider:
                    if egg_mut in ast.literal_eval(v['egg_muts']):
                        titer_diff = v['egg_titer'] - v['pair_titer']
                        titer_diffs_singlemut.append({'mutation': egg_mut, 'titer_diff': titer_diff, 'assay': assay})
            else:
                diff = v['egg_titer'] - v['pair_titer']
                titer_diffs.append({'mutation':'overall', 'titer_diff': diff, 'assay': assay})
                for egg_mut in muts_to_consider:
                    if egg_mut in ast.literal_eval(v['egg_muts']):
                        titer_diff = v['egg_titer'] - v['pair_titer']
                        titer_diffs.append({'mutation': egg_mut, 'titer_diff': titer_diff, 'assay': assay})

    titer_diffs_df = pd.DataFrame(titer_diffs)
    titer_diffs_singlemut_df = pd.DataFrame(titer_diffs_singlemut)

    #Total antigenic change in viruses that have these mutations
    fig, ax = plt.subplots()
    fig = sns.barplot(x='mutation', y='titer_diff', hue='assay', order=['overall', 'T160K', 'L194P', 'G186V'], data=titer_diffs_df, palette= ['#d73027', '#fdae61'])
    ax.set(xlabel = '', ylabel='log2(egg titer)-log2(paired titer)')
    fig.get_figure().savefig('plots/'+str(prefix)+'/titer_diffs_'+str(pre_prefix)+'.pdf', bbox_inches='tight')

    #Antigenic effect of single mutation
    fig2, ax2 = plt.subplots()
    fig2 = sns.barplot(x='mutation', y='titer_diff', hue='assay', order=['overall', 'T160K', 'L194P', 'G186V'], data=titer_diffs_singlemut_df, palette= ['#d73027', '#fdae61'])
    ax2.set(xlabel = '', ylabel='log2(egg titer)-log2(paired titer)')
    fig2.get_figure().savefig('plots/'+str(prefix)+'/titer_diffs_singlemut_'+str(pre_prefix)+'.pdf', bbox_inches='tight')
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Plot titer differences between paired strains")
    parser.add_argument('--in_file', help= "input dataframe file")
    args = parser.parse_args()

    plot_ctitersub(input_df = args.in_file)
