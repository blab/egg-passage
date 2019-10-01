"""
Find background genotype specifity of egg-passaging mutations
"""

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

def plot_kkclade_heatmap_mut(prefix):
    """
    Plot heatmap based on mutation, not site
    """
    df = pd.read_csv('dataframes/'+prefix+'_egg.csv')

    egg_mutations = ['A138S', 'H156R', 'H156Q', 'T160K', 'G186V', 'L194P', 'T203I', 'S219F', 'S219Y', 'D225G', 'N246H']

    for egg_mut in egg_mutations:
        if egg_mut[1:4] not in df.columns:
            egg_mutations.remove(egg_mut)

    #Remove clades with less than 5 egg seqs from analysis
    clade_size = df.groupby('kk_clade').size().reset_index().rename(columns={0:'kk_clade_total'})
    for k, v in clade_size.iterrows():
        if v['kk_clade_total']<=5:
            index_names = df[df['kk_clade']==v['kk_clade']].index
            df.drop(index_names, inplace=True)

    #Find enrichment of egg-passaged mutation at each site in each clade
    clade_enrichment = []
    for egg_mut in egg_mutations:
        egg_mut_site = egg_mut[1:-1]
        overall_pct = len(df[df['aa_mut'+str(egg_mut_site)]==egg_mut])/len(df)
        clade_count = df[df['aa_mut'+str(egg_mut_site)]==egg_mut].groupby('kk_clade').size()
        clade_pct = df[df['aa_mut'+str(egg_mut_site)]==egg_mut].groupby('kk_clade').size()/df.groupby('kk_clade').size()
        clade_count = clade_count.reset_index().rename(columns={0:'count'})
        #Fill 0-count clades with small number, for log
        clade_pct = clade_pct.reset_index().rename(columns={0:'pct'}).fillna(0.000001)
        site_clade = pd.merge(clade_count, clade_pct, how='right', on='kk_clade').fillna(0)
        site_clade['enrichment'] = site_clade['pct']/overall_pct
        site_clade['log_enrichment'] = np.log2(site_clade['enrichment'])
        #Replace small number with 0.0 for 0-count clades
        site_clade['pct'] = np.where(site_clade['pct']==0.000001, 0.0, site_clade['pct'])
        site_clade['site'] = (egg_mut+'\n ( '+str(round(overall_pct, 2))+' )')
        clade_enrichment.append(pd.DataFrame(site_clade))

    clade_enrichment = pd.concat(clade_enrichment, axis=0)
    #Reorder clades and mutations
    clade_enrichment.kk_clade = pd.CategoricalIndex(clade_enrichment.kk_clade,
                                                 categories= [v['kk_clade'] for k,v in (df.groupby('kk_clade')['date'].mean().reset_index().sort_values(by='date')).iterrows()])
    clade_enrichment.sort_values(by='kk_clade', inplace=True)

    site_order = sorted([x for x in list(clade_enrichment['site'].unique())], key = (lambda x: x[1:4]))
    clade_enrichment.site = pd.CategoricalIndex(clade_enrichment.site,
                                                 categories= site_order)
    clade_enrichment.sort_values(by='site', inplace=True)

    clade_heatmap = clade_enrichment.pivot('kk_clade', 'site', 'log_enrichment')
    heatmap_annotate = clade_enrichment.pivot('kk_clade', 'site', 'pct')

    fig, ax= plt.subplots(figsize=(10,20))
    sns.set(style="white")

    fig = sns.heatmap(clade_heatmap,
                      cmap=sns.diverging_palette(220, 20, sep=80, as_cmap=True),
                      vmax=3.25, vmin=-3.25, center = 0.0, cbar_kws={'label': 'log2 enrichment'})
    fig.set_xlabel('HA1 mutation \n (overall frequency in egg-passaged viruses)',fontsize=14)
    fig.set_ylabel('Clade',fontsize=14)
    fig.get_figure().savefig('plots/'+str(prefix)+'/testing_genetic_background_kkclade_heatmap_mut_'+str(prefix)+'.pdf', bbox_inches='tight')

def main(input_df):
    df_name = str.split(input_df, 'dataframes/')[1]
    prefix = str.split(df_name, '.csv')[0]
    plot_kkclade_heatmap_mut(prefix)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Determines epistasis between egg-specific mutations")
    parser.add_argument('--in_file', help= "input dataframe file")
    args = parser.parse_args()

    main(input_df = args.in_file)
