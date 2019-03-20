"""
Find background genotype specifity of egg-passaging mutations
"""

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

def plot_background_heatmap(prefix):
    """
    Plot heatmaps indicating the genotypic backgrounds each egg-passaging mutation prefentially occurs in. Heat maps indicate the log enrichment ratio, which compares the frequency of each egg-passaging mutation in each clade to its overall frequency
    """

    df = pd.read_csv('dataframes/'+prefix+'_egg.csv')

    sites = [col[3:] for col in df.columns if col[0:3]=='mut']
    sites = sorted(sites)

    #Remove clades with less than 5 seqs from analysis
    clade_size = df.groupby('clade').size().reset_index().rename(columns={0:'clade_total'})
    for k, v in clade_size.iterrows():
        if v['clade_total']<=5:
            index_names = df[df['clade']==v['clade']].index
            df.drop(index_names, inplace=True)

    #Find enrichment of egg-passaged mutation at each site in each clade
    clade_enrichment = []
    for site in sites:
        overall_pct = len(df[df['mut'+str(site)]==True])/len(df)
        clade_count = df[df['mut'+str(site)]==True].groupby('clade').size()
        clade_pct = df[df['mut'+str(site)]==True].groupby('clade').size()/df.groupby('clade').size()
        clade_count = clade_count.reset_index().rename(columns={0:'count'})
        #Fill 0-count clades with small number, for log
        clade_pct = clade_pct.reset_index().rename(columns={0:'pct'}).fillna(0.000001)
        site_clade = pd.merge(clade_count, clade_pct, how='right', on='clade').fillna(0)
        site_clade['enrichment'] = site_clade['pct']/overall_pct
        site_clade['log_enrichment'] = np.log2(site_clade['enrichment'])
        #Replace small number with 0.0 for 0-count clades
        site_clade['pct'] = np.where(site_clade['pct']==0.000001, 0.0, site_clade['pct'])
        site_clade['site'] = (site+'\n ( '+str(round(overall_pct, 2))+' )')
        clade_enrichment.append(pd.DataFrame(site_clade))

    clade_enrichment = pd.concat(clade_enrichment, axis=0)
    #Reorder clades to reflect tree topology
    clade_enrichment.clade = pd.CategoricalIndex(clade_enrichment.clade,
                                                 categories= ["unassigned", "3b", "3c", "3c3",
                                                 "3c3.A", "3c3.B", "3c2", "3c2.A", "A4", "A3",
                                                 "A2", "A2/re", "A1", "A1a", "A1b", "A1b/135N",
                                                 "A1b/135K", "A1b/131K"])
    clade_enrichment.sort_values(by='clade', inplace=True)

    clade_heatmap = clade_enrichment.pivot('clade', 'site', 'log_enrichment')
    heatmap_annotate = clade_enrichment.pivot('clade', 'site', 'pct')

    fig, ax= plt.subplots()
    sns.set(style="white")

    fig = sns.heatmap(clade_heatmap,
                      cmap=sns.diverging_palette(220, 20, sep=80, as_cmap=True),
                      vmax=3.25, vmin=-3.25, center = 0.0,
                      annot=heatmap_annotate, cbar_kws={'label': 'log2 enrichment'},
                      annot_kws={'fontsize':6})
    fig.set_xlabel('HA1 position \n (overall frequency in egg-passaged viruses)',fontsize=14)
    fig.set_ylabel('Clade',fontsize=14)
    fig.get_figure().savefig('plots/'+str(prefix)+'/genetic_background_heatmap_'+str(prefix)+'.pdf', bbox_inches='tight')

def plot_kkclade_heatmap(prefix):

    df = pd.read_csv('dataframes/'+prefix+'_egg.csv')

    sites = [col[3:] for col in df.columns if col[0:3]=='mut']
    sites = sorted(sites)
    #Remove clades with less than 5 egg seqs from analysis
    clade_size = df.groupby('kk_clade').size().reset_index().rename(columns={0:'kk_clade_total'})
    for k, v in clade_size.iterrows():
        if v['kk_clade_total']<=5:
            index_names = df[df['kk_clade']==v['kk_clade']].index
            df.drop(index_names, inplace=True)

    #Find enrichment of egg-passaged mutation at each site in each clade
    clade_enrichment = []
    for site in sites:
        overall_pct = len(df[df['mut'+str(site)]==True])/len(df)
        clade_count = df[df['mut'+str(site)]==True].groupby('kk_clade').size()
        clade_pct = df[df['mut'+str(site)]==True].groupby('kk_clade').size()/df.groupby('kk_clade').size()
        clade_count = clade_count.reset_index().rename(columns={0:'count'})
        #Fill 0-count clades with small number, for log
        clade_pct = clade_pct.reset_index().rename(columns={0:'pct'}).fillna(0.000001)
        site_clade = pd.merge(clade_count, clade_pct, how='right', on='kk_clade').fillna(0)
        site_clade['enrichment'] = site_clade['pct']/overall_pct
        site_clade['log_enrichment'] = np.log2(site_clade['enrichment'])
        #Replace small number with 0.0 for 0-count clades
        site_clade['pct'] = np.where(site_clade['pct']==0.000001, 0.0, site_clade['pct'])
        site_clade['site'] = (site+'\n ( '+str(round(overall_pct, 2))+' )')
        clade_enrichment.append(pd.DataFrame(site_clade))

    clade_enrichment = pd.concat(clade_enrichment, axis=0)
    #Reorder clades
    clade_enrichment.kk_clade = pd.CategoricalIndex(clade_enrichment.kk_clade,
                                                 categories= [v['kk_clade'] for k,v in (df.groupby('kk_clade')['date'].mean().reset_index().sort_values(by='date')).iterrows()])
    clade_enrichment.sort_values(by='kk_clade', inplace=True)

    clade_heatmap = clade_enrichment.pivot('kk_clade', 'site', 'log_enrichment')
    heatmap_annotate = clade_enrichment.pivot('kk_clade', 'site', 'pct')

    fig, ax= plt.subplots(figsize=(10,20))
    sns.set(style="white")

    fig = sns.heatmap(clade_heatmap,
                      cmap=sns.diverging_palette(220, 20, sep=80, as_cmap=True),
                      vmax=3.25, vmin=-3.25, center = 0.0,
                      annot=heatmap_annotate, cbar_kws={'label': 'log2 enrichment'},
                      annot_kws={'fontsize':6})
    fig.set_xlabel('HA1 position \n (overall frequency in egg-passaged viruses)',fontsize=14)
    fig.set_ylabel('Clade',fontsize=14)
    fig.get_figure().savefig('plots/'+str(prefix)+'/genetic_background_kkclade_heatmap_'+str(prefix)+'.pdf', bbox_inches='tight')

def main(input_df, clades_file):
    df_name = str.split(input_df, 'dataframes/')[1]
    prefix = str.split(df_name, '.csv')[0]
    if '6y' in prefix:
        plot_background_heatmap(prefix)
    plot_kkclade_heatmap(prefix)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Determines epistasis between egg-specific mutations")
    parser.add_argument('--in_file', help= "input dataframe file")
    parser.add_argument('--clades', help= "tsv containing clade defining mutations")
    args = parser.parse_args()

    main(input_df = args.in_file, clades_file = args.clades)
