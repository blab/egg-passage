"""
Find phenotypic correlations with egg-specific mutations
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def find_aa_titers(prefix):
    """
    Identify correlations between egg mutations and titer drops
    """

    df = pd.read_csv('data/'+prefix+'_df.csv')
    mut_sites = [col for col in df.columns if col[0:3]=='mut']

    aa_titers = []

    #Manually enter predominant egg mutation genotypes
    egg_muts = {'160':'K', '194': 'P', '186':'V', '225':'G', '219':['F','Y'], '203':'I', '156':['R','Q'], '138':'S'}

    for mut_site in mut_sites:
        site = str(mut_site[3:])
        titer_mean = df.groupby(site)['cTiterSub'].agg(['mean', 'size', 'std']).reset_index()

        for k, v in titer_mean.iterrows():
            if v[site] in egg_muts[site]:
                aa_titers.append({'site': site, 'aa': v[site], 'mean': v['mean'], 'count': v['size'], 'std': v['std'], 'egg_mut':'red'})
            else:
                aa_titers.append({'site': site, 'aa': v[site], 'mean': v['mean'], 'count': v['size'], 'std': v['std'], 'egg_mut':'white'})

    aa_titers_df = pd.DataFrame(aa_titers)

    sns.set(style="white")
    sites = list(np.unique(aa_titers_df['site']))
    palette=sns.color_palette("GnBu_d", len(sites))

    fig, ax = plt.subplots(1, len(sites), sharey=True)
    for site in sites:
        site_df = aa_titers_df[aa_titers_df.site == site]
        # sns.barplot(x='aa', y='mean', color= palette[sites.index(site)], data=site_df, ax = ax[sites.index(site)])
        ax[sites.index(site)].bar(site_df['aa'], site_df['mean'], color= palette[sites.index(site)], yerr= site_df['std'], edgecolor= site_df['egg_mut'], linewidth= 1.5, error_kw={'elinewidth':1})
        ax[sites.index(site)].set(xlabel = str(site))

    for axis in range(len(ax)):
        ax[axis].spines['right'].set_visible(False)
        ax[axis].spines['top'].set_visible(False)
        if axis == 0:
            ax[axis].set(ylabel='average cTiterSub')
        else:
            ax[axis].spines['left'].set_visible(False)

    fig.text(0.5, 0, 'HA position', ha='center', va='center')
    fig.text(0.5, 0.9, 'Titer drop for strains with given genotypes', ha='center', va='center')
    fig.text(0.85, 0.8, 'Prevalent \negg-mutations', color='red', size=8)
    fig.savefig('plots/titers_'+str(prefix)+'.pdf', bbox_inches='tight')


def find_titer_ratio(prefix):
    """
    Calculate ratio between titers of egg-passaged seqs with mutations and average titer of strains in that clade
    """

    df = pd.read_csv('data/'+prefix+'_df.csv')
    egg_df = df[df['passage']=='egg']

    clade_titer_df = df.groupby('clade')['cTiterSub'].agg(['mean', 'std', 'size']).reset_index()

    titer_ratios = []

    mut_sites = [col for col in egg_df.columns if col[0:3]=='mut']
    for k,v in egg_df.iterrows():
        for mut_site in mut_sites:
            if v[mut_site] == True:
                clade_mean = clade_titer_df[clade_titer_df['clade']==v['clade']]['mean']
                titer_ratios.append({'mut_site':str(mut_site[3:]), 'clade':v['clade'], 'titer_ratio': float(v['cTiterSub'])/float(clade_mean)})

    titer_ratios_df = pd.DataFrame(titer_ratios)

    means = titer_ratios_df.groupby('mut_site')['titer_ratio'].mean().values
    mean_labels = [str(np.round(x, 2)) for x in means.tolist()]
    counts = titer_ratios_df.groupby('mut_site')['titer_ratio'].size().values
    count_labels = ['n= '+str(x) for x in counts.tolist()]

    fig, ax = plt.subplots()
    sns.boxplot(x = 'mut_site', y = 'titer_ratio', data = titer_ratios_df, color='#33c7bc')
    ax.set(ylabel='Titer ratio (egg-passaging mutation/ clade average)')
    ax.set(xlabel='HA position')

    pos = range(len(means))
    for tick,label in zip(pos,ax.get_xticklabels()):
        ax.text(pos[tick], means[tick] + 0.55, mean_labels[tick]+'\n'+count_labels[tick], horizontalalignment='center', size='small', color='black')

    fig.savefig('plots/titer_ratios_'+str(prefix)+'.pdf')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Determines egg-specific mutations")
    parser.add_argument('--prefix', default= 'h3n2_6y_fra', help= "specify prefix for naming data files")
    args = parser.parse_args()

    find_aa_titers(prefix = args.prefix)
    find_titer_ratio(prefix = args.prefix)
