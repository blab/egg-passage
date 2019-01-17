"""
Find epistasis between egg-specific mutations
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def find_epistasis(prefix, threshold):
    df = pd.read_csv('data/'+prefix+'_egg_df.csv')

    mut_sites = [col for col in df.columns if col[0:3]=='mut']

    pos_interaction = []
    neg_interaction = []

    epistatic_interactions = []

    for site_1 in mut_sites:
        #find proportion of all egg seqs with mutation at this site
        prev_1 = float(len(df[df[site_1] == True]))/float(len(df))
        # egg_prev = {}
        # egg_prev['site']= site_1[3:]
        # egg_prev['prev']= prev_1
        # egg_prev['interaction']= 'overall'
        # epistatic_interactions.append(egg_prev)
        for site_2 in mut_sites:
            if site_1 != site_2:
                prev_2 = float(len(df[df[site_2] == True]))/float(len(df))
                epi_int = {}
                epi_int['site'] = site_2[3:]
                epi_int['prev'] = prev_2
                epi_int['interaction'] = 'overall'
                epistatic_interactions.append(epi_int)

                prev_2_in_1 = float(len(df[(df[site_1] == True) & (df[site_2] == True)]))/ float(len(df[df[site_1] == True]))
                prev_1_in_2 = float(len(df[(df[site_1] == True) & (df[site_2] == True)]))/ float(len(df[df[site_2] == True]))
                if prev_2_in_1 > (threshold*prev_2) or prev_2_in_1 < (prev_2/threshold):
                    epi_int = {}
                    epi_int['site'] = site_2[3:]
                    epi_int['prev'] = prev_2_in_1
                    epi_int['interaction'] = site_1[3:]
                    epistatic_interactions.append(epi_int)
                if prev_2_in_1 > (threshold*prev_2):
                    pos_interaction += [str(site_1)+' and '+ str(site_2)]
                elif prev_2_in_1 < (prev_2/threshold):
                    neg_interaction += [str(site_1)+' and '+ str(site_2)]

    plot_df = pd.DataFrame(epistatic_interactions).drop_duplicates()
    #Remove all sites without interactions
    for site in mut_sites:
        if len(plot_df[plot_df.site == site[3:]]) == 1:
            plot_df = plot_df[plot_df.site != site[3:]]

    cmap = {'160':'blue', '194':'green', '186':'orange', '203':'red', '225':'yellow', '156':'purple', '219':'black'}
    plot_df['color'] = plot_df['site'].map(cmap)

    sns.set(style="whitegrid")
    sites = list(np.unique(plot_df['site']))
    fig, ax = plt.subplots(1, len(sites), sharey=True)
    for site in sites:
        sns.barplot(x='site', y='prev', hue='interaction', palette=sns.color_palette("hls", 8), data=plot_df[plot_df.site == site], ax = ax[sites.index(site)])
    for axis in ax.flat:
        for p in axis.patches:
            h = p.get_height()
            x = p.get_x()+p.get_width()/2.
            axis.annotate("%g" % p.get_height(), xy=(x,h), xytext=(0,4), rotation=90, size=6,
                       textcoords="offset points", ha="center", va="bottom")
        axis.margins(0.2)
        axis.set(ylabel='prevalence of mutation')
        axis.xaxis.label.set_visible(False)
        axis.legend(loc= 9,prop={'size': 5})
        axis.label_outer()

    fig.text(0.5, 0.04, 'HA position of mutation', ha='center', va='center')
    fig.subplots_adjust(wspace=0.4)
    fig.savefig('plots/egg_mutation_epistasis_'+str(prefix)+'.pdf')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Determines epistasis between egg-specific mutations")
    parser.add_argument('--prefix', default= 'h3n2_6y_notiter', help= "specify prefix for naming data files")
    parser.add_argument('--threshold', default= 3.0, type= int, help= "threshold for fold-change in mutation prevalence in order to determine positive and negative epistatic interactions")
    args = parser.parse_args()

    find_epistasis(prefix = args.prefix, threshold = args.threshold)
