"""
Find epistasis between egg-specific mutations
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_heatmap(prefix):
    """
    Plot heatmaps indicating epistatic interactions between genotypes. Heat maps indicate whether the genotype of certain HA sites are correlated via a log enrichment ratio, that compares the frequency of observing genotype 1 at site 1 AND genotype 2 at site 2 versus the expected frequency (based on overall prevalence of the genotypes)
    """

    df = pd.read_csv('dataframes/'+prefix+'_egg.csv')

    mut_sites = [col for col in df.columns if col[0:3]=='mut']
    sites = [str(x[3:]) for x in mut_sites]


    aa_epistasis = []
    aa_alone = []

    for site1 in sites:
        other_sites = [x for x in sites if site1 not in x]
        #Group egg sequences by their genotype at site1
        site1_counts = df.groupby(site1).size().reset_index().rename(columns={0:'count'})

        for k, v in site1_counts.iterrows():
            #Find the prevalence of each genotype at site1
            overall_proportion = float(v['count'])/float(len(df))
            aa_alone.append({'site': site1, 'aa': v[site1], 'overall_proportion': overall_proportion})

        for site2 in other_sites:
            #Group egg sequences by site1 and all other sites
            both_counts = df.groupby([site1, site2]).size().reset_index().rename(columns={0:'count'})

            for i, r in both_counts.iterrows():
                #Find prevalence of each site1 genotype co-occuring with each site2 genotype
                proportion = float(r['count'])/float(len(df))
                aa_epistasis.append({'site1': site1, 'site1_aa': r[site1], 'site2': site2, 'site2_aa': r[site2], 'count': r['count'], 'proportion': proportion})


    aa_epistasis_df = pd.DataFrame(aa_epistasis)
    aa_alone_df = pd.DataFrame(aa_alone)

    #set threshold to exclude groups with only a few sequences
    #note: this means proportions won't exactly add to 1- doesn't matter for comparison to exp.
    aa_epistasis_df = aa_epistasis_df[aa_epistasis_df['count'] >= 5]


    fig, ax = plt.subplots(len(sites), len(sites), sharex='col', sharey='row')

    #Make list of lists to store arrays of enrichment ratios, like matrix
    #Row1, column1 stores interactions with sites[1], etc
    plot_arrays = [[0 for x in sites] for y in sites]

    #Add matrix with number of amino acid genotypes as dimensions
    for site_x in sites:
        for site_y in sites:
            # plot_arrays[sites.index(site_x)][sites.index(site_y)] = np.zeros(shape = (len(list(aa_epistasis_df[aa_epistasis_df['site1']==site_x]['site1_aa'].unique())), len(list(aa_epistasis_df[aa_epistasis_df['site2']==site_y]['site2_aa'].unique()))))
            plot_arrays[sites.index(site_x)][sites.index(site_y)] = np.empty(shape = (len(list(aa_epistasis_df[aa_epistasis_df['site1']==site_x]['site1_aa'].unique())), len(list(aa_epistasis_df[aa_epistasis_df['site2']==site_y]['site2_aa'].unique()))))
            #Fill all with -2.0, so if observed count is 0, log enrichment will be -2.0
            plot_arrays[sites.index(site_x)][sites.index(site_y)].fill(-2.0)


    #Fill in plot_arrays with enrichment ratios
    for k, v in aa_epistasis_df.iterrows():
        obs = v['proportion']
        exp = float(aa_alone_df[(aa_alone_df['site'] == v['site1']) & (aa_alone_df['aa'] == v['site1_aa'])]['overall_proportion']) * float(aa_alone_df[(aa_alone_df['site'] == v['site2']) & (aa_alone_df['aa'] == v['site2_aa'])]['overall_proportion'])
        enrichment = float(obs)/exp
        log_enrichment = np.log2(enrichment)
        site1_list = list(aa_epistasis_df[aa_epistasis_df['site1']==v['site1']]['site1_aa'].unique())
        site2_list = list(aa_epistasis_df[aa_epistasis_df['site2']==v['site2']]['site2_aa'].unique())
        plot_arrays[sites.index(v['site1'])][sites.index(v['site2'])][site1_list.index(v['site1_aa'])][site2_list.index(v['site2_aa'])]= log_enrichment

    #Manually enter predominant egg mutation genotypes
    egg_muts = {'160':'K', '194': 'P', '186':'V', '225':'G', '219':['F','Y'], '203':'I', '156':['R','Q'], '138':'S', '246':['H'], '183':['L']}

    for site1 in range(len(sites)):
        for site2 in range(len(sites)):
            site1_list = list(aa_epistasis_df[aa_epistasis_df['site1']==str(sites[site1])]['site1_aa'].unique())
            site2_list = list(aa_epistasis_df[aa_epistasis_df['site2']==str(sites[site2])]['site2_aa'].unique())
            cmap = sns.diverging_palette(220, 20, sep=10, as_cmap=True)
            heatmap = ax[site1, site2].imshow(plot_arrays[site1][site2], cmap=cmap, vmin= -2.0, vmax=2.0, aspect='auto')

            ax[site1, site2].tick_params(axis='both', which='both', length=0)
            ax[site1, site2].set_xticks([p for p in range(len(site2_list))])
            ax[site1, site2].set_yticks([p for p in range(len(site1_list))])
            ax[site1, site2].tick_params(labelbottom=False)
            ax[site1, site2].set_xticks(np.arange(len(site2_list)+1)-.5, minor=True)
            ax[site1, site2].set_yticks(np.arange(len(site1_list)+1)-.5, minor=True)
            ax[site1, site2].grid(which='minor', color='white', linestyle='-', linewidth=1)

            if site1 > (site2-1):
                ax[site1, site2].set_visible(False)

            if site1 == 0:
                ax[site1, site2].xaxis.set_label_position('top')
                ax[site1, site2].xaxis.set_ticks_position('top')
                ax[site1, site2].set(xlabel = str(sites[site2]))
                ax[site1, site2].set_xticklabels([str(s) for s in site2_list])
                colors= ['red' if str(s) in egg_muts[str(sites[site2])] else 'black' for s in site2_list]
                for xtick, color in zip(ax[site1, site2].get_xticklabels(), colors):
                    xtick.set_color(color)

            if site2 == (site1+1):
                ax[site1, site2].set(ylabel = str(sites[site1]))
                ax[site1, site2].yaxis.set_ticks_position('left')
                ax[site1, site2].set_yticklabels([str(s) for s in site1_list])
                colors= ['red' if str(s) in egg_muts[str(sites[site1])] else 'black' for s in site1_list]
                for ytick, color in zip(ax[site1, site2].get_yticklabels(), colors):
                    ytick.set_color(color)

    cbar_ax = fig.add_axes([0.95, 0.2, 0.05, 0.7])
    colorbar = fig.colorbar(heatmap, cax=cbar_ax)
    colorbar.ax.get_yaxis().labelpad = 15
    colorbar.ax.set_ylabel('log2 enrichment ratio', rotation=270)
    fig.suptitle('Epistasis between HA sites in egg-passaged influenza H3N2', fontsize=12, y=1.05, x=0.6)
    fig.savefig('plots/'+str(prefix)+'/epistasis_heatmap_'+str(prefix)+'.pdf', bbox_inches='tight')

def main(input_df):
    df_name = str.split(input_df, 'dataframes/')[1]
    prefix = str.split(df_name, '.csv')[0]
    plot_heatmap(prefix)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Determines epistasis between egg-specific mutations")
    parser.add_argument('--in_file', help= "input dataframe file")
    args = parser.parse_args()

    main(input_df = args.in_file)
