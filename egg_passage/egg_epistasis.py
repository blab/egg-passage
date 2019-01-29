"""
Find epistasis between egg-specific mutations
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def find_epistasis(prefix, threshold):
    """
    Identify sites that mutate more or less often when another site is also mutated.
    """
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

    sns.set(style="white")
    sites = list(np.unique(plot_df['site']))

    fig, ax = plt.subplots(1, len(sites), sharey=True)
    for site in sites:
        sns.barplot(x='site', y='prev', hue='interaction', palette=sns.color_palette("GnBu_d", 5), data=plot_df[plot_df.site == site], ax = ax[sites.index(site)])
    for axis in ax.flat:
        for p in axis.patches:
            h = p.get_height()
            x = p.get_x()+p.get_width()/2.
            axis.annotate("%g" % p.get_height(), xy=(x,h), xytext=(0,4), rotation=90, size=6,
                       textcoords="offset points", ha="center", va="bottom")
        axis.margins(0.2)
        axis.xaxis.label.set_visible(False)
        axis.yaxis.label.set_visible(False)
        axis.legend(loc= 9,prop={'size': 5})

    for axis in range(len(ax)):
        ax[axis].spines['right'].set_visible(False)
        ax[axis].spines['top'].set_visible(False)
        if axis == 0:
            ax[axis].set(ylabel='Prevalence of mutation')
        else:
            ax[axis].spines['left'].set_visible(False)

    fig.text(0.5, 0.04, 'HA position of mutation', ha='center', va='center')
    fig.subplots_adjust(wspace=0.4)
    fig.savefig('plots/egg_mutation_epistasis_'+str(prefix)+'.pdf')

def find_aa_epistasis(prefix):
    """
    Calculate ratio between titers of egg-passaged seqs with mutations and average titer of strains in that clade
    """

    df = pd.read_csv('data/'+prefix+'_df.csv')

    mut_sites = [col for col in df.columns if col[0:3]=='mut']
    sites = [str(x[3:]) for x in mut_sites]


    aa_epistasis = []

    for mut_site1 in mut_sites:
        site1 = str(mut_site1[3:])
        for mut_site2 in mut_sites:
            site2 = str(mut_site2[3:])
            if site1 != site2:
                aa_group = df.groupby([site1, site2]).size().reset_index().rename(columns={0:'count'})
                bottom_count = 0
                bottom_proportion = 0.0
                stack = 0
                site1_stack = 'X'
                for k, v in aa_group.iterrows():

                    if str(v[site1]) == site1_stack:
                        proportion = float(v['count'])/float(len(df[df[site1]==str(v[site1])]))
                        aa_epistasis.append({'site1': site1, 'site1_aa': v[site1], 'site2': site2, 'site2_aa': v[site2], 'count': v['count'], 'bottom_count': bottom_count, 'stack': stack, 'proportion': proportion, 'bottom_proportion': bottom_proportion})
                        bottom_count += v['count']
                        bottom_proportion += proportion
                        stack +=1

                    else:
                        bottom_count = 0
                        bottom_proportion = 0.0
                        stack = 0
                        proportion = float(v['count'])/float(len(df[df[site1]==str(v[site1])]))
                        aa_epistasis.append({'site1': site1, 'site1_aa': v[site1], 'site2': site2, 'site2_aa': v[site2], 'count': v['count'], 'bottom_count': bottom_count, 'stack': stack, 'proportion': proportion, 'bottom_proportion': bottom_proportion})
                        site1_stack = str(v[site1])
                        bottom_count += v['count']
                        bottom_proportion += proportion
                        stack +=1

    aa_epistasis_df = pd.DataFrame(aa_epistasis)

    #Dict of most prevalent egg mutations, for outlining bars
    egg_muts = {'160':['K'], '194': ['P'], '186':['V'], '225':['G'], '219':['F','Y'], '203':['I'], '156':['R','Q'], '138':['S']}

    for site in sites:
        site_df = aa_epistasis_df[aa_epistasis_df['site1']==site]

        other_sites = [x for x in sites if site not in x]

        fig, ax = plt.subplots(1, len(other_sites), sharey=True)

        width = 1.4

        for site2 in other_sites:
            #Add x position for grouped bar chart
            pmap = {}
            x_pos = 2
            for pos in site_df['site1_aa'].unique():
                pmap[pos] = x_pos
                x_pos += 2

            plot_df = site_df[site_df['site2']==site2]

            site2aas = list(plot_df['site2_aa'].unique())
            colors = sns.color_palette('GnBu', len(site2aas))
            cmap = dict(zip(site2aas, colors))

            ax[other_sites.index(site2)].bar(plot_df['site1_aa'].map(pmap), plot_df['proportion'], width, bottom=plot_df['bottom_proportion'], color= plot_df['site2_aa'].map(cmap))
            for i, r in plot_df.iterrows():
                ax[other_sites.index(site2)].annotate(str(len(df[df[site]==r['site1_aa']])), xy=((pmap[r['site1_aa']]), 1.05), color="black", va='center', ha='center', rotation=90, size=6)
                if r['site2_aa'] in egg_muts[site2]:
                    label_color = 'red'
                else:
                    label_color = 'black'
                if r['proportion']>= 0.1:
                    ax[other_sites.index(site2)].annotate(r['site2_aa'], xy=((pmap[r['site1_aa']]), (r['bottom_proportion']+0.5*r['proportion'])), color= label_color, va='center', ha='center')
                else:
                    ax[other_sites.index(site2)].annotate(r['site2_aa'], xy=((pmap[r['site1_aa']]), (r['bottom_proportion']+0.5*r['proportion'])), color= label_color, va='center', ha='center', size=4)
            ax[other_sites.index(site2)].set(xlabel = str(site2))
            ax[other_sites.index(site2)].set_xticks([p for p in pmap.values()])
            ax[other_sites.index(site2)].set_xticklabels([s for s in pmap.keys()])
            tick_colors = ['red' if s in egg_muts[site] else 'black' for s in pmap.keys()]
            [t.set_color(i) for (i,t) in zip(tick_colors, ax[other_sites.index(site2)].xaxis.get_ticklabels())]

            for axis in range(len(ax)):
                ax[axis].spines['right'].set_visible(False)
                ax[axis].spines['top'].set_visible(False)
                if axis == 0:
                    ax[axis].set(ylabel='Proportion of H3N2 sequences')
                else:
                    ax[axis].spines['left'].set_visible(False)
                    ax[axis].tick_params(axis='y', length=0)

            plt.subplots_adjust(wspace=0.4)
            fig.text(0.075, 0.075, str(site)+':', ha='center', va='center', size=8)
            fig.text(0.92, 0.8, 'Prevalent \negg-mutations', color='red', ha='left', va='center', size=8)
            fig.text(0.075, 0.9, 'Number of \nstrains', ha='center', va='center', size=8)
            fig.suptitle('Epistasis between %s and other HA sites' %str(site), fontsize=12)
            fig.savefig('plots/epistasis/%s_epistasis_'%str(site)+str(prefix)+'.pdf', bbox_inches='tight')

def plot_heatmap(prefix):
    """
    Plot heatmaps indicating epistatic interactions between genotypes. Heat maps indicate whether the genotype of certain HA sites are correlated via a log enrichment ratio, that compares the frequency of observing genotype 1 at site 1 AND genotype 2 at site 2 versus the expected frequency (based on overall prevalence of the genotypes)
    """

    df = pd.read_csv('data/'+prefix+'_egg_df.csv')

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
    aa_epistasis_df = aa_epistasis_df[aa_epistasis_df['count'] >= 10]


    fig, ax = plt.subplots(len(sites), len(sites), sharex='col', sharey='row')

    #Make list of lists to store arrays of enrichment ratios, like matrix
    #Row1, column1 stores interactions with sites[1], etc
    plot_arrays = [[0 for x in sites] for y in sites]

    #Add matrix with number of amino acid genotypes as dimensions
    for site_x in sites:
        for site_y in sites:
            plot_arrays[sites.index(site_x)][sites.index(site_y)] = np.zeros(shape = (len(list(aa_epistasis_df[aa_epistasis_df['site1']==site_x]['site1_aa'].unique())), len(list(aa_epistasis_df[aa_epistasis_df['site2']==site_y]['site2_aa'].unique()))))


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
    egg_muts = {'160':'K', '194': 'P', '186':'V', '225':'G', '219':['F','Y'], '203':'I', '156':['R','Q'], '138':'S'}

    for site1 in range(len(sites)):
        for site2 in range(len(sites)):
            site1_list = list(aa_epistasis_df[aa_epistasis_df['site1']==str(sites[site1])]['site1_aa'].unique())
            site2_list = list(aa_epistasis_df[aa_epistasis_df['site2']==str(sites[site2])]['site2_aa'].unique())
            cmap = sns.diverging_palette(220, 20, sep=10, as_cmap=True)
            heatmap = ax[site1, site2].imshow(plot_arrays[site1][site2], cmap=cmap, vmin= -2.0, vmax=2.0, aspect='auto')
            # heatmap = ax[site1, site2].imshow(plot_arrays[site1][site2], vmin= -2.0, vmax=2.0, cmap='RdBu', aspect='auto')

            ax[site1, site2].tick_params(axis='both', which='both', length=0)
            # ax[site1, site2].tick_params(axis='x', length=0)
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
    fig.savefig('plots/epistasis_heatmap.pdf', bbox_inches='tight')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Determines epistasis between egg-specific mutations")
    parser.add_argument('--prefix', default= 'h3n2_6y_notiter', help= "specify prefix for naming data files")
    parser.add_argument('--threshold', default= 3.0, type= int, help= "threshold for fold-change in mutation prevalence in order to determine positive and negative epistatic interactions")
    args = parser.parse_args()

    find_epistasis(prefix = args.prefix, threshold = args.threshold)
    find_aa_epistasis(prefix = args.prefix)
    plot_heatmap(prefix = args.prefix)
