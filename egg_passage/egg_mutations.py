"""
Find egg-passaging-specific mutations
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency

def find_mutation_prevalence(prefix):
    """
    Find the overall proportion of egg-, cell- or un-passaged sequences that have a mutation at a given site. Test whether mutation prevalence is statistically different in egg-passaged sequences by chi-square test
    """

    df = pd.read_csv('data/'+prefix+'_df.csv')

    mut_sites = [col for col in df.columns if col[0:3]=='mut']

    #Find prevalence of mutation at each site in egg-passaged vs. other viruses (as proportion of egg-passaged sequences that have that mutation)
    mut_prev = {}
    chi_square = {}
    for mut_site in mut_sites:
        #make dict to store prev of each mutation
        mut_prev_site = {}
        #make dict to store observations for chi-squared test
        chi_obs = {}
        #make dict to store results of chi-squared test
        chi_pvalues = {}

        for pas_type in df['passage'].unique():
            mut_count = float(len(df[(df.passage == pas_type) & (df[mut_site] == True)]))/float(len(df[df.passage == pas_type]))
            mut_prev_site[pas_type] = mut_count

            mut_obs = len(df[(df.passage == pas_type) & (df[mut_site] == True)])
            nomut_obs = len(df[(df.passage == pas_type) & (df[mut_site] == False)])
            chi_obs[pas_type] = [mut_obs, nomut_obs]

        egg_v_cell_obs = np.array([chi_obs['egg'], chi_obs['cell']])
        ec_chi2, ec_p, ec_dof, ec_expected = chi2_contingency(egg_v_cell_obs)
        if ec_p > 0.000001:
            print('%s is not not significantly more mutated in egg-passaged vs. cell-passaged strains' %str(mut_site))
        egg_v_unpassaged_obs = np.array([chi_obs['egg'], chi_obs['0']])
        eu_chi2, eu_p, eu_dof, eu_expected = chi2_contingency(egg_v_unpassaged_obs)
        if ec_p > 0.000001:
            print('%s is not not significantly more mutated in egg-passaged vs. unpassaged strains' %str(mut_site))
        chi_pvalues['egg_v_cell'] = ec_p
        chi_pvalues['egg_v_unpassaged'] = eu_p

        mut_prev[mut_site] = mut_prev_site
        chi_square[mut_site] = chi_pvalues


    mut_prev_df = pd.DataFrame(mut_prev)
    chi_square_df = pd.DataFrame(chi_square)

    return mut_prev_df, chi_square_df


def plot_mutation_prev(prefix):
    """
    Plot the prevalence of mutations at identified HA sites in egg-passaged, cell-passaged and unpassaged viruses
    """

    mut_prev_df, chi_square_df = find_mutation_prevalence(prefix)

    #Make dictionary to store significance of chi-squared tests at each site
    site_sig = {}
    for site, stats in chi_square_df.items():
        stat_sig = 0 #statistical significance (0 if neither egg_v_cell or egg_v_unpassaged is significant, 2 if both are, 1 if only 1 is)
        for stat in stats:
            if stat < 0.00001:
                stat_sig+=1
        site_sig[site[3:]] = stat_sig

    plot_df = mut_prev_df.unstack().reset_index().rename(columns={'level_0':'site', 'level_1':'virus_passage', 0:'prevalence'})
    plot_df['site'] = plot_df['site'].str[3:]
    plot_df['virus_passage'] = np.where(plot_df['virus_passage']=='0', 'unpassaged', plot_df['virus_passage'])

    sns.set(style="white")
    barplot = sns.barplot(x='site', y='prevalence', hue='virus_passage', data=plot_df)
    barplot.set(xlabel='HA position', ylabel='prevalence of mutation')
    barplot.get_figure().savefig('plots/egg_mutation_prevalence_'+str(prefix)+'.pdf')

def plot_overall_aas(prefix):
    """
    Plot overall amino diversity at each site before egg-passaging. Indicate the proportion of each genotype that mutated during egg-passaging
    """

    df = pd.read_csv('data/'+prefix+'_egg_df.csv')

    mut_sites = [col for col in df.columns if col[0:3]=='mut']

    before_aas = []

    for mut_site in mut_sites:
        bottom_proportion = 0.0
        stack = 0

        unpassaged_group = df.groupby('circulating'+str(mut_site[3:])).size().reset_index(name='count')

        for k, v in unpassaged_group.iterrows():

            site_group = df[df['circulating'+str(mut_site[3:])] == v['circulating'+str(mut_site[3:])]].groupby(mut_site).size().reset_index(name='count')

            for i, r in site_group.iterrows():
                proportion = float(r['count'])/float(len(df))
                before_aas.append({'site': mut_site, 'aa': v['circulating'+str(mut_site[3:])], 'mutated': r['mut'+str(mut_site[3:])], 'proportion': proportion, 'bottom_proportion': bottom_proportion, 'stack': stack})
                bottom_proportion += proportion
                stack +=1

    before_aas_df = pd.DataFrame(before_aas)

    #Add x position for grouped bar chart
    pmap = {}
    x_pos = 2
    for pos in before_aas_df['site'].unique():
        pmap[pos] = x_pos
        x_pos += 2

    #cmap based off 'stack' (0/1 are first genotype, 2/3 are second)
    cmap = {0:'#c2efec', 1:'#33c7bc', 2:'#e5c1f0', 3:'#b045d3'}

    width = 1.2

    fig, ax = plt.subplots()

    plt.bar(before_aas_df['site'].map(pmap), before_aas_df['proportion'], width, bottom = before_aas_df['bottom_proportion'], color = before_aas_df['stack'].map(cmap))
    for i, r in before_aas_df[before_aas_df['mutated']==False].iterrows():
        plt.annotate(r['aa'], xy=((pmap[r['site']]), (r['bottom_proportion']+0.27*r['proportion'])), color="black", ha='center')


    ax.set(xlabel = 'HA site', ylabel='Genotype')
    ax.set_xticks([p for p in pmap.values()])
    ax.set_xticklabels([s[3:] for s in pmap.keys()])
    plt.text(18, 0.97, 'Mutated during \negg-passaging', color='black', size=8, fontweight='bold')
    plt.text(18, 0.85, 'False', color='#c2efec', size=8)
    plt.text(19.1, 0.85, 'True', color='#33c7bc', size=8)
    plt.text(18, 0.77, 'False', color='#e5c1f0', size=8)
    plt.text(19.1, 0.77, 'True', color='#b045d3', size=8)
    plt.text(6, 1.08, 'H3N2 strain genotypes', color='black', size=12)


    fig.savefig('plots/genotypes_'+str(prefix)+'.pdf', bbox_inches='tight')

def find_mutation_aas(prefix):
    """
    Determine which mutations happen at each site. For each site, plot strains that mutated during egg passaging, showing the prevalence of each genotype before and after egg-passaging
    """

    df = pd.read_csv('data/'+prefix+'_egg_df.csv')

    mut_sites = [col for col in df.columns if col[0:3]=='mut']

    mut_aas = []

    for mut_site in mut_sites:

        unpassaged_group = df[df[str(mut_site)] == True].groupby('circulating'+str(mut_site[3:])).size().reset_index(name='count')
        mutation_group = df[df[str(mut_site)] == True].groupby(str(mut_site[3:])).size().reset_index(name='count')

        for aa_count in range(len(unpassaged_group)):

            if aa_count == 0:
                mut_aas.append({'site': mut_site, 'aa': unpassaged_group['circulating'+str(mut_site[3:])][aa_count], 'proportion': float(unpassaged_group['count'][aa_count])/float(len(df[df[str(mut_site)] == True])), 'stack': aa_count, 'strain': 'unpassaged', 'bottom_proportion': 0.0})
            #Find height where stacked bar should start
            else:
                mut_aas.append({'site': mut_site, 'aa': unpassaged_group['circulating'+str(mut_site[3:])][aa_count], 'proportion': float(unpassaged_group['count'][aa_count])/float(len(df[df[str(mut_site)] == True])), 'stack': aa_count, 'strain': 'unpassaged', 'bottom_proportion': float(unpassaged_group['count'][aa_count-1])/float(len(df[df[str(mut_site)] == True]))})

        for aa_count in range(len(mutation_group)):

            if aa_count == 0:
                mut_aas.append({'site': mut_site, 'aa': mutation_group[str(mut_site[3:])][aa_count], 'proportion': float(mutation_group['count'][aa_count])/float(len(df[df[str(mut_site)] == True])), 'stack': aa_count, 'strain': 'egg', 'bottom_proportion': 0.0})
            #Find height where stacked bar should start
            elif aa_count == 1:
                mut_aas.append({'site': mut_site, 'aa': mutation_group[str(mut_site[3:])][aa_count], 'proportion': float(mutation_group['count'][aa_count])/float(len(df[df[str(mut_site)] == True])), 'stack': aa_count, 'strain': 'egg', 'bottom_proportion': float(mutation_group['count'][aa_count-1])/float(len(df[df[str(mut_site)] == True]))})
            elif aa_count == 2:
                mut_aas.append({'site': mut_site, 'aa': mutation_group[str(mut_site[3:])][aa_count], 'proportion': float(mutation_group['count'][aa_count])/float(len(df[df[str(mut_site)] == True])), 'stack': aa_count, 'strain': 'egg', 'bottom_proportion': float(mutation_group['count'][aa_count-1]+mutation_group['count'][aa_count-2])/float(len(df[df[str(mut_site)] == True]))})
            elif aa_count == 3:
                mut_aas.append({'site': mut_site, 'aa': mutation_group[str(mut_site[3:])][aa_count], 'proportion': float(mutation_group['count'][aa_count])/float(len(df[df[str(mut_site)] == True])), 'stack': aa_count, 'strain': 'egg', 'bottom_proportion': float(mutation_group['count'][aa_count-1]+mutation_group['count'][aa_count-2]+mutation_group['count'][aa_count-3])/float(len(df[df[str(mut_site)] == True]))})
            elif aa_count == 4:
                mut_aas.append({'site': mut_site, 'aa': mutation_group[str(mut_site[3:])][aa_count], 'proportion': float(mutation_group['count'][aa_count])/float(len(df[df[str(mut_site)] == True])), 'stack': aa_count, 'strain': 'egg', 'bottom_proportion': float(mutation_group['count'][aa_count-1]+mutation_group['count'][aa_count-2]+mutation_group['count'][aa_count-3]+mutation_group['count'][aa_count-4])/float(len(df[df[str(mut_site)] == True]))})

    mut_aas_df = pd.DataFrame(mut_aas)

    #Add fake entries so that all stacks have the same number of rows, so seaborn doesn't complain
    #Not needed with matplotlib
    # for stack in range(len(mut_aas_df['stack'].unique())):
    #     for mut_site in mut_sites:
    #         for strain in mut_aas_df['strain'].unique():
    #             if len(mut_aas_df[(mut_aas_df['site'] == str(mut_site)) & (mut_aas_df['strain'] == str(strain)) & (mut_aas_df['stack'] == int(stack))]) == 0:
    #                 mut_aas_df = mut_aas_df.append({'site': mut_site, 'aa': '', 'proportion': 0.0, 'bottom_proportion': 0.0, 'stack':stack, 'strain': strain}, ignore_index=True)

    #Add x position for grouped bar chart
    pmap = {}
    x_pos = 2
    for pos in mut_aas_df['site'].unique():
        pmap[pos] = x_pos
        x_pos += 2

    cmap_u = {0:'#33c7bc', 1:'#159794', 2:'#008080', 3:'#24aea7', 4:'#40e0d0'}
    cmap_e = {0:'#cb2f44', 1:'#ffbd84', 2:'#8b0000', 3:'#f47461' , 4:'#FFE4B5'}
    width = 0.6

    fig, ax = plt.subplots()

    for strain in mut_aas_df['strain'].unique():
        for stack in range(len(mut_aas_df['stack'].unique())):
            stack_df = mut_aas_df[(mut_aas_df['stack']==int(stack)) & (mut_aas_df['strain']==str(strain))]
            if strain=='unpassaged':
                plt.bar(stack_df['site'].map(pmap), stack_df['proportion'], width, bottom = stack_df['bottom_proportion'], color = stack_df['stack'].map(cmap_u))
                for i, r in stack_df.iterrows():
                    plt.annotate(r['aa'], xy=((pmap[r['site']]), (r['bottom_proportion']+0.3*r['proportion'])), color="black", ha='center')
            else:
                plt.bar((stack_df['site'].map(pmap)+width+0.05), stack_df['proportion'], width, bottom = stack_df['bottom_proportion'], color = stack_df['stack'].map(cmap_e))
                for i, r in stack_df.iterrows():
                    plt.annotate(r['aa'], xy=((pmap[r['site']]+width+0.5), (r['bottom_proportion']+0.3*r['proportion'])), color="black", size=6, ha='center')


    ax.set(xlabel = 'HA site', ylabel='Genotype')
    ax.set_xticks([p for p in pmap.values()])
    ax.set_xticklabels([s[3:] for s in pmap.keys()])
    plt.text(18, 0.97, 'BEFORE \n egg-passaging', color='#33c7bc', size=8, fontweight='bold')
    plt.text(18, 0.85, 'AFTER \n egg-passaging', color='#cb2f44', size=8, fontweight='bold')
    plt.text(3, 1.08, 'H3N2 strains that mutated during egg-passaging', color='black', size=12)


    fig.savefig('plots/before_after_eggpassaging_'+str(prefix)+'.pdf', bbox_inches='tight')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Determines egg-specific mutations")
    parser.add_argument('--prefix', default= 'h3n2_6y_notiter', help= "specify prefix for naming data files")
    args = parser.parse_args()

    plot_mutation_prev(prefix = args.prefix)
    find_mutation_aas(prefix = args.prefix)
    plot_overall_aas(prefix = args.prefix)
