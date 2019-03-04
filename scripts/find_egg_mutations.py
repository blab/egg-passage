"""
Find egg-passaging-specific mutations
"""

import argparse
import ast
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency
from Bio import SeqIO
from collections import Counter

def find_mutation_prevalence(prefix):
    """
    Find the overall proportion of egg-, cell- or un-passaged sequences that have a mutation at a given site. Test whether mutation prevalence is statistically different in egg-passaged sequences by chi-square test
    """

    df = pd.read_csv('dataframes/'+prefix+'.csv')

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
        egg_v_unpassaged_obs = np.array([chi_obs['egg'], chi_obs['unpassaged']])
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


def plot_mutation_site(prefix):
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

    sns.set(style="white")
    fig, ax = plt.subplots()
    passage_palette = {'unpassaged': '#5c3d46', 'cell': '#f8c968', 'egg': '#99bfaa'}
    fig = sns.barplot(x='site', y='prevalence', hue='virus_passage', data=plot_df, palette= passage_palette)
    fig.set(xlabel='HA position', ylabel='prevalence of mutation')
    fig.get_figure().savefig('plots/'+str(prefix)+'/egg_mutation_site_prevalence_'+str(prefix)+'.pdf')

    # return plot_df

def plot_mutation_aa(prefix):
    """
    Plot the prevalence of the specific amino acid mutations identified as enriched in egg-passaged sequences. Plot the prevelence in egg-passaged, cell-passaged and unpassaged viruses
    """
    tip_df = pd.read_csv('dataframes/'+prefix+'_tidy.csv')
    df = pd.read_csv('dataframes/'+prefix+'.csv')

    top_muts = {}
    for pas_type in tip_df['passage'].unique():
        top = (tip_df[tip_df.passage==pas_type].groupby('mutation')['mutation']
                ).count().sort_values(ascending=False)[:10]
        top_muts[pas_type] = list((g_name, g) for g_name, g in top.iteritems())
    egg_top_muts = [x[0] for x in top_muts['egg']]

    top_aa_muts = {}
    for egg_top_mut in egg_top_muts:
        site = egg_top_mut.split('HA1')[1][1:-1]
        from_aa = egg_top_mut.split('HA1')[1][0]
        to_aa = egg_top_mut.split('HA1')[1][-1]
        aa_mut_count = {}

        for pas_type in df['passage'].unique():
            mut_count = float(len(df[(df[str(site)+'_lastnode'] == from_aa) & (df[str(site)] == to_aa) & (df['passage'] == pas_type)])) / float(len(df[df.passage == pas_type]))
            aa_mut_count[pas_type] = mut_count

        top_aa_muts[egg_top_mut] = aa_mut_count

    aa_muts_df = pd.DataFrame(top_aa_muts)
    plot_aa_muts_df = aa_muts_df.unstack().reset_index().rename(columns={'level_0':'mutation', 'level_1':'virus_passage', 0:'prevalence'})

    aa_barplot, ax = plt.subplots()
    sns.set(style="white")
    passage_palette = {'unpassaged': '#5c3d46', 'cell': '#f8c968', 'egg': '#99bfaa'}
    aa_barplot = sns.barplot(x= 'mutation', y= 'prevalence', hue= 'virus_passage', data= plot_aa_muts_df, palette = passage_palette)
    aa_barplot.set(xlabel='HA1 mutation', ylabel='prevalence of mutation')
    plt.xticks(rotation=45)
    aa_barplot.get_figure().savefig('plots/'+str(prefix)+'/egg_mutation_aa_prevalence_'+str(prefix)+'.pdf', bbox_inches='tight')

def plot_overall_aas(prefix):
    """
    Plot overall amino diversity at each site before egg-passaging. Indicate the proportion of each genotype that mutated during egg-passaging
    """

    df = pd.read_csv('dataframes/'+prefix+'_egg.csv')

    mut_sites = [col for col in df.columns if col[0:3]=='mut']

    before_aas = []

    for mut_site in mut_sites:
        bottom_proportion = 0.0
        stack = 0

        unpassaged_group = df.groupby(str(mut_site[3:])+'_lastnode').size().reset_index(name='count')

        for k, v in unpassaged_group.iterrows():

            site_group = df[df[str(mut_site[3:])+'_lastnode'] == v[str(mut_site[3:])+'_lastnode']].groupby(mut_site).size().reset_index(name='count')

            if len(site_group)!=2:
                site_group.loc[1] = [True, 0.0]

            for i, r in site_group.iterrows():
                proportion = float(r['count'])/float(len(df))
                before_aas.append({'site': mut_site, 'aa': v[str(mut_site[3:])+'_lastnode'], 'mutated': r[mut_site], 'proportion': proportion, 'bottom_proportion': bottom_proportion, 'stack': stack})
                bottom_proportion += proportion
                stack +=1

    before_aas_df = pd.DataFrame(before_aas)

    #Add x position for grouped bar chart
    pmap = {}
    x_pos = 2
    for pos in sorted(before_aas_df['site'].unique()):
        pmap[pos] = x_pos
        x_pos += 2

    cmap = {0:'#BA7C8E', 1:'#5c3d46', 2:'#C8EAD7', 3:'#99bfaa', 4:'#fbe3b1', 5:'#f8c968', 6: '#cdddf9', 7:'#75a1ef', 8: 'red', 9:'pink'}

    width = 1.2

    fig, ax = plt.subplots()

    plt.bar(before_aas_df['site'].map(pmap), before_aas_df['proportion'], width, bottom = before_aas_df['bottom_proportion'], color = before_aas_df['stack'].map(cmap))
    for i, r in before_aas_df[before_aas_df['mutated']==False].iterrows():
        if r['proportion'] >= 0.02:
            plt.annotate(r['aa'], xy=((pmap[r['site']]), (r['bottom_proportion']+0.5*r['proportion'])), color="black", ha='center', va='center', size=10)
        elif r['proportion'] < 0.02:
            plt.annotate(r['aa'], xy=((pmap[r['site']]+width/1.5), (r['bottom_proportion']+r['proportion'])), color="black", ha='center', va='center', size=6)


    ax.set(xlabel = 'HA site', ylabel='Genotype')
    ax.set_xticks([p for p in pmap.values()])
    ax.set_xticklabels([s[3:] for s in pmap.keys()])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.text(18, 0.97, 'Mutated during \negg-passaging', color='black', size=8, fontweight='bold')
    plt.text(18, 0.87, 'False', color='#BA7C8E', size=10, fontweight='bold')
    plt.text(19.2, 0.87, 'True', color='#5c3d46', size=10, fontweight='bold')
    plt.text(18, 0.79, 'False', color='#C8EAD7', size=10, fontweight='bold')
    plt.text(19.2, 0.79, 'True', color='#99bfaa', size=10, fontweight='bold')
    plt.text(18, 0.71, 'False', color='#fbe3b1', size=10, fontweight='bold')
    plt.text(19.2, 0.71, 'True', color='#f8c968', size=10, fontweight='bold')
    plt.text(18, 0.63, 'False', color='#cdddf9', size=10, fontweight='bold')
    plt.text(19.2, 0.63, 'True', color='#75a1ef', size=10, fontweight='bold')
    plt.text(6, 1.08, 'H3N2 strain genotypes', color='black', size=12)


    fig.savefig('plots/'+str(prefix)+'/genotypes_'+str(prefix)+'.pdf', bbox_inches='tight')

def find_mutation_aas(prefix):
    """
    Determine which mutations happen at each site. For each site, plot strains that mutated during egg passaging, showing the prevalence of each genotype before and after egg-passaging
    """

    df = pd.read_csv('dataframes/'+prefix+'_egg.csv')

    mut_sites = [col for col in df.columns if col[0:3]=='mut']

    mut_aas = []

    for mut_site in mut_sites:

        unpassaged_group = df[df[str(mut_site)] == True].groupby(str(mut_site[3:])+'_lastnode').size().sort_values(ascending=False).reset_index(name='count')
        mutation_group = df[df[str(mut_site)] == True].groupby(str(mut_site[3:])).size().sort_values(ascending=False).reset_index(name='count')

        for aa_count in range(len(unpassaged_group)):

            if aa_count == 0:
                mut_aas.append({'site': mut_site, 'aa': unpassaged_group[str(mut_site[3:])+'_lastnode'][aa_count], 'proportion': float(unpassaged_group['count'][aa_count])/float(len(df[df[str(mut_site)] == True])), 'stack': aa_count, 'strain': 'unpassaged', 'bottom_proportion': 0.0})
            #Find height where stacked bar should start
            else:
                prev_bars = unpassaged_group['count'][0]
                for n in range(1, aa_count):
                    prev_bars += unpassaged_group['count'][aa_count-n]
                mut_aas.append({'site': mut_site, 'aa': unpassaged_group[str(mut_site[3:])+'_lastnode'][aa_count], 'proportion': float(unpassaged_group['count'][aa_count])/float(len(df[df[str(mut_site)] == True])), 'stack': aa_count, 'strain': 'unpassaged', 'bottom_proportion': float(prev_bars)/float(len(df[df[str(mut_site)] == True]))})

        for aa_count in range(len(mutation_group)):

            if aa_count == 0:
                mut_aas.append({'site': mut_site, 'aa': mutation_group[str(mut_site[3:])][aa_count], 'proportion': float(mutation_group['count'][aa_count])/float(len(df[df[str(mut_site)] == True])), 'stack': aa_count, 'strain': 'egg', 'bottom_proportion': 0.0})
            #Find height where stacked bar should start
            else:
                prev_bars = mutation_group['count'][0]
                for n in range(1, aa_count):
                    prev_bars += mutation_group['count'][aa_count-n]
                mut_aas.append({'site': mut_site, 'aa': mutation_group[str(mut_site[3:])][aa_count], 'proportion': float(mutation_group['count'][aa_count])/float(len(df[df[str(mut_site)] == True])), 'stack': aa_count, 'strain': 'egg', 'bottom_proportion': float(prev_bars)/float(len(df[df[str(mut_site)] == True]))})

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
    for pos in sorted(mut_aas_df['site'].unique()):
        pmap[pos] = x_pos
        x_pos += 2

    stacks_u = mut_aas_df[mut_aas_df['strain']=='unpassaged']['stack'].unique()
    palette_u = dict(zip(stacks_u, sns.blend_palette(['#009a9a', '#00daa4'], n_colors = len(stacks_u))))

    stacks_e = mut_aas_df[mut_aas_df['strain']=='egg']['stack'].unique()
    palette_e = dict(zip(stacks_e, sns.blend_palette(['#f86968', '#f8c868', '#fada99'], n_colors= len(stacks_e))))
    # palette_e = dict(zip(stacks_e, sns.color_palette("OrRd", len(stacks_e))))

    width = 0.6

    fig, ax = plt.subplots()

    for strain in mut_aas_df['strain'].unique():
        for stack in range(len(mut_aas_df['stack'].unique())):
            stack_df = mut_aas_df[(mut_aas_df['stack']==int(stack)) & (mut_aas_df['strain']==str(strain))]
            if strain=='unpassaged':
                plt.bar(stack_df['site'].map(pmap), stack_df['proportion'], width, bottom = stack_df['bottom_proportion'], color = stack_df['stack'].map(palette_u))
                for i, r in stack_df.iterrows():
                    plt.annotate(r['aa'], xy=((pmap[r['site']]), (r['bottom_proportion']+0.3*r['proportion'])), color="black", ha='center')
            else:
                plt.bar((stack_df['site'].map(pmap)+width+0.05), stack_df['proportion'], width, bottom = stack_df['bottom_proportion'], color = stack_df['stack'].map(palette_e))
                for i, r in stack_df.iterrows():
                    plt.annotate(r['aa'], xy=((pmap[r['site']]+width+0.5), (r['bottom_proportion']+0.3*r['proportion'])), color="black", size=6, ha='center')


    ax.set(xlabel = 'HA site', ylabel='Genotype')
    ax.set_xticks([p for p in pmap.values()])
    ax.set_xticklabels([s[3:] for s in pmap.keys()])
    plt.text(20.5, 0.97, 'BEFORE \n egg-passaging', color='#009a9a', size=8, fontweight='bold')
    plt.text(20.5, 0.85, 'AFTER \n egg-passaging', color='#f86968', size=8, fontweight='bold')
    plt.text(3, 1.08, 'H3N2 strains that mutated during egg-passaging', color='black', size=12)


    fig.savefig('plots/'+str(prefix)+'/before_after_eggpassaging_'+str(prefix)+'.pdf', bbox_inches='tight')

def main(input_df):
    df_name = str.split(input_df, 'dataframes/')[1]
    prefix = str.split(df_name, '.csv')[0]
    plot_mutation_site(prefix)
    plot_mutation_aa(prefix)
    find_mutation_aas(prefix)
    plot_overall_aas(prefix)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Determines egg-specific mutations")
    parser.add_argument('--in_file', help= "input dataframe file")
    args = parser.parse_args()

    main(input_df = args.in_file)
