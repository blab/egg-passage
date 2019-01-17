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
        egg_v_unpassaged_obs = np.array([chi_obs['egg'], chi_obs['0']])
        eu_chi2, eu_p, eu_dof, eu_expected = chi2_contingency(egg_v_unpassaged_obs)
        chi_pvalues['egg_v_cell'] = ec_p
        chi_pvalues['egg_v_unpassaged'] = eu_p

        mut_prev[mut_site] = mut_prev_site
        chi_square[mut_site] = chi_pvalues


    mut_prev_df = pd.DataFrame(mut_prev)
    chi_square_df = pd.DataFrame(chi_square)

    return mut_prev_df, chi_square_df


def plot_mutation_prev(prefix):

    mut_prev_df, chi_square_df = find_mutation_prevalence(prefix)

    #Make dictionary to store significance of chi-squared tests at each site
    site_sig = {}
    for site, stats in chi_square_df.items():
        stat_sig = 0 #statistical significance (0 if neither egg_v_cell or egg_v_unpassaged is significant, 2 if both are, 1 if only 1 is)
        for stat in stats:
            if stat < 0.0001:
                stat_sig+=1
        site_sig[site[3:]] = stat_sig

    plot_df = mut_prev_df.unstack().reset_index().rename(columns={'level_0':'site', 'level_1':'virus_passage', 0:'prevalence'})
    plot_df['site'] = plot_df['site'].str[3:]
    plot_df['virus_passage'] = np.where(plot_df['virus_passage']=='0', 'unpassaged', plot_df['virus_passage'])

    sns.set(style="whitegrid")
    barplot = sns.barplot(x='site', y='prevalence', hue='virus_passage', data=plot_df)
    barplot.set(xlabel='HA position', ylabel='prevalence of mutation')
    barplot.get_figure().savefig('plots/egg_mutation_prevalence_'+str(prefix)+'.png')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Determines egg-specific mutations")
    parser.add_argument('--prefix', default= 'h3n2_6y_notiter', help= "specify prefix for naming data files")
    args = parser.parse_args()

    plot_mutation_prev(prefix = args.prefix)
