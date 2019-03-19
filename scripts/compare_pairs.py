"""
Compare methods for determining egg-passaging mutations
"""

import argparse, json
import pandas as pd
import numpy as np

def find_paired_mutations(prefix, output):
    """
    Analysis of strains where sequences are available for the virus before and after egg-passaging (or after cell-passaging and after egg-passaging)
    """

    #initialize dict to write json
    pairs_json = {'paired_egg_viruses': {}, 'total_egg_viruses': {}}
    all_df = pd.read_csv('dataframes/'+prefix+'.csv')
    #filter data for only paired sequences
    df = all_df[all_df['pair_id']!=0]

    positions = [col[3:] for col in df.columns if col[0:3]=='mut']
    mut_sites = [col for col in df.columns if col[0:3]=='mut']

    #Re-organize DF to one row per pair
    sub_egg = df[df['passage']=='egg'][['source'] + [str(pos) for pos in positions] + [('mut'+str(pos)) for pos in positions]].rename(columns = dict((str(pos), (str(pos)+'_egg')) for pos in positions))
    sub_u = df[df['passage']=='unpassaged'][['source'] + [str(pos) for pos in positions]].rename(columns = dict((str(pos), (str(pos)+'_u')) for pos in positions))
    sub_cell = df[df['passage']=='cell'][['source'] + [str(pos) for pos in positions]].rename(columns = dict((str(pos), (str(pos)+'_cell')) for pos in positions))

    pairs_u_df = sub_egg.merge(sub_u)
    pairs_cell_df = sub_egg.merge(sub_cell)
    pairs_cell_u_df = sub_u.merge(sub_cell)
    pairs_df = pairs_u_df.merge(pairs_cell_df, how='outer')

    #Find number of mutations mis-called by phylogenetic method
    total_muts_called_pairs = 0
    for mut_site in mut_sites:
        total_muts_called_pairs+=len(pairs_df[pairs_df[mut_site]==True])

    miscalled_muts = 0
    pairs_json['paired_egg_viruses']['miscalled_muts'] = []
    for k,v in pairs_df.iterrows():
        for pos in positions:
            if v['mut'+str(pos)]==True:
                if v[str(pos)+'_egg'] == v[str(pos)+'_u']:
                    miscalled_muts+=1
                    pairs_json['paired_egg_viruses']['miscalled_muts'].append({'strain': v['source'], 'HA1_position': pos, 'genotype': v[str(pos)+'_egg'], 'paired_virus_passage': 'unpassaged'})
                elif v[str(pos)+'_egg'] == v[str(pos)+'_cell']:
                    miscalled_muts+=1
                    pairs_json['paired_egg_viruses']['miscalled_muts'].append({'strain': v['source'], 'HA1_position': pos, 'genotype': v[str(pos)+'_egg'], 'paired_virus_passage': 'cell'})

    pairs_json['paired_egg_viruses']['num_inferred_muts'] = total_muts_called_pairs
    pairs_json['paired_egg_viruses']['num_miscalled_muts'] = miscalled_muts

    #Estimate number of mis-called mutations overall
    total_muts_called = 0
    for mut_site in mut_sites:
        total_muts_called+=len(all_df[all_df[mut_site]==True])

    miscalled_muts_est = total_muts_called*(miscalled_muts/total_muts_called_pairs)

    pairs_json['total_egg_viruses']['num_inferred_muts'] = total_muts_called
    pairs_json['total_egg_viruses']['estimated_num_miscalled_muts'] = miscalled_muts_est


    with open(output, 'w') as fh:
        json.dump(pairs_json, fh, indent=4)


def main(input_df):
    df_name = str.split(input_df, 'dataframes/')[1]
    prefix = str.split(df_name, '.csv')[0]
    output = "plots/egg_results/egg_mutation_accuracy_" + str(prefix) + ".json"
    find_paired_mutations(prefix, output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Compare methods for identifying egg-passaged\
                                                    mutations: 1) inferred by phylogeny vs. \
                                                    2) direct comparison of paired sequences")
    parser.add_argument('--in_file', help= "input dataframe file")
    args = parser.parse_args()

    main(input_df = args.in_file)
