"""
Create .csv file containing relevant information from augur output. Find most prevalent tip
mutations in egg-passaged strains and extract the genotype of each sequence at these positions
"""

import argparse
import re
import json
import ast
import pandas as pd
import numpy as np
from Bio import SeqIO


def organize_egg_data_csv(prefix, positions, tree_path, sequences_path):

    #Load tree file
    with open(tree_path, 'r') as jsonfile:
        tree = json.load(jsonfile)

    positions = ast.literal_eval(positions)

    tip_muts = {}

    #Get root sequence
    root_seq = tree['translations']['HA1']

    def traverse(branch, seq, pos_list):

        #keep track of mutations at internal nodes
        if 'children' in branch.keys():
            for child in branch['children']:
                if 'aa_muts' in child.keys():
                    traverse_level.append(child['aa_muts']['HA1'])
                    traverse(child, seq, pos_list)
                    traverse_level.remove(child['aa_muts']['HA1'])
                else:
                    traverse(child, seq, pos_list)

        elif 'children' not in branch.keys():
            if 'aa_muts' in branch.keys():
                traverse_level.append(branch['aa_muts']['HA1'])

            muts_list = [str(mut) for sublist in traverse_level for mut in sublist]

            tip_sequence = seq
            for mut in muts_list:
                internal_mut_pos = int(re.findall('\d+', mut)[0])
                internal_mut_aa = mut[-1:]
                tip_sequence = tip_sequence[:internal_mut_pos-1] + internal_mut_aa + tip_sequence[internal_mut_pos:]

            tip_muts[branch['strain']]=[branch['aa_muts']['HA1'], branch['aa_muts']['HA2'],
                                        branch['aa_muts']['SigPep'],branch['attr']['num_date'],
                                        (branch['attr']['dTiterSub'] if 'dTiterSub' in branch['attr'] else None),
                                        (branch['attr']['cTiterSub'] if 'cTiterSub' in branch['attr'] else None),
                                        branch['attr']['clade_membership']] + [tip_sequence[pos-1] for pos in pos_list]

            if 'aa_muts' in branch.keys():
                traverse_level.remove(branch['aa_muts']['HA1'])

    traverse_level = []
    traverse(tree, root_seq, positions)

    df = pd.DataFrame(tip_muts).T
    df.reset_index(inplace=True)
    df.columns = ['strain', 'tip_HA1_muts', 'tip_HA2_muts', 'tip_SigPep_muts', 'date','dTiterSub','cTiterSub', 'clade'] + positions
    df['dTiterSub'], df['cTiterSub']= df['dTiterSub'].astype(float, inplace=True), df['cTiterSub'].astype(float, inplace=True)
    df['passage'] = np.select((df.strain.str.contains('egg'), df.strain.str.contains('cell')), ('egg', 'cell'))
    #Identify pairs where strain sequence exists for multiple passage types
    df['source'] = np.select((df.passage=='egg', df.passage=='cell', df.passage=='0'),
                             (df.strain.str.replace('-egg',''), df.strain.str.replace('-cell',''), df.strain))
    e_u_df = df[(df['passage']=='egg') | (df['passage']=='0')]
    pairs_u = e_u_df[e_u_df.duplicated(subset='source', keep=False)]
    e_c_df = df[(df['passage']=='egg') | (df['passage']=='cell')]
    pairs_c = e_c_df[e_c_df.duplicated(subset='source', keep=False)]
    pairs = pd.concat([pairs_u, pairs_c])
    pair_ids = dict(zip(list(pairs['source'].unique()),[[n+1] for n in range(len(pairs['source'].unique()))]))
    pair_ids = pd.DataFrame(pair_ids).T.reset_index().rename(columns={'index':'source', 0:'pair_id'})
    df = df.merge(pair_ids, on='source', how='left')

    #!!!! Must first check all sites of insterest to make sure clades have predeominantly one genotype!!!!
    #At each position, find predominant genotype of circulating virus by clade
    #Use this to determine whether egg-passaged strains have mutated
    #!!!! Must first check all sites of insterest to make sure clades have predeominantly one genotype!!!!

    clade_gtype = {}
    for p in positions:
        clade_gtype_pos = {}
        for c_name, clade in df[df['passage']=='0'].groupby('clade'):
            clade_gtype_pos[c_name] = str(clade[p].value_counts().nlargest(1))[0]
        clade_gtype[p] = clade_gtype_pos

    for p in positions:
        df['circulating'+str(p)] = df['clade'].map(clade_gtype[p])
    #Determine whether there sequence has mutated relative to clade background, at each position
    for p in positions:

        df['mut'+str(p)] = np.select(
        (df[p]==df['circulating'+str(p)], df[p]!=df['circulating'+str(p)]),
        (False, True))

    #Save organized data to a .csv
    df.to_csv('data/'+prefix+'_df.csv')

    #Make egg-seq only DF and save to .csv
    egg_df = df[df['passage']=='egg']
    egg_df.to_csv('data/'+prefix+'_egg_df.csv')



    #make tidy version of df where each mutation gets a row
    mut_df = pd.DataFrame(columns=['mutation']+ list(df.columns))

    count=0
    for i, r in df.iterrows():

        for ha1_mut in r['tip_HA1_muts']:
            mut_df.loc[count]= ['HA1'+str(ha1_mut)] + list(df.loc[i])
            count+=1

        for ha2_mut in r['tip_HA2_muts']:
            mut_df.loc[count]= ['HA2'+str(ha2_mut)] + list(df.loc[i])
            count+=1

        for sp_mut in r['tip_SigPep_muts']:
            mut_df.loc[count]= ['SP'+str(sp_mut)] + list(df.loc[i])
            count+=1

    #Save organized data to a .csv
    mut_df.to_csv('data/'+prefix+'_df_tidy.csv')

def find_tip_mutations(prefix):
    """
    Find egg-passaging-specific mutations and re-run organize_egg_data_csv with these positions
    """
    tip_df = pd.read_csv('data/'+prefix+'_df_tidy.csv')

    top_muts = {}
    for pas_type in tip_df['passage'].unique():
        top = (tip_df[tip_df.passage==pas_type].groupby('mutation')['mutation']
                ).count().sort_values(ascending=False)[:10]
        top_muts[pas_type] = list((g_name, g) for g_name, g in top.iteritems())

    #Extract positions of most prevalent egg-passaged tip mutations
    egg_tip_mutations = str(list(set([int(x[0].split('HA1')[1][1:-1]) for x in top_muts['egg']])))

    return egg_tip_mutations

def main(prefix, positions, tree_path, sequences_path, tip_mutations):
    organize_egg_data_csv(prefix, positions, tree_path, sequences_path)

    if tip_mutations == True:
        egg_tip_mutations = find_tip_mutations(prefix)
        if sorted(egg_tip_mutations) != sorted(positions):
            organize_egg_data_csv(prefix= 'tip_refine_'+str(prefix), positions= egg_tip_mutations, tree_path= tree_path, sequences_path= sequences_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description= "Organize Augur output into .csv to analyze egg-specific mutations")
    parser.add_argument('--prefix', default= 'h3n2_6y_hi', help= "specify prefix for naming data files")
    parser.add_argument('-pos', '--positions', default = '[160, 194, 186, 225, 219, 203, 156, 138]', help="specify a list of HA1 positions to analyze")
    parser.add_argument('-tree','--tree_path', default= 'augur/hi/flu_seasonal_h3n2_ha_6y_hi_tree.json', help= "specify the filepath to _tree.json file")
    parser.add_argument('-seqs','--sequences_path', default= 'augur/hi/flu_seasonal_h3n2_ha_6y_hi_sequences.json', help= "specify the filepath to _sequences.json file")
    parser.add_argument('--tip_mutations', default= True, help= "refine list of positions to consider based on the locations of the most prevalent tip mutations in egg-passaged sequences")
    args = parser.parse_args()

    main(prefix = args.prefix, positions = args.positions, tree_path = args.tree_path, sequences_path=args.sequences_path, tip_mutations = args.tip_mutations)
