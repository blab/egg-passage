import argparse
import json
import re
import ast
import pandas as pd
import numpy as np
from Bio import SeqIO

def organize_output(tree_path, seq_path, root_path, positions, prefix):

    #Load tree file
    with open(tree_path, 'r') as tree_json:
        tree = json.load(tree_json)

    #Load sequences file
    seqs = SeqIO.to_dict(SeqIO.parse(seq_path, "fasta"))

    #Load root seq
    with open(root_path, 'r') as root_json:
        root_seq = json.load(root_json)

    positions = ast.literal_eval(positions)

    tip_muts = {}
    def traverse(branch, seq, root, pos_list):

        #keep track of mutations at internal nodes
        if 'children' in branch.keys():
            for child in branch['children']:
                if 'HA1' in child['aa_muts']:
                    traverse_aa.append({str(child['clade']):child['aa_muts']['HA1']})
                    # traverse_aa.append(child['aa_muts']['HA1'])
                    # aa_mut_clade.append({str(child['clade']):child['aa_muts']['HA1']})
                    traverse(child, seq, root, pos_list)
                    traverse_aa.remove({str(child['clade']):child['aa_muts']['HA1']})
                    # traverse_aa.remove(child['aa_muts']['HA1'])
                    # aa_mut_clade.remove({str(child['clade']):child['aa_muts']['HA1']})

                else:
                    traverse(child, seq, root, pos_list)

        elif 'children' not in branch.keys():

            # muts_list = [str(mut) for sublist in traverse_aa for mut in sublist]
            # aa_mut_clade_list = [str(mut) for mut in aa_mut_clade]
            # last_node = [str(mut) for sublist in traverse_aa[:-1] for mut in sublist]

            aa_mut_clade_list = [str(mut) for mut in traverse_aa]
            muts_list = [str(mut) for sublist in [list(ast.literal_eval(x).values())[0]
                                                  for x in aa_mut_clade_list] for mut in sublist]
            last_node = muts_list

            if 'HA1' in branch['aa_muts']:
                branch_tip_muts = len(branch['aa_muts']['HA1'])
                last_node = last_node[:-(branch_tip_muts)]

            #Find sequence of tip and sequence one branch in
            tip_sequence = seq[branch['strain']].seq
            last_node_sequence = root_seq['HA1']

            for mut in last_node:
                internal_mut_pos = int(re.findall('\d+', mut)[0])
                internal_mut_aa = mut[-1:]
                last_node_sequence = last_node_sequence[:internal_mut_pos-1] + internal_mut_aa + last_node_sequence[internal_mut_pos:]

            tip_muts[branch['strain']]=([(branch['aa_muts']['HA1']
                                        if 'HA1' in branch['aa_muts'] else []),
                                        (branch['aa_muts']['HA2']
                                        if 'HA2' in branch['aa_muts'] else []),
                                        (branch['aa_muts']['SigPep']
                                        if 'SigPep' in branch['aa_muts'] else []),
                                        branch['attr']['num_date'],
                                        (branch['muts'] if 'muts' in branch else []),
                                        (branch['attr']['dTiterSub'] if 'dTiterSub' in branch['attr'] else None),
                                        (branch['attr']['cTiterSub'] if 'cTiterSub' in branch['attr'] else None),
                                        branch['attr']['clade_membership'],
                                        branch['attr']['kk_clade'],
                                        aa_mut_clade_list] +
                                        [tip_sequence[pos-1] for pos in pos_list] +
                                        [last_node_sequence[pos-1] for pos in pos_list])

    traverse_aa = []
    traverse(tree, seqs, root_seq, positions)

    df = pd.DataFrame(tip_muts).T
    df.reset_index(inplace=True)
    df.columns = ['strain', 'tip_HA1_muts', 'tip_HA2_muts', 'tip_SigPep_muts', 'date', 'tip_nt_muts', 'dTiterSub','cTiterSub', 'clade', 'kk_clade', 'aa_mut_list'] + positions + [str(x)+'_lastnode' for x in positions]
    df['dTiterSub'], df['cTiterSub']= df['dTiterSub'].astype(float, inplace=True), df['cTiterSub'].astype(float, inplace=True)
    df['passage'] = np.select((df.strain.str.contains('egg'), df.strain.str.contains('cell')), ('egg', 'cell'))
    df['passage'] = np.where(df['passage']=='0', 'unpassaged', df['passage'])
        #Identify pairs where strain sequence exists for multiple passage types
    df['source'] = np.select((df.passage=='egg', df.passage=='cell', df.passage=='unpassaged'),
                             (df.strain.str.replace('-egg',''), df.strain.str.replace('-cell',''), df.strain))
    e_u_df = df[(df['passage']=='egg') | (df['passage']=='unpassaged')]
    pairs_u = e_u_df[e_u_df.duplicated(subset='source', keep=False)]
    e_c_df = df[(df['passage']=='egg') | (df['passage']=='cell')]
    pairs_c = e_c_df[e_c_df.duplicated(subset='source', keep=False)]
    pairs = pd.concat([pairs_u, pairs_c])
    pair_ids = dict(zip(list(pairs['source'].unique()),[[n+1] for n in range(len(pairs['source'].unique()))]))
    pair_ids = pd.DataFrame(pair_ids).T.reset_index().rename(columns={'index':'source', 0:'pair_id'})
    df = df.merge(pair_ids, on='source', how='left')
    df['pair_id']= df['pair_id'].fillna(0)
    df['pair_id'] = df['pair_id'].astype(int, inplace=True)

    #Determine whether there sequence has mutated relative to ancestor 1 branch in, at each position
    for p in positions:
        df['mut'+str(p)] = np.select(
        (df[p]==df[str(p)+'_lastnode'], df[p]!=df[str(p)+'_lastnode']),
        (False, True))
    for p in positions:
        df['aa_mut'+str(p)] = np.where(df['mut'+str(p)]==1, df[str(p)+'_lastnode']+str(p)+df[p], None)

    #Find clusters of egg-passaged sequences,
    #allow mutations shared by these clusters to be called as mutations in egg strains
    #If multiple mutations occur at the same site within cluster, most recent mutation should be taken
    #Tips=cluster of size 1, so tip mutations override ancestral
    df['egg_muts'] = np.empty((len(df), 0)).tolist()
    max_internal_length=df['aa_mut_list'].map(len).max()

    for internal_branch in range(0,max_internal_length):
        sub_df = df[df['aa_mut_list'].map(len) > internal_branch]

        group= sub_df.groupby((sub_df.aa_mut_list.apply(lambda col: col[0:(internal_branch+1)])).map(tuple))
        for k, v in group:
            if len(v[v['passage']=='egg']) != 0:
                if len(v.groupby('passage')) == 1:
                    recent_muts = list(ast.literal_eval(k[-1]).values())[0]
                    for k_strain, v_strain in v.iterrows():
                        df.at[k_strain, 'egg_muts']+=recent_muts

                        #Find mutation(s) at specified positions
                        for recent_mut in recent_muts:
                            site = int(re.findall('\d+', recent_mut)[0])
                            if site in positions:
                                df.at[k_strain, 'mut' + str(site)] = 1
                                df.at[k_strain, 'aa_mut' + str(site)] = recent_mut
                                df.at[k_strain, str(site) + '_lastnode'] = recent_mut[0]

    #Save organized data to a .csv
    df.to_csv('dataframes/'+prefix+'.csv')

    #Make egg-seq only DF and save to .csv
    egg_df = df[df['passage']=='egg']
    egg_df.to_csv('dataframes/'+prefix+'_egg.csv')

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
    mut_df.to_csv('dataframes/'+prefix+'_tidy.csv')

def find_tip_mutations(prefix):
    """
    Find egg-passaging-specific mutations and re-run organize_output with these positions
    """
    tip_df = pd.read_csv('dataframes/'+prefix+'_tidy.csv')

    top_muts = {}
    for pas_type in tip_df['passage'].unique():
        top = (tip_df[tip_df.passage==pas_type].groupby('mutation')['mutation']
                ).count().sort_values(ascending=False)[:10]
        top_muts[pas_type] = list((g_name, g) for g_name, g in top.iteritems())

    #Extract positions of most prevalent egg-passaged tip mutations
    egg_tip_mutations = str(list(set([int(x[0].split('HA1')[1][1:-1]) for x in top_muts['egg']])))

    return egg_tip_mutations

def main(tree_path, seq_path, root_path, positions, tip_mutations):
    lineage = str.split(tree_path, '_')[2]
    segment = str.split(tree_path, '_')[3]
    resolution = str.split(tree_path, '_')[4]
    assay = str.split(tree_path, '_')[5]
    prefix = lineage + '_' + segment + '_' + resolution + '_' + assay
    organize_output(tree_path, seq_path, root_path, positions, prefix)

    if tip_mutations == True:
        egg_tip_mutations = find_tip_mutations(prefix)
        if sorted(ast.literal_eval(egg_tip_mutations)) != sorted(sorted(ast.literal_eval(positions))):
            organize_output(tree_path= tree_path, seq_path= seq_path, root_path= root_path, positions= egg_tip_mutations, prefix= prefix)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Extract and organize relevant data and results from augur output to analyze egg-passaging mutations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', help= "path to _tree.json file")
    parser.add_argument('--seqs', help= "path to aa-seq_ fasta file")
    parser.add_argument('--root_seq', help= "path to sequence of tree root")
    parser.add_argument('--positions', default = '[160, 194, 186, 225, 219, 203, 156, 138]', help="list of HA1 positions to analyze")
    parser.add_argument('--tip_mutations', default= True, help= "refine list of positions to consider based on the locations of the most prevalent tip mutations in egg-passaged sequences")

    args = parser.parse_args()

    main(tree_path = args.tree, seq_path = args.seqs, root_path = args.root_seq, positions = args.positions, tip_mutations= args.tip_mutations)
