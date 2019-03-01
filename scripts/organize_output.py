import argparse
import json
import re
import ast
import pandas as pd
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
                # if 'aa_muts' in child.keys():
                if child['mut']:
                    traverse_aa.append(child['aa_muts']['HA1'])
                    aa_mut_clade.append({str(child['clade']):child['aa_muts']['HA1']})
                    traverse(child, seq, root, pos_list)
                    traverse_aa.remove(child['aa_muts']['HA1'])
                    aa_mut_clade.remove({str(child['clade']):child['aa_muts']['HA1']})

                else:
                    #Append place holder for branches with no mutations
                    traverse_aa.append([])
                    traverse(child, seq, root, pos_list)
                    traverse_aa.remove([])

        elif 'children' not in branch.keys():

            muts_list = [str(mut) for sublist in traverse_aa for mut in sublist]
            aa_mut_clade_list = [str(mut) for mut in aa_mut_clade]
            last_node = [str(mut) for sublist in traverse_aa[:-1] for mut in sublist]

            #Find sequence of tip and sequence one branch in
            tip_sequence = seq[branch['strain']]
            last_node_sequence = root_seq['HA1']

            for mut in last_node:
                internal_mut_pos = int(re.findall('\d+', mut)[0])
                internal_mut_aa = mut[-1:]
                last_node_sequence = last_node_sequence[:internal_mut_pos-1] + internal_mut_aa + last_node_sequence[internal_mut_pos:]

            tip_muts[branch['strain']]=[branch['aa_muts']['HA1'],
                                        branch['aa_muts']['HA2'],
                                        branch['aa_muts']['SigPep'],
                                        branch['attr']['num_date'],
                                        (branch['muts'] if 'muts' in branch else None),
                                        (branch['attr']['dTiterSub'] if 'dTiterSub' in branch['attr'] else None),
                                        (branch['attr']['cTiterSub'] if 'cTiterSub' in branch['attr'] else None),
                                        branch['attr']['clade_membership'],
                                        aa_mut_clade_list] +
                                        [tip_sequence[pos-1] for pos in pos_list] +
                                        [last_node_sequence[pos-1] for pos in pos_list]

    traverse_aa = []
    aa_mut_clade = []
    traverse(tree, seqs, root_seq, positions)

    df = pd.DataFrame(tip_muts).T
    df.reset_index(inplace=True)
    df.columns = ['strain', 'tip_HA1_muts', 'tip_HA2_muts', 'tip_SigPep_muts', 'date', 'tip_nt_muts', 'dTiterSub','cTiterSub', 'clade', 'aa_mut_clade'] + positions + [str(x)+'_lastnode' for x in positions]
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

    #Determine whether there sequence has mutated relative to ancestor 1 branch in, at each position
    for p in positions:
        df['mut'+str(p)] = np.select(
        (df[p]==df[str(p)+'_lastnode'], df[p]!=df[str(p)+'_lastnode']),
        (False, True))
    for p in positions:
        df['aa_mut'+str(p)] = np.where(df['mut'+str(p)]==1, df[str(p)+'_lastnode']+str(p)+df[p], None)


def find_tip_mutations(prefix):
    """
    Find egg-passaging-specific mutations and re-run organize_output with these positions
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

def main(tree_path, seq_path, root_path, prefix, positions, tip_mutations):
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
    parser.add_argument('--prefix', help= "prefix for the csv output files")
    parser.add_argument('--positions', default = '[160, 194, 186, 225, 219, 203, 156, 138]', help="list of HA1 positions to analyze")
    parser.add_argument('--tip_mutations', default= True, help= "refine list of positions to consider based on the locations of the most prevalent tip mutations in egg-passaged sequences")

    args = parser.parse_args()

    main(tree_path = args.tree, seq_path = args.seqs, root_path = args.root_seq, prefix = args.prefix, positions = args.positions, tip_mutations= args.tip_mutations)
