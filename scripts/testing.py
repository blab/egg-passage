"""
Compare methods for determining egg-passaging mutations
"""

import argparse, json, ast, re
import pandas as pd
import numpy as np
from Bio import SeqIO

def find_paired_mutations(prefix, output, seq_path):
    """
    Analysis of strains where sequences are available for the virus before and after egg-passaging (or after cell-passaging and after egg-passaging)
    """

    #Load dataframe, sequence file
    seqs = SeqIO.to_dict(SeqIO.parse(seq_path, "fasta"))
    all_df = pd.read_csv('dataframes/'+prefix+'.csv')

    positions = [col[3:] for col in all_df.columns if col[0:3]=='mut']

    #filter data for only paired sequences
    df = all_df[all_df['pair_id']!=0]
    groups = all_df.groupby('pair_id')

    #Re-organize DF to one row per pair
    sub_egg = df[df['passage']=='egg'][['source', 'egg_muts']]
    sub_u = df[df['passage']=='unpassaged'][['source', 'strain']].rename(columns = {'strain':'unpassaged_pair'})
    sub_u['unpassaged_pair'] = True
    sub_cell = df[df['passage']=='cell'][['source', 'strain']].rename(columns = {'strain':'cell_pair'})
    sub_cell['cell_pair'] = True

    pairs_u_df = sub_egg.merge(sub_u)
    pairs_cell_df = sub_egg.merge(sub_cell)
    pairs_cell_u_df = sub_u.merge(sub_cell)
    pairs_df = pairs_u_df.merge(pairs_cell_df, how='outer')

    #Find number of false positive inferred muts (as determined by paired sequences)
    num_muts = 0
    num_miscall = 0
    num_muts_pos = 0
    num_miscall_pos = 0

    for k,v in pairs_df.iterrows():
        for egg_mut in ast.literal_eval(v['egg_muts']):

            egg_mut_pos = int(re.findall('\d+', egg_mut)[0])
            egg_aa = seqs[v['source']+'-egg'][(egg_mut_pos-1)]

            num_muts+=1
            if str(egg_mut_pos) in positions:
                num_muts_pos+=1

            if v['unpassaged_pair']==True:
                u_aa = seqs[v['source']][(egg_mut_pos-1)]
                if u_aa == egg_aa:
                    num_miscall+=1
                    if str(egg_mut_pos) in positions:
                        num_miscall_pos+=1

            if v['cell_pair']==True:
                cell_aa = seqs[v['source']+'-cell'][(egg_mut_pos-1)]
                if cell_aa == egg_aa:
                    num_miscall+=1
                    if str(egg_mut_pos) in positions:
                        num_miscall_pos+=1
    print(num_miscall)
    print(num_miscall_pos)


    #     print(v['source'])
    #     tip_sequence = seqs[(v['source']+'-egg')].seq
    #     print(tip_sequence[159])
    #     for pos in positions:
    #         if v['mut'+str(pos)]==True:
    #             if v[str(pos)+'_egg'] == v[str(pos)+'_u']:
    #                 miscalled_muts+=1
    #                 pairs_json['paired_egg_viruses']['miscalled_muts'].append({'strain': v['source'], 'HA1_position': pos, 'genotype': v[str(pos)+'_egg'], 'paired_virus_passage': 'unpassaged'})
    #             elif v[str(pos)+'_egg'] == v[str(pos)+'_cell']:
    #                 miscalled_muts+=1
    #                 pairs_json['paired_egg_viruses']['miscalled_muts'].append({'strain': v['source'], 'HA1_position': pos, 'genotype': v[str(pos)+'_egg'], 'paired_virus_passage': 'cell'})
    #
    # pairs_json['paired_egg_viruses']['num_inferred_muts'] = total_muts_called_pairs
    # pairs_json['paired_egg_viruses']['num_miscalled_muts'] = miscalled_muts
    #
    # #Estimate number of mis-called mutations overall
    # total_muts_called = 0
    # for mut_site in mut_sites:
    #     total_muts_called+=len(all_df[all_df[mut_site]==True])
    #
    # miscalled_muts_est = total_muts_called*(miscalled_muts/total_muts_called_pairs)
    #
    # pairs_json['total_egg_viruses']['num_inferred_muts'] = total_muts_called
    # pairs_json['total_egg_viruses']['estimated_num_miscalled_muts'] = miscalled_muts_est
    #
    #
    # with open(output, 'w') as fh:
    #     json.dump(pairs_json, fh, indent=4)


def main(input_df, seq_path):
    df_name = str.split(input_df, 'dataframes/')[1]
    prefix = str.split(df_name, '.csv')[0]
    output = "plots/egg_results/egg_mutation_accuracy_" + str(prefix) + ".json"
    find_paired_mutations(prefix, output, seq_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Compare methods for identifying egg-passaged\
                                                    mutations: 1) inferred by phylogeny vs. \
                                                    2) direct comparison of paired sequences")
    parser.add_argument('--in_file', help= "input dataframe file")
    parser.add_argument('--seqs', help= "path to aa-seq_ fasta file")
    args = parser.parse_args()

    main(input_df = args.in_file, seq_path = args.seqs)
