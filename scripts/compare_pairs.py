"""
Use pairs of strains that were sequenced before/after egg-passaging to assess error-rate of mutation inference method
"""

import argparse, json, ast, re
import pandas as pd
import numpy as np
from Bio import SeqIO

def find_paired_mutations(prefix, output, seq_path):
    """
    Analysis of strains where sequences are available for the virus before and after egg-passaging.
    """

    #initialize dict to write json
    pairs_json = {'paired_egg_viruses': {}, 'total_egg_viruses': {}}

    #Load dataframe, sequence file
    seqs = SeqIO.to_dict(SeqIO.parse(seq_path, "fasta"))
    all_df = pd.read_csv('dataframes/'+prefix+'.csv')
    egg_df = all_df[all_df['passage']=='egg']

    #filter data for only paired sequences
    df = all_df[all_df['pair_id']!=0]

    positions = [col[3:] for col in all_df.columns if col[0:3]=='mut']

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

    #Find false positives (mutation inferred, but strain is not mutated)
    num_muts_inferred = 0
    num_false_pos = 0
    false_pos_strain = []
    num_muts_inferred_limitsites = 0
    num_false_pos_limitsites = 0

    for k,v in pairs_df.iterrows():
        for egg_mut in ast.literal_eval(v['egg_muts']):

            egg_mut_pos = int(re.findall('\d+', egg_mut)[0])
            egg_aa = seqs[v['source']+'-egg'][(egg_mut_pos-1)]

            num_muts_inferred+=1
            if str(egg_mut_pos) in positions:
                num_muts_inferred_limitsites+=1

            if v['unpassaged_pair']==True:
                u_aa = seqs[v['source']][(egg_mut_pos-1)]
                if u_aa == egg_aa:
                    num_false_pos+=1
                    if str(egg_mut_pos) in positions:
                        num_false_pos_limitsites+=1
                        false_pos_strain.append({'strain':(v['source']+'-egg'), 'paired_strain': v['source'], 'mutation': egg_mut})

            if v['cell_pair']==True:
                cell_aa = seqs[v['source']+'-cell'][(egg_mut_pos-1)]
                if cell_aa == egg_aa:
                    num_false_pos+=1
                    if str(egg_mut_pos) in positions:
                        num_false_pos_limitsites+=1
                        false_pos_strain.append({'strain':(v['source']+'-egg'), 'paired_strain': (v['source']+'-cell'), 'mutation': egg_mut})

    pairs_json['paired_egg_viruses']['num_strains_with_pair'] = len(pairs_df)
    pairs_json['paired_egg_viruses']['all_ha1_mutations'] = []
    pairs_json['paired_egg_viruses']['top_sites_only'] = []
    pairs_json['paired_egg_viruses']['all_ha1_mutations'].append({'num_inferred_mutations':num_muts_inferred, 'num_false_positives':num_false_pos, 'false_positive_rate': num_false_pos/num_muts_inferred})
    pairs_json['paired_egg_viruses']['top_sites_only'].append({'num_inferred_mutations':num_muts_inferred_limitsites, 'num_false_positives':num_false_pos_limitsites, 'false_positive_rate': num_false_pos_limitsites/num_muts_inferred_limitsites, 'false_positive_strains': false_pos_strain})

    #Extrapolate false positive numbers to ALL egg strains
    pairs_json['total_egg_viruses']['total_num_egg_strains'] = len(egg_df)

    total_muts_inferred = 0
    total_muts_inferred_limitsites = 0

    for k,v in egg_df.iterrows():
        for egg_mut in ast.literal_eval(v['egg_muts']):
            egg_mut_pos = int(re.findall('\d+', egg_mut)[0])
            total_muts_inferred+=1
            if str(egg_mut_pos) in positions:
                total_muts_inferred_limitsites+=1

    est_false_pos = total_muts_inferred*(num_false_pos/num_muts_inferred)
    est_false_pos_limitsites = total_muts_inferred_limitsites*(num_false_pos_limitsites/num_muts_inferred_limitsites)

    pairs_json['total_egg_viruses']['all_ha1_mutations'] = []
    pairs_json['total_egg_viruses']['top_sites_only'] = []
    pairs_json['total_egg_viruses']['all_ha1_mutations'].append({'num_inferred_mutations':total_muts_inferred, 'estimated_num_false_positives':est_false_pos})
    pairs_json['total_egg_viruses']['top_sites_only'].append({'num_inferred_mutations':total_muts_inferred_limitsites, 'estimated_num_false_positives':est_false_pos_limitsites})

    #Find false negatives (strain is mutated, but mutation not inferred)
    num_muts_direct = 0
    num_false_neg = 0
    num_muts_direct_limitsites = 0
    num_false_neg_limitsites = 0

    for k,v in pairs_df.iterrows():

        egg_ha1 = seqs[v['source']+'-egg']

        for residue in range(len(egg_ha1)):
            egg_aa = egg_ha1[int(residue)-1]

            if v['unpassaged_pair']==True:
                u_ha1 = seqs[v['source']]
                u_aa = u_ha1[int(residue)-1]

                #Only want egg muts, not unpassaged muts
                u_tip_muts = df[df['strain']==(v['source'])]['tip_HA1_muts'].item()
                if str(residue) not in u_tip_muts:
                    if egg_aa != u_aa:
                        num_muts_direct+=1
                        u_mutation = u_aa + str(residue) + egg_aa
                        if u_mutation not in ast.literal_eval(v['egg_muts']):
                            num_false_neg+=1

                        if str(residue) in positions:
                            num_muts_direct_limitsites+=1
                            if u_mutation not in ast.literal_eval(v['egg_muts']):
                                num_false_neg_limitsites+=1


            if v['cell_pair']==True:
                cell_ha1 = seqs[v['source']+'-cell']
                cell_aa = cell_ha1[int(residue)-1]

                #Only want egg muts, not cell muts
                cell_tip_muts = df[df['strain']==(v['source']+'-cell')]['tip_HA1_muts'].item()
                if str(residue) not in cell_tip_muts:
                    if egg_aa != cell_aa:
                        num_muts_direct+=1
                        cell_mutation = cell_aa + str(residue) + egg_aa
                        if cell_mutation not in ast.literal_eval(v['egg_muts']):
                            num_false_neg+=1

                        if str(residue) in positions:
                            num_muts_direct_limitsites+=1
                            if cell_mutation not in ast.literal_eval(v['egg_muts']):
                                num_false_neg_limitsites+=1

    pairs_json['paired_egg_viruses']['all_ha1_mutations'].append({'real_num_mutations':num_muts_direct, 'num_false_negatives':num_false_neg, 'false_negative_rate': num_false_neg/num_muts_direct})
    pairs_json['paired_egg_viruses']['top_sites_only'].append({'real_num_mutations':num_muts_direct_limitsites, 'num_false_negatives':num_false_neg_limitsites, 'false_negative_rate': num_false_neg_limitsites/num_muts_direct_limitsites})

    #Estimate number false negatives in ALL egg strains
    num_paired_egg_seqs = len(pairs_df)
    num_total_egg_seqs = len(egg_df)

    est_total_mutations = (num_total_egg_seqs/num_paired_egg_seqs) * num_muts_direct
    est_total_mutations_limitsites = (num_total_egg_seqs/num_paired_egg_seqs) * num_muts_direct_limitsites
    est_false_neg = (num_total_egg_seqs/num_paired_egg_seqs)* num_false_neg
    est_false_neg_limitsites = (num_total_egg_seqs/num_paired_egg_seqs)* num_false_neg_limitsites

    pairs_json['total_egg_viruses']['all_ha1_mutations'].append({'estimated_num_mutations':est_total_mutations, 'estimated_num_false_negatives':est_false_neg})
    pairs_json['total_egg_viruses']['top_sites_only'].append({'estimated_num_mutations':est_total_mutations_limitsites, 'estimated_num_false_negatives':est_false_neg_limitsites})


    with open(output, 'w') as fh:
        json.dump(pairs_json, fh, indent=4)


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
