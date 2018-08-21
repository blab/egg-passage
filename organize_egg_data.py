import pandas as pd
import numpy as np
import json
import ast
from Bio import SeqIO
import argparse

def organize_egg_data_csv(positions= '[160,194,186,225,219,203,156,138]', tree_path='auspice/flu_seasonal_h3n2_ha_6y_tree.json', sequences_path= 'auspice/flu_seasonal_h3n2_ha_6y_sequences.json', flu_ref_seq='h3n2_outgroup.gb'):

    #Load tree and sequence files
    with open(tree_path, 'r') as jsonfile:
        tree = json.load(jsonfile)

    with open(sequences_path, 'r') as jsonfile:
        seqs = json.load(jsonfile)

    #Load genbank file with reference HA sequence
    ref_names = {'hemagglutinin':'hemagglutinin','Signal peptide':'SP_ref','HA1 protein':'HA1_ref', 'HA2 protein':'HA2_ref'}
    ref_seqs = {}

    for seq_record in SeqIO.parse(flu_ref_seq, 'genbank'):
        for feature in seq_record.features:
            if feature.type=='CDS':
                ref_seqs[str(ref_names[feature.qualifiers['product'][0]])]= (feature.location.extract(seq_record).seq.translate())


    #Input desired 0-based HA1 amino acid residue positions, find reference amino acid
    positions = ast.literal_eval(positions)
    position_refaa = []

    for pos in positions:
        position_refaa.append((str(pos-1), ref_seqs['HA1_ref'][pos-1]))

    tip_muts = {}

    def traverse(branch, seq, pos_list):

        if 'children' not in branch.keys():

            tip_muts[branch['strain']]=[branch['aa_muts']['HA1'], branch['aa_muts']['HA2'],
                                        branch['aa_muts']['SigPep'],branch['attr']['num_date'],
                                        branch['attr']['dTiterSub'],
                                        branch['attr']['cTiterSub'],
                                        branch['attr']['clade_membership']] + [str(seq[str(branch['clade'])]['HA1'][pos]) if pos in seq[str(branch['clade'])]['HA1'] else str(ref_aa)
                                        for pos, ref_aa in pos_list]

        else:
            for child in branch['children']:
                traverse(child, seq, pos_list)


    traverse(tree, seqs, position_refaa)

    #Organize data in a DF
    df = pd.DataFrame(tip_muts).T
    df.reset_index(inplace=True)
    df.columns = ['strain', 'tip_HA1_muts', 'tip_HA2_muts', 'tip_SigPep_muts', 'date','dTiterSub','cTiterSub', 'clade']+positions
    df['passage'] = np.select((df.strain.str.contains('egg'), df.strain.str.contains('cell')), ('egg', 'cell'))
    df['source'] = np.select((df.passage=='egg', df.passage=='cell', df.passage=='0'),
                             (df.strain.str.replace('-egg',''), df.strain.str.replace('-cell',''), df.strain))
    df['dTiterSub'], df['cTiterSub']= df['dTiterSub'].astype(float, inplace=True), df['cTiterSub'].astype(float, inplace=True)

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
    df.to_csv('dataframes/auspice_df.csv')

    #Make egg-seq only DF and save to .csv
    egg_df = df[df['passage']=='egg']
    egg_df.to_csv('dataframes/egg_df.csv')



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
    mut_df.to_csv('dataframes/auspice_by_tipmut.csv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Organize Augur output into .csv to analyze egg-specific mutations")
    parser.add_argument('-p', '--positions', default = '[160,194,186,225,219,203,156,138]', help="specify a list of HA1 positions to analyze")
    parser.add_argument('-tree','--tree_path', default= 'auspice/flu_seasonal_h3n2_ha_6y_tree.json', help="specify the filepath to _tree.json file")
    parser.add_argument('-seqs','--sequences_path', default= 'auspice/flu_seasonal_h3n2_ha_6y_sequences.json', help="specify the filepath to _sequences.json file")
    parser.add_argument('--flu_ref_seq', default='h3n2_outgroup.gb', help="specify the filepath to genbank reference sequence")
    args = parser.parse_args()
    organize_egg_data_csv(positions = args.positions, tree_path = args.tree_path, sequences_path=args.sequences_path, flu_ref_seq=args.flu_ref_seq)
