import argparse
import json
import pandas as pd

def assign_clades(tree_path, output_csv, method):

    with open(tree_path, 'r') as tree_json:
        tree = json.load(tree_json)

    leaf_paths = {}
    branch_traits = {}

    def find_path(branch):

        node_path.append(branch['clade'])
        if branch['clade'] not in branch_traits.keys():
            if method == 'clock_length':
                branch_traits[str(branch['clade'])] = branch['attr']['clock_length']
            elif method == 'mutations':
                branch_traits[str(branch['clade'])] = branch['muts']
        if 'children' in branch.keys():
            for child in branch['children']:
                find_path(child)
                node_path.remove(child['clade'])


        elif 'children' not in branch.keys():

            node_path_list = [str(node) for node in node_path]
            leaf_paths[branch['strain']] = [branch['attr']['clade_membership'], node_path_list]

    node_path = []
    find_path(tree)


    df = pd.DataFrame(leaf_paths).T
    df.reset_index(inplace=True)
    df.columns = ['strain', 'clade_membership', 'path']
    df['clade'] = 'unassigned'

    #Find clades
    max_path_length=df['path'].map(len).max()

    current_clade = 0
    assigned_clades = {}

    for internal_branch in reversed(range(0,max_path_length)):
        exclude_assigned = df[df['clade']=='unassigned']
        sub_df = exclude_assigned[exclude_assigned['path'].map(len) > internal_branch]

        path_group= sub_df.groupby((sub_df.path.apply(lambda col: col[0:(internal_branch+1)])).map(tuple))
        for k, v in path_group:
            if len(v)>=100:
                current_clade+=1
                df.at[v.index, 'clade']=current_clade
                assigned_clades[current_clade] = {'clade_mrca':k[-1]}

            elif len(v)>=50:
                if method == 'clock_length':
                    if branch_traits[str(k[-1])] >= 0.0008:
                        current_clade+=1
                        df.at[v.index, 'clade']=current_clade
                        assigned_clades[current_clade] = {'clade_mrca':k[-1]}
                elif method == 'mutations':
                    if len(branch_traits[str(k[-1])]) >= 1:
                        current_clade+=1
                        df.at[v.index, 'clade']=current_clade
                        assigned_clades[current_clade] = {'clade_mrca':k[-1]}

            elif len(v)>=20:
                if method == 'clock_length':
                    if branch_traits[str(k[-1])] >= 0.003:
                        current_clade+=1
                        df.at[v.index, 'clade']=current_clade
                        assigned_clades[current_clade] = {'clade_mrca':k[-1]}
                elif method == 'mutations':
                    if len(branch_traits[str(k[-1])]) >= 4:
                        current_clade+=1
                        df.at[v.index, 'clade']=current_clade
                        assigned_clades[current_clade] = {'clade_mrca':k[-1]}

    df = df.set_index('strain')

    #Find mutations that occur on most recent common ascestor of kk_clades
    def find_defining_genotypes(branch):

        if 'children' in branch.keys():
            for child in branch['children']:
                for k,v in assigned_clades.items():
                    if str(branch['clade']) == v['clade_mrca']:
                        assigned_clades[k]['aa_muts'] = branch['aa_muts']
                        assigned_clades[k]['nt_muts'] = branch['muts']
                find_defining_genotypes(child)

    find_defining_genotypes(tree)
    assigned_clades_df = pd.DataFrame(assigned_clades).T.reset_index().rename(columns={'index':'kk_clade'})
    assigned_clades_df = assigned_clades_df[['kk_clade', 'clade_mrca', 'aa_muts', 'nt_muts']]
    assigned_clades_df.to_csv(output_csv, index=False)

    def add_clades(branch, clade_assignments):
        if 'children' in branch.keys():
            for child in branch['children']:
                add_clades(child, df)

        elif 'children' not in branch.keys():
            assigned_clade = df.loc[branch['strain']]['clade']
            branch['attr']['kk_clade'] = assigned_clade

    add_clades(tree, df)

    with open(tree_path, 'w') as outfile:
        json.dump(tree, outfile, indent=2, sort_keys = True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assign clades by grouping viruses into groups of at least 100 strains")

    parser.add_argument('--tree', help= "path to _tree.json file")
    parser.add_argument('--kk-clades-file', help= "path to output file listing clades and clade-defining mutations")
    parser.add_argument('--method', choices=["clock_length", "mutations"], default= "clock_length", help= "determine clades based on clade size and clock length of common ancestor OR number of mutations in common ancestor")

    args = parser.parse_args()

    assign_clades(tree_path = args.tree, output_csv = args.kk_clades_file, method = args.method)
