import argparse
import json
import pandas as pd
from ete3 import Tree

def assign_clades(tree_path, newick_path, output_csv, method):

    with open(tree_path, 'r') as tree_json:
        tree = json.load(tree_json)

    t = Tree(newick_path, format=1)

    leaf_paths = {}
    branch_traits = {}

    def find_path(branch):

        node_path.append(branch['strain'])
        if branch['strain'] not in branch_traits.keys():
            if method == 'clock_length':
                branch_traits[str(branch['strain'])] = branch['attr']['clock_length']
            elif method == 'mutations':
                branch_traits[str(branch['strain'])] = branch['aa_muts']
        if 'children' in branch.keys():
            for child in branch['children']:
                #Find correlation between num aa mutations and clock length
                num_aa_muts = sum([len(child['aa_muts'][x]) for x in child['aa_muts']])
                if num_aa_muts in clocklength_per_aamut.keys():
                    clocklength_per_aamut[num_aa_muts].append(child['attr']['clock_length'])
                else:
                    clocklength_per_aamut[num_aa_muts] = [child['attr']['clock_length']]

                find_path(child)
                node_path.remove(child['strain'])


        elif 'children' not in branch.keys():

            node_path_list = [str(node) for node in node_path]
            leaf_paths[branch['strain']] = [branch['attr']['clade_membership'], node_path_list]

    node_path = []
    clocklength_per_aamut = {}
    find_path(tree)


    df = pd.DataFrame(leaf_paths).T
    df.reset_index(inplace=True)
    df.columns = ['strain', 'clade_membership', 'path']
    df['kk_clade'] = 'unassigned'

    #Find correlation between num aa mutations and clock length
    avg_clocklength = {k: (float(sum(v))/len(v))
                       for k,v in clocklength_per_aamut.items() if len(v)>1}

    #Find clades
    max_path_length=df['path'].map(len).max()

    current_clade = 0
    assigned_clades = {}

    for internal_branch in reversed(range(0,max_path_length)):
        exclude_assigned = df[df['kk_clade']=='unassigned']
        sub_df = exclude_assigned[exclude_assigned['path'].map(len) > internal_branch]

        path_group= sub_df.groupby((sub_df.path.apply(lambda col: col[0:(internal_branch+1)])).map(tuple))
        for k, v in path_group:
            if method == 'clock_length':
                if len(v)>=100:
                    if branch_traits[str(k[-1])] >= 0.8*avg_clocklength[1]:
                        current_clade+=1
                        df.at[v.index, 'kk_clade'] = 'c'+str(current_clade)
                        assigned_clades['c'+str(current_clade)] = {'clade_mrca':k[-1]}
                elif len(v)>=50:
                    if branch_traits[str(k[-1])] >= 0.8*avg_clocklength[2]:
                        current_clade+=1
                        df.at[v.index, 'kk_clade'] = 'c'+str(current_clade)
                        assigned_clades['c'+str(current_clade)] = {'clade_mrca':k[-1]}
                elif len(v)>=20:
                    if branch_traits[str(k[-1])] >= 0.8*avg_clocklength[3]:
                        current_clade+=1
                        df.at[v.index, 'kk_clade'] = 'c'+str(current_clade)
                        assigned_clades['c'+str(current_clade)] = {'clade_mrca':k[-1]}

            elif method == 'mutations':
                if len(branch_traits[str(k[-1])]) >= 3:
                    if len(v)>=20:
                        current_clade+=1
                        df.at[v.index, 'kk_clade'] = 'c'+str(current_clade)
                        assigned_clades['c'+str(current_clade)] = {'clade_mrca':k[-1]}
                if len(branch_traits[str(k[-1])]) >= 2:
                    if len(v)>=50:
                        current_clade+=1
                        df.at[v.index, 'kk_clade'] = 'c'+str(current_clade)
                        assigned_clades['c'+str(current_clade)] = {'clade_mrca':k[-1]}
                elif len(branch_traits[str(k[-1])]) >= 1:
                    if len(v)>=100:
                        current_clade+=1
                        df.at[v.index, 'kk_clade'] = 'c'+str(current_clade)
                        assigned_clades['c'+str(current_clade)] = {'clade_mrca':k[-1]}

    df = df.set_index('strain')

    #Find mutations that occur on most recent common ascestor of kk_clades
    def find_defining_genotypes(branch):

        if 'children' in branch.keys():
            for child in branch['children']:
                for k,v in assigned_clades.items():
                    if str(branch['strain']) == v['clade_mrca']:
                        assigned_clades[k]['aa_muts'] = branch['aa_muts']
                        assigned_clades[k]['nt_muts'] = branch['muts']
                find_defining_genotypes(child)

    find_defining_genotypes(tree)
    assigned_clades_df = pd.DataFrame(assigned_clades).T.reset_index().rename(columns={'index':'kk_clade'})
    assigned_clades_df = assigned_clades_df[['kk_clade', 'clade_mrca', 'aa_muts', 'nt_muts']]
    assigned_clades_df.to_csv(output_csv, index=False)

    kkclade_nodes = {}
    for k,v in assigned_clades.items():
        for n in t.traverse(strategy='postorder'):
            if n.name == v['clade_mrca']:
                descendents = n.get_descendants()
                for d_node in descendents:
                    if d_node.name not in kkclade_nodes.keys():
                        kkclade_nodes[d_node.name] = k

    def add_clades(branch):
        if branch['strain'] not in kkclade_nodes.keys():
            assigned_clade = 'unassigned'
        elif branch['strain'] in kkclade_nodes.keys():
            assigned_clade = kkclade_nodes[branch['strain']]

        branch['attr']['kk_clade'] = assigned_clade

        if 'children' in branch.keys():
            for child in branch['children']:
                add_clades(child)

    add_clades(tree)

    with open(tree_path, 'w') as outfile:
        json.dump(tree, outfile, indent=2, sort_keys = True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assign clades by grouping viruses into groups of at least 100 strains")

    parser.add_argument('--tree', help= "path to _tree.json file")
    parser.add_argument('--newick-path', help= "path to .newick tree file")
    parser.add_argument('--kk-clades-file', help= "path to output file listing clades and clade-defining mutations")
    parser.add_argument('--method', choices=["clock_length", "mutations"], default= "clock_length", help= "determine clades based on clade size and clock length of common ancestor OR number of mutations in common ancestor")

    args = parser.parse_args()

    assign_clades(tree_path = args.tree, newick_path=args.newick_path, output_csv = args.kk_clades_file, method = args.method)
