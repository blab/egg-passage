import argparse
import json
import pandas as pd

def assign_clades(tree_path):

    with open(tree_path, 'r') as tree_json:
        tree = json.load(tree_json)

    leaf_paths = {}
    clock_lengths = {}

    def find_path(branch):

        node_path.append(branch['clade'])
        if branch['clade'] not in clock_lengths.keys():
            clock_lengths[str(branch['clade'])] = branch['attr']['clock_length']
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

    for internal_branch in reversed(range(0,max_path_length)):
        exclude_assigned = df[df['clade']=='unassigned']
        sub_df = exclude_assigned[exclude_assigned['path'].map(len) > internal_branch]

        path_group= sub_df.groupby((sub_df.path.apply(lambda col: col[0:(internal_branch+1)])).map(tuple))
        for k, v in path_group:
            if len(v)>=100:
                current_clade+=1
                df.at[v.index, 'clade']=current_clade

            elif len(v)>=50:
                if clock_lengths[str(k[-1])] >= 0.0008:
                    current_clade+=1
                    df.at[v.index, 'clade']=current_clade

            elif len(v)>=20:
                if clock_lengths[str(k[-1])] >= 0.003:
                    current_clade+=1
                    df.at[v.index, 'clade']=current_clade

    df = df.set_index('strain')

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

    args = parser.parse_args()

    assign_clades(tree_path = args.tree)
