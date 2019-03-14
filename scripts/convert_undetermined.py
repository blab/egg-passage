import argparse
import json

def convert_undetermined(tree_path):

    with open(tree_path, 'r') as tree_json:
        tree = json.load(tree_json)

    def traverse(branch):
        if 'children' in branch.keys():
            for child in branch['children']:
                traverse(child)

        elif 'children' not in branch.keys():
            if branch['attr']['passage']=='undetermined':
                branch['attr']['passage']='unpassaged'

    traverse(tree)

    with open(tree_path, 'w') as outfile:
        json.dump(tree, outfile, indent=2, sort_keys = True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert undetermined passage type to unpassaged")

    parser.add_argument('--tree', help= "path to _tree.json file")

    args = parser.parse_args()

    convert_undetermined(tree_path = args.tree)
