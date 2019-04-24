import argparse
import json, ast
import pandas as pd
from Bio import SeqIO
from ete3 import Tree

def main(tree_path, meta_path, seq_path, newick_tree, config_path, dataframe, sites, manual_colors=True):

    with open(tree_path, 'r') as tree_json:
        tree = json.load(tree_json)

    with open(meta_path, 'r') as meta_json:
        meta = json.load(meta_json)

    t = Tree(newick_tree, format=1)

    seqs = SeqIO.to_dict(SeqIO.parse(seq_path, "fasta"))

    df = pd.read_csv(dataframe)

    with open(config_path, 'r') as config_json:
        config = json.load(config_json)

    sites = ast.literal_eval(sites)
    for site in sites:
        site_gt = 'ha'+str(site)

        genotype_nodes = {}
        for n in t.traverse(strategy='postorder'):
            site_genotype = seqs[n.name][int(site)-1]
            genotype_nodes[n.name] = site_genotype

        def add_genotypes(branch, site):
            if site_gt not in branch['attr'].keys():
                genotype = genotype_nodes[branch['strain']]

                branch['attr'][site_gt] = genotype

                if 'children' in branch.keys():
                    for child in branch['children']:
                        add_genotypes(child, site)

        add_genotypes(tree, site)

        if manual_colors == True:
            color_map = [['K', '#60AA9E'],
                         ['R', '#5097BA'],
                         ['T', '#D9AD3D'],
                         ['A', '#E67030'],
                         ['I', '#8EBC66']]

        else:
            site_genotypes = list(df[str(160)].unique())
            auspice_colors = ["#60AA9E", "#D9AD3D", "#5097BA", "#E67030", "#8EBC66", "#E59637", "#AABD52", "#DF4327", "#C4B945", "#75B681"]
            color_map = [[gt, auspice_colors[site_genotypes.index(gt)]] for gt in site_genotypes]

        if site_gt not in config['color_options'].keys():
            config['color_options'][site_gt] = {"color_map": color_map,
                                                "menuItem": site_gt,
                                                "legendTitle": site_gt,
                                                "type": "discrete",
                                                "key": site_gt}
        if site_gt not in meta['color_options'].keys():
            meta['color_options'][site_gt] = {"color_map": color_map,
                                            "menuItem": site_gt,
                                            "legendTitle": site_gt,
                                            "type": "discrete",
                                            "key": site_gt}

    with open(tree_path, 'w') as outfile:
        json.dump(tree, outfile, indent=2, sort_keys = True)

    with open(config_path, 'w') as config_out:
        json.dump(config, config_out, indent=2, sort_keys = True)

    with open(meta_path, 'w') as meta_out:
        json.dump(meta, meta_out, indent=2, sort_keys = True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Work around to recolor specific genotypes for auspice",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', help= "path to _tree.json file")
    parser.add_argument('--meta', help= "path to _meta.json file")
    parser.add_argument('--seqs', help= "path to aa-seq_ fasta file")
    parser.add_argument('--newick', help= "path to newick tree")
    parser.add_argument('--dataframe', help= "path to dataframe")
    parser.add_argument('--config', help= "path to auspice config file")
    parser.add_argument('--sites', default = '[160]', help="list of HA1 positions to manually color")
    parser.add_argument('--manual_colors', default = True, help="hard code colors or automatically assign from list of auspice colors")

    args = parser.parse_args()

    main(tree_path = args.tree, meta_path = args.meta, seq_path = args.seqs, newick_tree = args.newick, dataframe = args.dataframe, config_path = args.config, sites = args.sites, manual_colors = args.manual_colors)
