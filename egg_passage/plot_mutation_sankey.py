import os
import argparse
import plotly.plotly as py
import plotly.io as pio
import pandas as pd
import matplotlib.pyplot as plt

def plot_sankey(prefix):
    df = pd.read_csv('data/'+prefix+'_egg_df.csv')
    positions = [col[3:] for col in df.columns if col[0:3]=='mut']


    def sankey(site):
        before_after = []

        group = df.groupby(['circulating'+site, site]).size().reset_index(name='count')
        before_aas = list(group['circulating'+site].unique())
        after_aas = list(group[site].unique())

        colors = ['rgba(73, 148, 206, 1)', 'rgba(242, 116, 32, 1)', 'rgba(127, 194, 65, 1)',
                  'rgba(211, 211, 211, 1)', 'rgba(250, 188, 19, 1)',
                  'rgba(138, 89, 136,1)', 'rgba(68, 158, 158,1)', 'rgba(158, 68, 68,1)']
        colors_transparent = ['rgba(73, 148, 206, 0.25)', 'rgba(242, 116, 32, 0.25)', 'rgba(127, 194, 65, 0.25)',
                'rgba(211, 211, 211, 0.25)', 'rgba(250, 188, 19, 0.25)',
                  'rgba(138, 89, 136,0.25)', 'rgba(68, 158, 158,0.25)', 'rgba(158, 68, 68,0.25)']
        colors_node = ['rgba(73, 148, 206, 0.5)', 'rgba(242, 116, 32, 0.5)', 'rgba(127, 194, 65, 0.5)',
            'rgba(211, 211, 211, 0.5)', 'rgba(250, 188, 19, 0.5)',
              'rgba(138, 89, 136,0.5)', 'rgba(68, 158, 158,0.5)', 'rgba(158, 68, 68,0.5)']
        aa_colors = {}
        aa_colors_transparent = {}
        aa_colors_node = {}
        for aa in after_aas:
            aa_colors[aa] = colors[after_aas.index(aa)]
            aa_colors_transparent[aa] = colors_transparent[after_aas.index(aa)]
            aa_colors_node[aa] = colors_node[after_aas.index(aa)]

        before_id = {}
        after_id = {}
        labels = []

        for i in range(len(before_aas)):
            before_id[before_aas[i]] = i
            labels.append(before_aas[i])

        for i in range(len(after_aas)):
            after_id[after_aas[i]] = i + len(before_aas)
            labels.append(after_aas[i])

        label_colors = [aa_colors_node[label] for label in labels]

        for k, v in group.iterrows():
            if v['circulating'+site]== v[site]:
                link_color = aa_colors_transparent[v[site]]
            else:
                link_color = aa_colors[v[site]]
            before_after.append({'site': site, 'before': v['circulating'+site], 'source': before_id[v['circulating'+site]],
                                 'after': v[site], 'target': after_id[v[site]], 'count': v['count'],
                                 'link_color': link_color})

        before_after_df = pd.DataFrame(before_after)

        sankey_trace = dict(
            type='sankey',
            textfont=dict(
                size=28
            ),
            node = dict(
                pad = 20,
                thickness = 200,
                line = dict(
                    width = 0
                ),
                label = labels,
                color = label_colors
            ),
            link = dict(
                source =  before_after_df['source'],
                target =  before_after_df['target'],
                value =  before_after_df['count'],
                color = before_after_df['link_color']
            ))

        layout =  dict(
            annotations = [dict(
                x=0.49,
                y=1.2,
                font = dict(size = 24),
                showarrow=False,
                text= 'HA1 site '+str(site)),
                           dict(
                x=0.0,
                y=1.2,
                font = dict(size = 14),
                showarrow=False,
                text= 'Genotype<br>Before egg-passaging'),
                           dict(
                x=0.99,
                y=1.2,
                font = dict(size = 14),
                showarrow=False,
                text= 'Genotype<br>After egg-passaging')],
        )

        fig = dict(data=[sankey_trace], layout=layout)
        pio.write_image(fig, 'plots/mutation_sankey/mutations_'+str(site)+'.pdf')

    for position in positions:
        sankey(position)
    #
    # fig, ax = plt.subplots(4,2)
    # for position in positions:
    #     mut_sankey = plt.imread('plots/mutation_sankey/mutations_'+str(position)+'.jpeg')
    #     if positions.index(position) <= 3:
    #         ax[positions.index(position), 0].imshow(mut_sankey)
    #         ax[positions.index(position), 0].axis('off')
    #     else:
    #         ax[(positions.index(position)-4), 1].imshow(mut_sankey)
    #         ax[(positions.index(position)-4), 1].axis('off')
    # fig.savefig('plots/mutations_sankey_'+str(prefix)+'.pdf', bbox_inches='tight')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Plots sankey flow diagram of genotypes at each position before and after egg passaging")
    parser.add_argument('--prefix', default= 'h3n2_6y_hi', help= "specify prefix for naming data files")
    args = parser.parse_args()

    plot_sankey(prefix = args.prefix)
