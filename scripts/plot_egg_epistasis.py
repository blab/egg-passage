"""
Find epistasis between egg-specific mutations
"""

import argparse, itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.io as pio
import plotly.plotly as py
import plotly.graph_objs as go

def plot_heatmap(prefix, limit_clades):
    """
    Plot heatmaps indicating epistatic interactions between genotypes. Heat maps indicate whether the genotype of certain HA sites are correlated via a log enrichment ratio, that compares the frequency of observing genotype 1 at site 1 AND genotype 2 at site 2 versus the expected frequency (based on overall prevalence of the genotypes)
    """

    df = pd.read_csv('dataframes/'+prefix+'_egg.csv')

    #Analyze only viruses from clades where greater than 10% of egg-seqs have mutation at 186 AND 194
    if limit_clades==True:
        clade_pct_186 = (df[df['mut186']==True].groupby('clade').size()/
                         df.groupby('clade').size()).reset_index().rename(columns={0:'pct_186'})
        clade_pct_194 = (df[df['mut194']==True].groupby('clade').size()/
                         df.groupby('clade').size()).reset_index().rename(columns={0:'pct_194'})
        clade_pcts = pd.merge(clade_pct_186, clade_pct_194, on='clade').fillna(0.0)
        limited_clades = []
        for k, v in clade_pcts.iterrows():
            if v['pct_186']>0.1:
                if v['pct_194']>0.1:
                    limited_clades.append(v.clade)
        df = df[df['clade'].isin(limited_clades)]

    mut_sites = [col for col in df.columns if col[0:3]=='mut']

    #Exclude sites 225 and 246 because no epistatic interaction with 186 or 194
    if limit_clades==True:
        for x in ['mut138', 'mut225', 'mut246']:
            if x in mut_sites:
                mut_sites.remove(x)

    sites = [str(x[3:]) for x in mut_sites]

    aa_epistasis = []
    aa_alone = []

    for site1 in sites:
        other_sites = [x for x in sites if site1 not in x]
        #Group egg sequences by their genotype at site1
        site1_counts = df.groupby(site1).size().reset_index().rename(columns={0:'count'})

        for k, v in site1_counts.iterrows():
            #Find the prevalence of each genotype at site1
            overall_proportion = float(v['count'])/float(len(df))
            aa_alone.append({'site': site1, 'aa': v[site1], 'overall_proportion': overall_proportion})

        for site2 in other_sites:
            #Group egg sequences by site1 and all other sites
            both_counts = df.groupby([site1, site2]).size().reset_index().rename(columns={0:'count'})

            for i, r in both_counts.iterrows():
                #Find prevalence of each site1 genotype co-occuring with each site2 genotype
                proportion = float(r['count'])/float(len(df))
                aa_epistasis.append({'site1': site1, 'site1_aa': r[site1], 'site2': site2, 'site2_aa': r[site2], 'count': r['count'], 'proportion': proportion})


    aa_epistasis_df = pd.DataFrame(aa_epistasis)
    aa_alone_df = pd.DataFrame(aa_alone)

    #set threshold to exclude groups with only a few sequences
    #note: this means proportions won't exactly add to 1- doesn't matter for comparison to exp.
    aa_epistasis_df = aa_epistasis_df[aa_epistasis_df['count'] >= 5]

    fig, ax = plt.subplots(len(sites), len(sites), sharex='col', sharey='row')

    #Make list of lists to store arrays of enrichment ratios, like matrix
    #Row1, column1 stores interactions with sites[1], etc
    plot_arrays = [[0 for x in sites] for y in sites]

    #Add matrix with number of amino acid genotypes as dimensions
    for site_x in sites:
        for site_y in sites:
            # plot_arrays[sites.index(site_x)][sites.index(site_y)] = np.zeros(shape = (len(list(aa_epistasis_df[aa_epistasis_df['site1']==site_x]['site1_aa'].unique())), len(list(aa_epistasis_df[aa_epistasis_df['site2']==site_y]['site2_aa'].unique()))))
            plot_arrays[sites.index(site_x)][sites.index(site_y)] = np.empty(shape = (len(list(aa_epistasis_df[aa_epistasis_df['site1']==site_x]['site1_aa'].unique())), len(list(aa_epistasis_df[aa_epistasis_df['site2']==site_y]['site2_aa'].unique()))))
            #Fill all with -2.0, so if observed count is 0, log enrichment will be -2.0
            plot_arrays[sites.index(site_x)][sites.index(site_y)].fill(-2.0)


    #Fill in plot_arrays with enrichment ratios
    for k, v in aa_epistasis_df.iterrows():
        obs = v['proportion']
        exp = float(aa_alone_df[(aa_alone_df['site'] == v['site1']) & (aa_alone_df['aa'] == v['site1_aa'])]['overall_proportion']) * float(aa_alone_df[(aa_alone_df['site'] == v['site2']) & (aa_alone_df['aa'] == v['site2_aa'])]['overall_proportion'])
        enrichment = float(obs)/exp
        log_enrichment = np.log2(enrichment)
        site1_list = list(aa_epistasis_df[aa_epistasis_df['site1']==v['site1']]['site1_aa'].unique())
        site2_list = list(aa_epistasis_df[aa_epistasis_df['site2']==v['site2']]['site2_aa'].unique())
        plot_arrays[sites.index(v['site1'])][sites.index(v['site2'])][site1_list.index(v['site1_aa'])][site2_list.index(v['site2_aa'])]= log_enrichment

    #Manually enter predominant egg mutation genotypes
    egg_muts = {'160':'K', '194': 'P', '186':'V', '225':'G', '219':['F','Y'], '203':'I', '156':['R','Q'], '138':'S', '246':['H'], '183':['L']}

    for site1 in range(len(sites)):
        for site2 in range(len(sites)):
            site1_list = list(aa_epistasis_df[aa_epistasis_df['site1']==str(sites[site1])]['site1_aa'].unique())
            site2_list = list(aa_epistasis_df[aa_epistasis_df['site2']==str(sites[site2])]['site2_aa'].unique())
            cmap = sns.diverging_palette(220, 20, sep=10, as_cmap=True)
            heatmap = ax[site1, site2].imshow(plot_arrays[site1][site2], cmap=cmap, vmin= -2.0, vmax=2.0, aspect='auto')

            ax[site1, site2].tick_params(axis='both', which='both', length=0)
            ax[site1, site2].set_xticks([p for p in range(len(site2_list))])
            ax[site1, site2].set_yticks([p for p in range(len(site1_list))])
            ax[site1, site2].tick_params(labelbottom=False)
            ax[site1, site2].set_xticks(np.arange(len(site2_list)+1)-.5, minor=True)
            ax[site1, site2].set_yticks(np.arange(len(site1_list)+1)-.5, minor=True)
            ax[site1, site2].grid(which='minor', color='white', linestyle='-', linewidth=1)

            if site1 > (site2-1):
                ax[site1, site2].set_visible(False)

            if site1 == 0:
                ax[site1, site2].xaxis.set_label_position('top')
                ax[site1, site2].xaxis.set_ticks_position('top')
                ax[site1, site2].set(xlabel = str(sites[site2]))
                ax[site1, site2].set_xticklabels([str(s) for s in site2_list])
                colors= ['red' if str(s) in egg_muts[str(sites[site2])] else 'black' for s in site2_list]
                for xtick, color in zip(ax[site1, site2].get_xticklabels(), colors):
                    xtick.set_color(color)

            if site2 == (site1+1):
                ax[site1, site2].set(ylabel = str(sites[site1]))
                ax[site1, site2].yaxis.set_ticks_position('left')
                ax[site1, site2].set_yticklabels([str(s) for s in site1_list])
                colors= ['red' if str(s) in egg_muts[str(sites[site1])] else 'black' for s in site1_list]
                for ytick, color in zip(ax[site1, site2].get_yticklabels(), colors):
                    ytick.set_color(color)

    cbar_ax = fig.add_axes([0.95, 0.2, 0.05, 0.7])
    colorbar = fig.colorbar(heatmap, cax=cbar_ax)
    colorbar.ax.get_yaxis().labelpad = 15
    colorbar.ax.set_ylabel('log2 enrichment ratio', rotation=270)

    if limit_clades==False:
        fig.suptitle('Epistasis between HA sites in egg-passaged influenza H3N2', fontsize=12, y=1.05, x=0.6)
        fig.savefig('plots/'+str(prefix)+'/epistasis_heatmap_'+str(prefix)+'.pdf',
                    bbox_inches='tight')
    elif limit_clades==True:
        fig.suptitle('Epistasis between HA sites in egg-passaged H3N2 viruses \n(from clades without a 186 vs. 194 mutation preference)', fontsize=12, y=1.05, x=0.6)
        fig.savefig('plots/'+str(prefix)+'/epistasis_heatmap_'+str(prefix)+'_limit_clades.pdf',
                    bbox_inches='tight')

def plot_chord(prefix, limit_clades):
    #Adapted from https://nbviewer.jupyter.org/github/empet/Plotly-plots/blob/master/Chord-diagram.ipynb?flush_cache=true

    df = pd.read_csv('dataframes/'+prefix+'_egg.csv')

    #Analyze only viruses from clades where greater than 10% of egg-seqs have mutation at 186 AND 194
    if limit_clades==True:
        clade_pct_186 = (df[df['mut186']==True].groupby('clade').size()/
                         df.groupby('clade').size()).reset_index().rename(columns={0:'pct_186'})
        clade_pct_194 = (df[df['mut194']==True].groupby('clade').size()/
                         df.groupby('clade').size()).reset_index().rename(columns={0:'pct_194'})
        clade_pcts = pd.merge(clade_pct_186, clade_pct_194, on='clade').fillna(0.0)
        limited_clades = []
        for k, v in clade_pcts.iterrows():
            if v['pct_186']>0.1:
                if v['pct_194']>0.1:
                    limited_clades.append(v.clade)
        df = df[df['clade'].isin(limited_clades)]

    sites = ['186','194','138','156','203','219','225','246']
    for site in sites:
        if site not in df.columns:
            sites.remove(site)
    #Exclude sites 225 and 246 because no epistatic interaction with 186 or 194
    if limit_clades==True:
        for x in ['138', '225', '246']:
            if x in sites:
                sites.remove(x)


    egg_muts = {'186':['V'], '225':['G'], '219':['F','Y'],
                '203':['I'], '156':['R','Q'], '138':['S'],
                '246':['H'], '183':['L'], '194': ['P'], '246':['H']}
    matrix = np.zeros((len(sites), len(sites)), dtype=int)
    for site_a in sites:
        for site_b in sites:
            matrix_entry = 0
            #Set diagonal to be number of viruses with ONLY that mutation
            if site_a == site_b:
                only_sitea = 0
                other_sites = sites.copy()
                other_sites.remove(site_a)
                for egg_genotype_a in egg_muts[site_a]:
                    mask = df[site_a] == egg_genotype_a
                    for other_site in other_sites:
                        for other_genotype in egg_muts[other_site]:
                            res = df[other_site] != other_genotype
                            mask &= res
                    only_sitea += len(df[mask])
                matrix_entry= only_sitea
            #Find number of viruses with both egg mutations
            else:
                for egg_genotype_a in egg_muts[site_a]:
                    for egg_genotype_b in egg_muts[site_b]:
                        matrix_entry+= len(df[(df[site_a]==egg_genotype_a)&(df[site_b]==egg_genotype_b)])
            matrix[sites.index(site_a)][sites.index(site_b)] = matrix_entry

    labels = []
    for site in sites:
        site_label = str(site)+egg_muts[site][0]
        for x in range(1,len(egg_muts[site])):
            site_label+='/'+str(egg_muts[site][x])

        labels.append(site_label)

    ideo_colors = ['rgba(26,152,80, 0.75)',
                   'rgba(215,48,39, 0.75)',
                 'rgba(244,109,67, 0.75)',
                 'rgba(253,174,97, 0.75)',
                 'rgba(254,224,139, 0.75)',
                 'rgba(217,239,139, 0.75)',
                 'rgba(166,217,106, 0.75)',
                 'rgba(102,189,99, 0.75)']

    def check_data(data_matrix):
        L, M=data_matrix.shape
        if L!=M:
            raise ValueError('Data array must have (n,n) shape')
        return L

    L=check_data(matrix)

    radii_sribb=[0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]

    def moduloAB(x, a, b): #maps a real number onto the unit circle identified with
                           #the interval [a,b), b-a=2*PI
            if a>= b:
                raise ValueError('Incorrect interval ends')
            y = (x-a) % (b-a)
            return y+b if y < 0 else y+a

    def test_2PI(x):
        return True
    #      return 0 <= x < 2*pi

    pi = np.pi

    row_sum = [np.sum(matrix[k,:]) for k in range(L)]

    #set the gap between two consecutive ideograms
    gap = 2*pi*0.005
    ideogram_length = 2*pi * np.asarray(row_sum) / sum(row_sum) - gap*np.ones(L)

    def get_ideogram_ends(ideogram_len, gap):
        ideo_ends = []
        left = 0
        for k in range(len(ideogram_len)):
            right = left + ideogram_len[k]
            ideo_ends.append([left, right])
            left = right + gap
        return ideo_ends

    ideo_ends = get_ideogram_ends(ideogram_length, gap)

    def make_ideogram_arc(R, phi, a=100):
        # R is the circle radius
        # phi is the list of  angle coordinates of an arc ends
        # a is a parameter that controls the number of points to be evaluated on an arc
        if not test_2PI(phi[0]) or not test_2PI(phi[1]):
            phi = [moduloAB(t, 0, 2*pi) for t in phi]
        length = (phi[1]-phi[0]) % 2*pi
        nr = 18 if length <= pi/4 else int(a*length/pi)

        if phi[0] < phi[1]:
            theta = np.linspace(phi[0], phi[1], nr)
        else:
            phi = [moduloAB(t, -pi, pi) for t in phi]
            theta = np.linspace(phi[0], phi[1], nr)
        return R * np.exp(1j*theta)

    make_ideogram_arc(1.3, [11*pi/6, pi/17])

    def map_data(data_matrix, row_value, ideogram_length):
        mapped = np.zeros(data_matrix.shape)
        for j  in range(L):
            mapped[:, j] = ideogram_length * data_matrix[:,j] / row_value
        return mapped

    mapped_data = map_data(matrix, row_sum, ideogram_length)
    mapped_data

    idx_sort = np.argsort(mapped_data, axis=1)
    idx_sort

    def make_ribbon_ends(mapped_data, ideo_ends,  idx_sort):
        L = mapped_data.shape[0]
        ribbon_boundary = np.zeros((L,L+1))
        for k in range(L):
            start = ideo_ends[k][0]
            ribbon_boundary[k][0] = start
            for j in range(1,L+1):
                J = idx_sort[k][j-1]
                ribbon_boundary[k][j] = start + mapped_data[k][J]
                start = ribbon_boundary[k][j]
        return [[(ribbon_boundary[k][j], ribbon_boundary[k][j+1] ) for j in range(L)] for k in range(L)]

    ribbon_ends = make_ribbon_ends(mapped_data, ideo_ends,  idx_sort)

    def control_pts(angle, radius):
        #angle is a  3-list containing angular coordinates of the control points b0, b1, b2
        #radius is the distance from b1 to the  origin O(0,0)

        if len(angle) != 3:
            raise InvalidInputError('angle must have len =3')
        b_cplx = np.array([np.exp(1j*angle[k]) for k in range(3)])
        b_cplx[1] = radius * b_cplx[1]
        return list(zip(b_cplx.real, b_cplx.imag))

    def ctrl_rib_chords(l, r, radius):
        # this function returns a 2-list containing control poligons of the two quadratic Bezier
        #curves that are opposite sides in a ribbon
        #l (r) the list of angular variables of the ribbon arc ends defining
        #the ribbon starting (ending) arc
        # radius is a common parameter for both control polygons
        if len(l) != 2 or len(r) != 2:
            raise ValueError('the arc ends must be elements in a list of len 2')
        return [control_pts([l[j], (l[j]+r[j])/2, r[j]], radius) for j in range(2)]

    #Define ribbon color
    ribbon_color = [L * ['rgba(175,175,175,0.5)'] for k in range(L)]

    #Change 186 and 194 ribbon colors
    for k_186 in range(len(ribbon_color[0])):
        ribbon_color[0][k_186] = ideo_colors[0]
    for k_194 in range(len(ribbon_color[1])):
        ribbon_color[1][k_194] = ideo_colors[1]


    def make_q_bezier(b):# defines the Plotly SVG path for a quadratic Bezier curve defined by the
                         #list of its control points
        if len(b) != 3:
            raise valueError('control poligon must have 3 points')
        A, B, C = b
        return f'M {A[0]}, {A[1]} Q {B[0]}, {B[1]} {C[0]}, {C[1]}'

    b=[(1,4), (-0.5, 2.35), (3.745, 1.47)]
    make_q_bezier(b)

    def make_ribbon_arc(theta0, theta1):

        if test_2PI(theta0) and test_2PI(theta1):
            if theta0 < theta1:
                theta0 = moduloAB(theta0, -pi, pi)
                theta1 = moduloAB(theta1, -pi, pi)
                if theta0  *theta1 > 0:
                    raise ValueError('incorrect angle coordinates for ribbon')

            nr = int(40 * (theta0 - theta1) / pi)
            if nr <= 2: nr = 3
            theta = np.linspace(theta0, theta1, nr)
            pts=np.exp(1j*theta)# points in polar complex form, on the given arc

            string_arc = ''
            for k in range(len(theta)):
                string_arc += f'L {pts.real[k]}, {pts.imag[k]} '
            return   string_arc
        else:
            raise ValueError('the angle coordinates for an arc side of a ribbon must be in [0, 2*pi]')

    make_ribbon_arc(np.pi/3, np.pi/6)

    def make_layout(title, plot_size):

        return dict(title=title,
                    xaxis=dict(visible=False),
                    yaxis=dict(visible=False),
                    showlegend=False,
                    width=plot_size,
                    height=plot_size,
                    margin=dict(t=25, b=25, l=25, r=25),
                    hovermode=False,
                     )

    def make_ideo_shape(path, line_color, fill_color):
        #line_color is the color of the shape boundary
        #fill_color is the color assigned to an ideogram

        return  dict(line=dict(color=line_color,
                               width=0.45),
                     path=path,
                     layer='below',
                     type='path',
                     fillcolor=fill_color)

    def make_ribbon(l, r, line_color, fill_color, radius=0.2):
        #l=[l[0], l[1]], r=[r[0], r[1]]  represent the opposite arcs in the ribbon
        #line_color is the color of the shape boundary
        #fill_color is the fill color for the ribbon shape

        poligon = ctrl_rib_chords(l,r, radius)
        b, c = poligon

        return  dict(line=dict(color=line_color,
                                 width=0.5),
                     path=make_q_bezier(b) + make_ribbon_arc(r[0], r[1])+
                             make_q_bezier(c[::-1]) + make_ribbon_arc(l[1], l[0]),
                     type='path',
                     layer='below',
                     fillcolor = fill_color,
            )

    def make_self_rel(l, line_color, fill_color, radius):
        #radius is the radius of Bezier control point b_1

        b = control_pts([l[0], (l[0]+l[1])/2, l[1]], radius)

        return  dict(line = dict(color=line_color,
                                 width=0.5),
                     path =  make_q_bezier(b)+make_ribbon_arc(l[1], l[0]),
                     type = 'path',
                     layer = 'below',
                     fillcolor = fill_color
                    )

    def invPerm(perm):
        # function that returns the inverse of a permutation, perm
        inv = [0] * len(perm)
        for i, s in enumerate(perm):
            inv[s] = i
        return inv

    if limit_clades==False:
        layout=make_layout('Pairwise interactions between common egg mutation genotypes', 800)
    elif limit_clades==True:
        layout=make_layout('Pairwise interactions between common egg mutation genotypes <br> (in clades without a 186 vs. 194 mutation preference)', 800)

    ribbon_info=[]
    shapes=[]
    annotations = []
    for k in range(L):

        sigma = idx_sort[k]
        sigma_inv = invPerm(sigma)
        for j in range(k, L):
            if matrix[k][j] == 0 and matrix[j][k]==0: continue
            eta = idx_sort[j]
            eta_inv = invPerm(eta)
            l = ribbon_ends[k][sigma_inv[j]]

            if j == k:
                shapes.append(make_self_rel(l, 'rgba(175,175,175,0.5)',
                                            'rgba(175,175,175,0.5)', radius=radii_sribb[k]))
    #             shapes.append(make_self_rel(l, ideo_colors[k] ,
    #                         ideo_colors[k], radius=radii_sribb[k]))

                z = 0.9*np.exp(1j*(l[0]+l[1])/2)


                ribbon_info.append(go.Scatter(x=[z.real],
                                           y=[z.imag],
                                           mode='markers',
                                           marker=dict(size=0.5, color=ideo_colors[k])
                                           )
                                  )
            else:
                r = ribbon_ends[j][eta_inv[k]]
                zi = 0.9 * np.exp(1j*(l[0]+l[1])/2)
                zf = 0.9 * np.exp(1j*(r[0]+r[1])/2)
                #texti and textf are the strings that will be displayed when hovering the mouse
                #over the two ribbon ends
                texti = f'{labels[k]}'
                textf = f'{labels[j]}'

                ribbon_info.append(go.Scatter(x=[zi.real],
                                              y=[zi.imag],
                                              mode='markers',
                                              marker=dict(size=0.5, color=ribbon_color[k][j]),
                                              text=texti
                                           )
                                  ),
                ribbon_info.append(go.Scatter(x=[zf.real],
                                              y=[zf.imag],
                                              mode='markers',
                                              marker=dict(size=0.5, color=ribbon_color[k][j]),
                                              text=textf
                                           )
                                  )
                r = (r[1], r[0]) # IMPORTANT!!!  Reverse these arc ends because otherwise you get
                              # a twisted ribbon
                #append the ribbon shape
                shapes.append(make_ribbon(l, r, ribbon_color[k][j] , ribbon_color[k][j]))



    ideograms = []
    for k in range(len(ideo_ends)):
        z =  make_ideogram_arc(1.1, ideo_ends[k])
        zi = make_ideogram_arc(1.0, ideo_ends[k])
        m = len(z)
        n = len(zi)
        ideograms.append(go.Scatter(x=z.real,
                                    y=z.imag,
                                    mode='lines',
                                    line=dict(color=ideo_colors[k], shape='spline', width=0.25)
    #                                 text=f'{labels[k]} <br>{int(row_sum[k])} viruses'
                                 )
                         )


        path = 'M '
        for s in range(m):
            path += f'{z.real[s]}, {z.imag[s]} L '

        Zi = np.array(zi.tolist()[::-1])

        for s in range(m):
            path += f'{Zi.real[s]}, {Zi.imag[s]} L '
        path += f'{z.real[0]} ,{z.imag[0]}'

        shapes.append(make_ideo_shape(path,ideo_colors[k] , ideo_colors[k]))

        z_text =  make_ideogram_arc(1.2, ideo_ends[k])
        Z_text = np.array(z_text.tolist()[::-1])
        annotations.append(dict(x=Z_text.real[int(len(Z_text)/2)],
                                y=Z_text.imag[int(len(Z_text)/2)],
                                showarrow=False,
                                align='left',
                                valign='bottom',
                                text=f'{labels[k]}'))


    data = ideograms + ribbon_info
    layout['shapes'] = shapes
    layout['annotations'] = annotations
    fig = go.Figure(data=data, layout=layout)
    if limit_clades==False:
        pio.write_image(fig, 'plots/'+str(prefix)+'/epistasis_chord_diagram_'+str(prefix)+'.pdf')
    elif limit_clades==True:
        pio.write_image(fig, 'plots/'+str(prefix)+'/epistasis_chord_diagram_'+str(prefix)+'_limit_clades.pdf')

def main(input_df):
    df_name = str.split(input_df, 'dataframes/')[1]
    prefix = str.split(df_name, '.csv')[0]
    plot_heatmap(prefix, limit_clades=False)
    plot_chord(prefix, limit_clades=False)
    if '6y' in prefix:
        plot_heatmap(prefix, limit_clades=True)
        plot_chord(prefix, limit_clades=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Determines epistasis between egg-specific mutations")
    parser.add_argument('--in_file', help= "input dataframe file")
    args = parser.parse_args()

    main(input_df = args.in_file)
