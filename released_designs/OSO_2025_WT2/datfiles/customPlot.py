# This cell does plotting things that you don't need to use.  Skip to the next cell

import matplotlib.pyplot as plt
import copy

def custom_plot(gcf,fname,*args,**kwargs):
    if 'legend_dictionary' in kwargs:
        legend_dict = kwargs.pop('legend_dictionary')
        if 'facecolor' not in legend_dict:
            legend_dict['facecolor'] = 'black'
        if 'labelcolor' not in legend_dict:
            legend_dict['labelcolor'] = 'white'
        if 'edgecolor' not in legend_dict:
            legend_dict['edgecolor'] = 'white'
    else:
        legend_dict = {'facecolor':'black','labelcolor':'white', 'edgecolor':'white'}

    # Check for one period in the file name
    assert fname.count('.')==1, 'File name must contain exactly one period (.) to separate the file name and extension.'

    fname_write = fname.split('.')[0]
    ext = fname.split('.')[-1]
    gcf.savefig(f'{fname_write}.svg', *args, **kwargs)
    gcf.savefig(f'{fname_write}.pdf', *args, **kwargs)
    gcf.savefig(f'{fname_write}.'+ext, *args, **kwargs)

    fig2 = copy.deepcopy(gcf)

    # Mismatch found in line 15:axes.edgecolor
    # White: black
    # Black: white
    # Mismatch found in line 16:axes.facecolor
    # White: white
    # Black: black
    # Mismatch found in line 26:axes.labelcolor
    # White: black
    # Black: white
    # Mismatch found in line 57:boxplot.boxprops.color
    # White: black
    # Black: white
    # Mismatch found in line 60:boxplot.capprops.color
    # White: black
    # Black: white
    # Mismatch found in line 63:boxplot.flierprops.color
    # White: black
    # Black: white
    # Mismatch found in line 67:boxplot.flierprops.markeredgecolor
    # White: black
    # Black: white
    # Mismatch found in line 89:boxplot.whiskerprops.color
    # White: black
    # Black: white
    # Mismatch found in line 116:figure.edgecolor
    # White: white
    # Black: black
    # Mismatch found in line 117:figure.facecolor
    # White: white
    # Black: black
    # Mismatch found in line 145:grid.color
    # White: #b0b0b0
    # Black: white
    # Mismatch found in line 196:lines.color
    # White: C0
    # Black: white
    # Mismatch found in line 225:patch.edgecolor
    # White: black
    # Black: white

    fig2.set_edgecolor('white')
    fig2.set_facecolor('black')

    for i, ax in enumerate(fig2.axes):
        ax.set_facecolor('black')
        if ax.xaxis.label.get_color() == 'black':
            ax.xaxis.label.set_color('white')
        if ax.yaxis.label.get_color() == 'black':
            ax.yaxis.label.set_color('white')

        if 'color' in ax.xaxis.get_tick_params():
            if ax.xaxis.get_tick_params()['color'] == 'black':
                ax.tick_params(axis='x', which='both', colors='white')
        else:
            ax.tick_params(axis='x', which='both', colors='white')

        if 'color' in ax.yaxis.get_tick_params():
            if ax.yaxis.get_tick_params()['color'] == 'black':
                ax.tick_params(axis='y', which='both', colors='white')
        else:
            ax.tick_params(axis='y', which='both', colors='white')

        # print(ax.yaxis.get_tick_params())
        # ax.tick_params(axis='both', which='both', colors='white')
        # print(ax.yaxis.get_tick_params())

        for spine in ax.spines.values():
            if spine.get_edgecolor()[0] == 0.0 and spine.get_edgecolor()[1] == 0.0 and spine.get_edgecolor()[2] == 0.0:
                spine.set_edgecolor('white')
        ax.title.set_color('white')
        
        # # ==================================================
        # # Done with AI 
        # # ==================================================
        # # Handle hatching colors for patches (like fill_between)
        # for patch in ax.patches:
        #     # If the patch has hatching and a black edge color, make it white
        #     if hasattr(patch, 'get_hatch') and patch.get_hatch() is not None:
        #         edge_color = patch.get_edgecolor()
        #         # Check if edge color is black (or very dark)
        #         if hasattr(edge_color, '__len__') and len(edge_color) >= 3:
        #             if edge_color[0] < 0.1 and edge_color[1] < 0.1 and edge_color[2] < 0.1:
        #                 patch.set_edgecolor('white')
        #         elif edge_color == 'black' or edge_color == 'k':
        #             patch.set_edgecolor('white')
        # # ==================================================

        if ax.get_legend() is not None:
            ax.legend(**legend_dict)

    fig2.savefig(f'{fname_write}_dark.svg', *args, **kwargs)
    fig2.savefig(f'{fname_write}_dark.pdf', *args, **kwargs)
    fig2.savefig(f'{fname_write}_dark.'+ext, *args, **kwargs)
    plt.close(fig2)