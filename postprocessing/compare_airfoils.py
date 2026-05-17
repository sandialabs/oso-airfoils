import numpy as np
import pandas as pd
from kulfan import Kulfan
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})
import matplotlib
colors = ['#0065cc', '#eea800', '#009e73', '#d55e00', '#7860aa', '#56b4ff', '#fca7c7', '#ede13f', '#5d5d5d', '#000000']
matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=colors)

import numpy as np 
import matplotlib
import matplotlib.pyplot as plt
import math


from xfoil_wrapper_noprint import run as run_xfoil
from neuralfoil_wrapper_noprint import run as run_neuralfoil
from kulfan import Kulfan

import matplotlib.pyplot as plt
import numpy as np
import json
import os
import numpy as np


def handleZeroDivide(num,dem):
    if dem == 0:
        return np.inf
    else:
        return num/dem
    
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def get_colors(N, selectedMap, lower = 0.15, upper = 0.75):
    cmap = plt.get_cmap(selectedMap)
    cmap = truncate_colormap(cmap,lower, upper)
    norm = plt.Normalize(0, N - 1)
    colors = cmap(norm(np.arange(N)))
    return colors

def get_fractional_color(frac, selectedMap, lower=0.15, upper=0.75):
    colors = get_colors(101, selectedMap, lower, upper)
    return(colors[int(math.floor(frac*100.0))])

def computeNormals(xdata, ydata):
    xmin_index = xdata.index(min(xdata))
    normals = []
    for i, x in enumerate(xdata):
        y = ydata[i]

        if i == 0:
            x_p1 = xdata[i+1]
            x_m1 = x
            y_p1 = ydata[i+1]
            y_m1 = y

            dxdy = handleZeroDivide(y_p1-y, x_p1-x)

        elif i == len(xdata)-1:
            x_p1 = x
            x_m1 = xdata[i-1]
            y_p1 = y
            y_m1 = ydata[i-1]

            dxdy = handleZeroDivide(y-y_m1, x-x_m1)

        else:
            x_p1 = xdata[i+1]
            x_m1 = xdata[i-1]
            y_p1 = ydata[i+1]
            y_m1 = ydata[i-1]

            dxdy_p1 = handleZeroDivide(y_p1-y, x_p1-x)
            dxdy_m1 = handleZeroDivide(y-y_m1, x-x_m1)

            if dxdy_p1!=dxdy_m1 and max([dxdy_p1,dxdy_m1])==np.inf:
                dxdy = min([dxdy_p1,dxdy_m1])
            else:
                dxdy = handleZeroDivide(y_p1-y_m1, x_p1-x_m1)

        if dxdy == 0.0:
            normal = [0,1]
        elif dxdy == np.inf:
            normal = [1,0]
        else:
            normal =  [ 1/(1+1/dxdy**2)**0.5, -1/dxdy/(1+1/dxdy**2)**0.5 ]

        if i<=xmin_index:
            # in upper
            if normal[1] < 0:
                normal = [-1*normal[0], -1*normal[1]]
        else:
            if x>1:
                #in wake
                if normal[1] < 0:
                    normal = [-1*normal[0], -1*normal[1]]
            else:
                #in lower
                if normal[1] > 0:
                    normal = [-1*normal[0], -1*normal[1]]

        if i < xmin_index+3 and i > xmin_index-3:
            # leading edge points need to point to left, but sometimes there are issues
            if normal[0] > 0:
                normal = [-1*normal[0], -1*normal[1]]

        normals.append(normal)

    return normals

def _arc_length_indices(x_vals, y_vals, n_marks):
    """Return a list of indices evenly spaced by arc length along (x_vals, y_vals)."""
    xrange = max(x_vals.max() - x_vals.min(), 1e-9)
    yrange = max(y_vals.max() - y_vals.min(), 1e-9)
    dx_n = np.diff(x_vals) / xrange
    dy_n = np.diff(y_vals) / yrange
    arc = np.concatenate([[0], np.cumsum(np.sqrt(dx_n**2 + dy_n**2))])
    targets = np.linspace(0, arc[-1], n_marks)
    return [int(np.argmin(np.abs(arc - t))) for t in targets]

def polarPlot(dataList, airfoil_coords=None, legend_entries=None,
              linecolors=None, linestyles=None, 
              coordinatecolors=None, coordinatestyles=None,
              width_ratios=None, wspace=0.12):
    """Plot polar data.

    Parameters
    ----------
    dataList : list of lists of dicts
        Polar data, one list per dataset.
    airfoil_coords : list of array-like, optional
        One entry per airfoil; each entry is a list/array of [x, y] pairs.
        Plotted in the upper-right panel.
    legend_entries : list of lists of dicts, optional
        Outer list = columns, inner list = rows within that column.
        Each dict may contain: 'text', 'linestyle', 'linecolor', 'markersize'.
        A markersize of 0 or None draws only the line.
        Displayed in the upper-left panel.
    linecolors : list of colors, optional
        Per-dataset line color.
    linestyles : list of str, optional
        Per-dataset line style.
    width_ratios : list of 5 floats, optional
        Relative widths of the 5 polar panels [cpmin, drag, CL/CM, L/D, xtr].
    wspace : float
        Horizontal spacing between the 5 polar panels.
    """
    
    for data in dataList:
        assert(isinstance(data, list))
        for de in data:
            assert(isinstance(de, dict))

    # colorCycle = ['Blues','Oranges','Greens','Purples','Reds']
    # assert(len(dataList) <= len(colorCycle))

    dataframeList = []
    for data in dataList:
        dataDict = {}
        dataDict['alpha']    = []
        dataDict['cl']       = []
        dataDict['cd']       = []
        dataDict['cm']       = []
        dataDict['cpmin']    = []
        dataDict['xtp_u']    = []
        dataDict['xtp_l']    = []
        dataDict['xtr_u']    = []
        dataDict['xtr_l']    = []
        dataDict['re']       = []
        dataDict['m']        = []
        dataDict['n_crit']   = []
        dataDict['n_panels'] = []

        for i, rdata in enumerate(data):
            if rdata is not None:
                for ky in dataDict.keys():
                    if ky in list(rdata.keys()):
                        vl = rdata[ky]
                    elif ky == 're':
                        vl = rdata['Re']
                    elif ky == 'm':
                        vl = rdata['M']
                    elif ky == 'n_crit':
                        vl = rdata['N_crit']
                    elif ky == 'xtr_u':
                        vl = rdata['xtr_top']
                    elif ky == 'xtr_l':
                        vl = rdata['xtr_bot']
                    elif ky == 'xtp_u':
                        vl = rdata['xtp_top']
                    elif ky == 'xtp_l':
                        vl = rdata['xtp_bot']
                    elif ky == 'n_panels':
                        vl = rdata['N_panels']
                    else:
                        raise ValueError('Could not find key: %s' % (ky))
                    dataDict[ky].append(vl)

        assert(len(np.unique(dataDict['m'])) == 1)
        assert(len(np.unique(dataDict['xtp_u'])) == 1)
        assert(len(np.unique(dataDict['xtp_l'])) == 1)
        assert(len(np.unique(dataDict['n_crit'])) == 1)

        df = pd.DataFrame.from_dict(dataDict)
        dataframeList.append(df)

    if width_ratios is None:
        width_ratios = [0.15, 0.2, 0.2, 0.25, 0.2]
    gs_width_ratios = [w / sum(width_ratios) * 5 for w in width_ratios]

    # ---- Figure and GridSpec layout ----
    show_top_row = (airfoil_coords is not None) or (legend_entries is not None)

    fig_width = 25
    bot_row_height_in = 7.0

    if show_top_row:
        if airfoil_coords is not None:
            # Compute bounding box across all provided airfoil coordinate sets
            all_coords = np.vstack([np.array(c) for c in airfoil_coords])
            x_range = max(all_coords[:, 0].max() - all_coords[:, 0].min(), 1e-9)
            y_range = max(all_coords[:, 1].max() - all_coords[:, 1].min(), 1e-9)
            # Airfoil panel spans cols 2–4; estimate physical width in inches
            # (0.80 accounts for wspace, constrained-layout margins, etc.)
            airfoil_col_frac = sum(width_ratios[2:]) / sum(width_ratios)
            airfoil_panel_w_in = fig_width * airfoil_col_frac * 0.80
            # Height required for equal aspect, plus room for axis labels
            top_row_height_in = max(airfoil_panel_w_in * (y_range / x_range) + 1.2, 3.0)
        else:
            top_row_height_in = 3.5

        fig_height = top_row_height_in + bot_row_height_in
        fig = plt.figure(figsize=(fig_width, fig_height), dpi=300, layout='constrained')
        outer_gs = GridSpec(2, 1, figure=fig, height_ratios=[top_row_height_in, bot_row_height_in])
        top_gs   = GridSpecFromSubplotSpec(1, 5, subplot_spec=outer_gs[0], width_ratios=gs_width_ratios, wspace=wspace)
        ax_legend  = fig.add_subplot(top_gs[0, 0:2])
        ax_airfoil = fig.add_subplot(top_gs[0, 2:5])
        bot_gs = GridSpecFromSubplotSpec(1, 5, subplot_spec=outer_gs[1],
                                         width_ratios=gs_width_ratios, wspace=wspace)
    else:
        fig = plt.figure(figsize=(fig_width, bot_row_height_in), dpi=300, layout='constrained')
        bot_gs = GridSpec(1, 5, figure=fig, width_ratios=gs_width_ratios, wspace=wspace)

    ax0  = fig.add_subplot(bot_gs[0, 0])
    ax1  = fig.add_subplot(bot_gs[0, 1])
    ax2  = fig.add_subplot(bot_gs[0, 2])
    ax2r = ax2.twinx()
    ax3  = fig.add_subplot(bot_gs[0, 3])
    ax4  = fig.add_subplot(bot_gs[0, 4])

    # ---- Top row: legend and airfoil shape ----
    if show_top_row:
        # Legend panel
        ax_legend.axis('off')
        if legend_entries is not None:
            # legend_entries is a list of columns; each column is a list of entry dicts.
            n_cols = len(legend_entries)
            for col_idx, col_entries in enumerate(legend_entries):
                handles = []
                for entry in col_entries:
                    ms = entry.get('markersize', None)
                    if ms is None or ms == 0:
                        marker = 'None'
                        ms_val = 0
                    else:
                        marker = 'o'
                        ms_val = ms
                    h = Line2D([0], [0],
                               color=entry.get('linecolor', 'black'),
                               linestyle=entry.get('linestyle', '-'),
                               marker=marker,
                               markersize=ms_val,
                               label=entry.get('text', ''))
                    handles.append(h)
                x_anchor = col_idx / max(n_cols, 1) + 0.02
                leg = ax_legend.legend(handles=handles, loc='upper left',
                                       frameon=False,
                                       bbox_to_anchor=(x_anchor, 0.95))
                if col_idx < n_cols - 1:
                    ax_legend.add_artist(leg)

        # Airfoil shape panel
        if airfoil_coords is not None:
            for j, coords in enumerate(airfoil_coords):
                coords_arr = np.array(coords)
                ac  = coordinatecolors[j]  if (coordinatecolors  is not None and j < len(coordinatecolors))  else colors[j % len(colors)]
                als = coordinatestyles[j]  if (coordinatestyles  is not None and j < len(coordinatestyles))  else '-'
                ax_airfoil.plot(coords_arr[:, 0], coords_arr[:, 1], color=ac, linestyle=als)
            ax_airfoil.set_aspect('equal')
            ax_airfoil.grid(True, linewidth=0.5, alpha=0.5)
            ax_airfoil.minorticks_on()
            ax_airfoil.grid(True, which='minor', linewidth=0.4, alpha=0.4)
            ax_airfoil.set_xlabel('$x/c$')
            ax_airfoil.set_ylabel('$y/c$')
        else:
            ax_airfoil.axis('off')

    clmin = np.inf
    clmax = -np.inf
    cmmin = np.inf
    cmmax = -np.inf
    # ---- Bottom row: polar plots ----
    for i, df in enumerate(dataframeList):
        if linecolors is not None and i < len(linecolors):
            plot_color = linecolors[i]
        else:
            plot_color = colors[i % len(colors)]

        if linestyles is not None and i < len(linestyles):
            linestyle = linestyles[i]
        else:
            linestyle = '-'

        dataEntry = df.sort_values('alpha')

        N_marks = 8

        clmin = min(clmin, min(dataEntry['cl']))
        clmax = max(clmax, max(dataEntry['cl']))
        cmmin = min(cmmin, min(dataEntry['cm']))
        cmmax = max(cmmax, max(dataEntry['cm']))

        # Plot 0: Cp_min vs CL
        ax0.plot(dataEntry['cpmin'], dataEntry['cl'], color=plot_color, linestyle=linestyle)
        ax0.set_xlabel('$C_{p,min}$')
        ax0.set_ylabel('$C_L$', labelpad=0)
        ax0.xaxis.set_inverted(True)
        ax0.grid(1)

        # Plot 1: drag polar
        ax1.plot(dataEntry['cd'], dataEntry['cl'], color=plot_color, linestyle=linestyle)
        ax1.set_xlim([0, 0.05])
        ax1.set_ylabel('$C_L$', labelpad=0)
        ax1.set_xlabel('$C_D$')
        ax1.grid(1)

        # Plot 2 (left axis): CL vs alpha
        ax2.plot(dataEntry['alpha'], dataEntry['cl'], color=plot_color, linestyle=linestyle)
        ax2.set_ylabel('$C_L$', labelpad=0)
        ax2.set_xlabel(r'$\alpha$')
        ax2.grid(1)

        # Plot 2 (right axis): CM vs alpha — arc-length-spaced open markers
        alpha_vals  = dataEntry['alpha'].values
        cm_vals     = dataEntry['cm'].values
        mark_idx_cm = _arc_length_indices(alpha_vals, cm_vals, N_marks)
        ax2r.plot(alpha_vals, cm_vals, 'o', linestyle=linestyle, color=plot_color,
                  markersize=3, markevery=mark_idx_cm)
        ax2r.set_ylabel('$C_M$', labelpad=+3)

        # Plot 3: L/D vs CL
        ax3.plot(dataEntry['cl'] / dataEntry['cd'], dataEntry['cl'],
                 color=plot_color, linestyle=linestyle)
        ax3.set_xlabel('$L/D$')
        ax3.set_ylabel('$C_L$', labelpad=0)
        ax3.grid(1)

        # Plot 4: transition location vs CL
        #   upper surface — solid line; lower surface — solid line + arc-length-spaced markers
        ax4.plot(dataEntry['xtr_u'], dataEntry['cl'], linestyle=linestyle, color=plot_color, markersize=3)
        xtr_l_vals  = dataEntry['xtr_l'].values
        cl_vals_arr = dataEntry['cl'].values
        mark_idx_l  = _arc_length_indices(xtr_l_vals, cl_vals_arr, N_marks)
        ax4.plot(xtr_l_vals, cl_vals_arr, 'o', linestyle=linestyle, color=plot_color,
                 markersize=3, markevery=mark_idx_l)
        ax4.set_xlabel('$x_{tr}/c$')
        ax4.set_ylabel('$C_L$', labelpad=0)
        ax4.set_xlim([0, 1.01])
        ax4.grid(1)

    ax0.set_xlim(left=0.0)

    # Align CM right axis to CL left axis so gridlines coincide
    cl_lo, cl_hi = ax2.get_ylim()
    mid_cm = 0.5 * (cmmin + cmmax)
    # round midcm to the nearest 0.04
    mid_cm_rounded = round(mid_cm / 0.04) * 0.04
    zero_cmplot = 2/3 * clmin + 1/3 * clmax
    zero_cmplot_rounded = round(zero_cmplot / 0.2) * 0.2 
    # land the cm_mid_rounded at the same gridline as zero_cmplot_rounded, but not too close to the edge of the plot


    ax2r.set_ylim((cl_lo - zero_cmplot_rounded + 0.2*mid_cm_rounded/0.04) / 5, (cl_hi - zero_cmplot_rounded + 0.2*mid_cm_rounded/0.04) / 5)

    # Major tick spacing
    for _ax in [ax0, ax1, ax2, ax3, ax4]:
        _ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax2r.yaxis.set_major_locator(MultipleLocator(0.04))

    # Minor gridlines
    for _ax in [ax0, ax1, ax2, ax2r, ax3, ax4]:
        _ax.minorticks_on()
        _ax.grid(True, which='minor', linewidth=0.4, alpha=0.4)

    # Legend on ax2/ax2r panel (attached to ax2r so framealpha works over twinx lines)
    ax2r.legend(handles=[
        Line2D([0], [0], color='black', linestyle='-', label='$C_L$'),
        Line2D([0], [0], color='black', linestyle='-', marker='o', markersize=3, label='$C_M$'),
    ], loc='upper left', framealpha=1.0)

    # Legend on ax4 panel
    ax4.legend(handles=[
        Line2D([0], [0], color='black', linestyle='-', label='Upper'),
        Line2D([0], [0], color='black', linestyle='-', marker='o', markersize=3, label='Lower'),
    ], loc='lower right', framealpha=1.0)

    return fig


def expand_polar_result(res):
    """Expand a single neuralfoil result dict (array-valued fields) into a list of
    per-alpha scalar dicts suitable for polarPlot."""
    scalar_keys = ['Re', 'M', 'N_crit', 'xtp_top', 'xtp_bot', 'N_panels']
    array_keys  = ['alpha', 'cl', 'cd', 'cm', 'cpmin', 'xtr_top', 'xtr_bot']
    entries = []
    for j in range(len(res['alpha'])):
        entry = {k: res[k] for k in scalar_keys}
        for k in array_keys:
            entry[k] = res[k][j]
        entries.append(entry)
    return entries


def compare_airfoils(afl_dict_input, reynolds_numbers, turb_cases, tools, figurePath, color_override=None, alpha_min=-5, alpha_max=20, alpha_step=0.25):
    """Run aerodynamic analyses and produce a polar comparison plot.

    Parameters
    ----------
    afl_dict_input : dict
        Mapping of airfoil name -> airfoil definition. Each value may be:
          - str ending in '.dat'          : path to a coordinate file
          - dict with keys 'K_upper', 'K_lower', 'TE_gap'  : Kulfan coefficients
          - Kulfan instance               : used directly
          - list/tuple of 2 lists         : [x_coords, y_coords], >= 50 points each
          - flat list of 2*N_k+1 numbers  : [K_upper..., K_lower..., TE_gap]
    reynolds_numbers : list of float
        Reynolds numbers to evaluate (e.g. [1.5e6, 3e6]).
        If more than one Re and more than 5 airfoils, an error is raised.
    turb_cases : list of [N_crit, xtp_u, xtp_l]
        One or two turbulence cases. Each entry is a 3-element list:
          [N_crit (float), xtp_upper (0-1), xtp_lower (0-1)]
        e.g. [[9, 1.0, 1.0], [3, 0.05, 0.05]]
    tools : list of str
        Solvers to use. Valid entries: 'neuralfoil', 'xfoil'.
        e.g. ['neuralfoil'] or ['neuralfoil', 'xfoil']
    color_override : dict, optional
        Mapping of airfoil name -> matplotlib color string.
        When a single Reynolds number is used, the named airfoil's lines and
        legend entry will use this color instead of the default palette.
        e.g. {'mhkf1-180': '#000000'}
    alpha_min : float, optional
        Minimum angle of attack (degrees). Default is -5.
    alpha_max : float, optional
        Maximum angle of attack (degrees). Default is 20.
    alpha_step : float, optional
        Step size for angle of attack (degrees). Default is 1.0.

    Returns
    -------
    matplotlib.figure.Figure
        The completed polar comparison figure. Also saved to 'polar_comparison.png'.
    """

    if color_override is None:
        color_override = {}

    if len(afl_dict_input) > 5 and len(reynolds_numbers) > 1:
        raise ValueError('Too many airfoils and Reynolds numbers may result in an overcrowded plot. Please reduce the number of airfoils or Reynolds numbers.')

    afl_dict = {}

    for name,afl_input in afl_dict_input.items():
        if isinstance(afl_input, str) and afl_input.endswith('.dat'):
            afl = Kulfan()
            afl.readFile(afl_input)

        elif isinstance(afl_input, dict) and 'K_upper' in afl_input and 'K_lower' in afl_input and 'TE_gap' in afl_input:
            afl = Kulfan(TE_gap=afl_input.get('TE_gap'))
            afl.upperCoefficients = afl_input['K_upper']
            afl.lowerCoefficients = afl_input['K_lower']
        
        elif isinstance(afl_input, Kulfan):
            afl = afl_input

        elif isinstance(afl_input, (list, tuple)):
            if len(afl_input) == 2:
                assert isinstance(afl_input[0], list) and isinstance(afl_input[1], list)
                assert len(afl_input[0]) == len(afl_input[1])
                assert(all(isinstance(x, (int, float)) for x in afl_input[0]))
                assert(all(isinstance(x, (int, float)) for x in afl_input[1]))
                if len(afl_input[0]) >= 50:
                    afl = Kulfan()
                    afl.fit2coordinates(afl_input[0], afl_input[1])
                else:
                    raise ValueError('Coordinate lists must have at least 50 points for fitting.')
            else:
                assert(len(afl_input) <= 20)
                assert(len(afl_input) % 2 == 1)
                assert(all(isinstance(x, (int, float)) for x in afl_input))
                N_k = (len(afl_input) - 1) // 2
                afl = Kulfan(TE_gap=afl_input[-1])
                afl.upperCoefficients = afl_input[0:N_k]
                afl.lowerCoefficients = afl_input[N_k:2*N_k]
        else:
            print(afl_input.keys() if isinstance(afl_input, dict) else str(type(afl_input)))
            raise ValueError('Invalid airfoil input format: %s'%(str(afl_input)))
        
        afl_dict[name] = afl


    re_min = min(reynolds_numbers)
    re_max = max(reynolds_numbers)

    colorCycle = ['Blues','Oranges','Greens','Purples','Reds']

    data_list = []
    coords_list = []
    leg = []
    linecolors = []
    linestyles = []
    coordinatecolors = []

    # airfoil names and tools
    legend_col_1 = []

    # Reynolds Numbers
    legend_col_2 = []

    # Turbulence Cases
    legend_col_3 = []

    re_list_legend = []

    is_rainbow_plot = False
    for i, (name, afl) in enumerate(afl_dict.items()):
        for re in reynolds_numbers:
            if re_min == re_max:
                if len(afl_dict) > 8:
                    # assume this is a rainbow plot
                    turbo_cmap = plt.get_cmap('turbo', len(afl_dict))
                    color = turbo_cmap(i)
                    is_rainbow_plot = True
                else:
                    colors = ['#0065cc', '#e69f00', '#009e73', '#d55e00', '#56b4ff', '#fca7c7', '#ede13f', '#666666', '#000000']
                    color = colors[i]
                if re not in re_list_legend:
                    legend_col_2.append({'text': 'Re=%.2e'%(re), 'linestyle': '-', 'linecolor': 'k', 'markersize': 0})
                    re_list_legend.append(re)
            else:
                color_frac = (np.log10(re) - np.log10(re_min)) / (np.log10(re_max)-np.log10(re_min))
                color = get_fractional_color(color_frac, colorCycle[i], lower=0.3, upper=0.8)
                grey_color = get_fractional_color(color_frac, 'Greys', lower=0.3, upper=0.8)
                if re not in re_list_legend:
                    legend_col_2.append({'text': 'Re=%.2e'%(re), 'linestyle': '-', 'linecolor': grey_color, 'markersize': 0})
                    re_list_legend.append(re)

            if re_min == re_max and name in color_override:
                color = color_override[name]
            linecolor = color

            if len(coordinatecolors) <= i and re_min != re_max:
                ccolor = get_fractional_color(1.0, colorCycle[i], lower=0.6, upper=0.6)
                coordinatecolors.append(ccolor)
                legend_entry = {'text': name, 'linestyle': '-', 'linecolor': ccolor, 'markersize': 0}
                legend_col_1.append(legend_entry)

                coords_list.append(list(zip(afl.xcoordinates, afl.ycoordinates)))

            elif re_min == re_max:
                legend_entry = {'text': name, 'linestyle': '-', 'linecolor': color, 'markersize': 0}
                legend_col_1.append(legend_entry)

                coordinatecolors.append(color)
                coords_list.append(list(zip(afl.xcoordinates, afl.ycoordinates)))
            else:
                # already handled
                pass
        
            # legend_entry = {'text': name, 'linestyle': '-', 'linecolor': color, 'markersize': 0}
            
            K_upper = afl.upperCoefficients
            K_lower = afl.lowerCoefficients
            # N_crit_clean = 9
            # xtp_u_clean = 1.0
            # xtp_l_clean = 1.0

            if turb_cases is not None:
                assert(len(turb_cases) == 1 or len(turb_cases) == 2)
            else:
                turb_cases = [[9.0, 1.0, 1.0]]

            for ii, turb_case in enumerate(turb_cases):
                assert(len(turb_case) == 3)
                N_crit = turb_case[0]
                xtp_u = turb_case[1]
                xtp_l = turb_case[2]        

                for tool in tools:
                    if tool == 'xfoil':
                        res = run_xfoil('alfa', K_upper, K_lower, [alpha_min, alpha_max, alpha_step], Re=re, N_crit=N_crit, xtp_u=xtp_u, xtp_l=xtp_l, TE_gap = afl.constants.TE_gap)
                        data_list.append(expand_polar_result(res))
                        if ii == 0:
                            linestyles.append('-')
                        else:
                            linestyles.append('--')
                        linecolors.append(linecolor)
                        # legend_entry['linecolor'] = linecolor
                    elif tool == 'neuralfoil':
                        res = run_neuralfoil('alfa', K_upper, K_lower, [alpha_min, alpha_max, alpha_step], Re=re, N_crit=N_crit, xtp_u=xtp_u, xtp_l=xtp_l, TE_gap = afl.constants.TE_gap)
                        data_list.append(expand_polar_result(res))
                        if 'xfoil' in tools:
                            if ii == 0:
                                linestyles.append('-.')
                            else:
                                linestyles.append(':')
                        else:
                            if ii == 0:
                                linestyles.append('-')
                            else:
                                linestyles.append('--')
                        linecolors.append(linecolor)
                        # legend_entry['linecolor'] = linecolor
                    else:
                        raise ValueError('Unknown tool: %s'%(tool))



    if 'xfoil' in tools and 'neuralfoil' in tools:
        legend_col_1.append( {'text': 'NeuralFoil', 'linestyle': '-.', 'linecolor': 'k', 'markersize': 0} )
    elif 'neuralfoil' in tools:
        legend_col_1.append( {'text': 'NeuralFoil', 'linestyle': '-', 'linecolor': 'k', 'markersize': 0} )
    elif 'xfoil' in tools:
        legend_col_1.append( {'text': 'XFOIL', 'linestyle': '-', 'linecolor': 'k', 'markersize': 0} )
    else:
        raise ValueError('Unknown tool combination in legend construction.')


    if tools == ['xfoil', 'neuralfoil'] or tools == ['neuralfoil', 'xfoil']:
        if len(turb_cases) == 1:
            legend_col_3.append( {'text': r'XF--$N_{crit}$: %.1f, $x_{tp,u}$: %.2f, $x_{tp,l}$: %.2f'%(turb_cases[0][0], turb_cases[0][1], turb_cases[0][2]), 'linestyle': '-', 'linecolor': 'k', 'markersize': 0} )
            # legend_col_3.append( {'text': r'XF--$N_{crit}$: %.1f, $x_{tp,u}$: %.2f, $x_{tp,l}$: %.2f'%(turb_cases[1][0], turb_cases[1][1], turb_cases[1][2]), 'linestyle': '--', 'linecolor': 'k', 'markersize': 0} )
            legend_col_3.append( {'text': r'NF--$N_{crit}$: %.1f, $x_{tp,u}$: %.2f, $x_{tp,l}$: %.2f'%(turb_cases[0][0], turb_cases[0][1], turb_cases[0][2]), 'linestyle': '-.', 'linecolor': 'k', 'markersize': 0} )
            # legend_col_3.append( {'text': r'NF--$N_{crit}$: %.1f, $x_{tp,u}$: %.2f, $x_{tp,l}$: %.2f'%(turb_cases[1][0], turb_cases[1][1], turb_cases[1][2]), 'linestyle': ':', 'linecolor': 'k', 'markersize': 0} )
        else:
            assert(len(turb_cases) == 2)
            legend_col_3.append( {'text': r'XF--$N_{crit}$: %.1f, $x_{tp,u}$: %.2f, $x_{tp,l}$: %.2f'%(turb_cases[0][0], turb_cases[0][1], turb_cases[0][2]), 'linestyle': '-', 'linecolor': 'k', 'markersize': 0} )
            legend_col_3.append( {'text': r'XF--$N_{crit}$: %.1f, $x_{tp,u}$: %.2f, $x_{tp,l}$: %.2f'%(turb_cases[1][0], turb_cases[1][1], turb_cases[1][2]), 'linestyle': '--', 'linecolor': 'k', 'markersize': 0} )
            legend_col_3.append( {'text': r'NF--$N_{crit}$: %.1f, $x_{tp,u}$: %.2f, $x_{tp,l}$: %.2f'%(turb_cases[0][0], turb_cases[0][1], turb_cases[0][2]), 'linestyle': '-.', 'linecolor': 'k', 'markersize': 0} )
            legend_col_3.append( {'text': r'NF--$N_{crit}$: %.1f, $x_{tp,u}$: %.2f, $x_{tp,l}$: %.2f'%(turb_cases[1][0], turb_cases[1][1], turb_cases[1][2]), 'linestyle': ':', 'linecolor': 'k', 'markersize': 0} )

    elif tools == ['xfoil']:
        if len(turb_cases) == 1:
            legend_col_3.append( {'text': r'$N_{crit}$: %.1f, $x_{tp,u}$: %.2f, $x_{tp,l}$: %.2f'%(turb_cases[0][0], turb_cases[0][1], turb_cases[0][2]), 'linestyle': '-', 'linecolor': 'k', 'markersize': 0} )
        else:
            assert(len(turb_cases) == 2)
            legend_col_3.append( {'text': r'$N_{crit}$: %.1f, $x_{tp,u}$: %.2f, $x_{tp,l}$: %.2f'%(turb_cases[0][0], turb_cases[0][1], turb_cases[0][2]), 'linestyle': '-', 'linecolor': 'k', 'markersize': 0} )
            legend_col_3.append( {'text': r'$N_{crit}$: %.1f, $x_{tp,u}$: %.2f, $x_{tp,l}$: %.2f'%(turb_cases[1][0], turb_cases[1][1], turb_cases[1][2]), 'linestyle': '--', 'linecolor': 'k', 'markersize': 0} )

    elif tools == ['neuralfoil']:
        if len(turb_cases) == 1:
            legend_col_3.append( {'text': r'$N_{crit}$: %.1f, $x_{tp,u}$: %.2f, $x_{tp,l}$: %.2f'%(turb_cases[0][0], turb_cases[0][1], turb_cases[0][2]), 'linestyle': '-', 'linecolor': 'k', 'markersize': 0} )
        else:
            assert(len(turb_cases) == 2)
            legend_col_3.append( {'text': r'$N_{crit}$: %.1f, $x_{tp,u}$: %.2f, $x_{tp,l}$: %.2f'%(turb_cases[0][0], turb_cases[0][1], turb_cases[0][2]), 'linestyle': '-', 'linecolor': 'k', 'markersize': 0} )
            legend_col_3.append( {'text': r'$N_{crit}$: %.1f, $x_{tp,u}$: %.2f, $x_{tp,l}$: %.2f'%(turb_cases[1][0], turb_cases[1][1], turb_cases[1][2]), 'linestyle': '--', 'linecolor': 'k', 'markersize': 0} )
    else:
        raise ValueError('Unknown tool combination in legend construction.')

    if is_rainbow_plot:
        leg = [[],[],[]]
        if len(tools) == 2:
            tool_entries = legend_col_1[-2:]
        else:
            tool_entries = legend_col_1[-1:]
        # tool_entry = legend_col_1[-1]
        for i in range(len(legend_col_1)-1):
            entry = legend_col_1[i]
            leg[i%3].append(entry)
        for te in tool_entries:
            leg[0].append(te)
        re_entry = legend_col_2[-1]
        leg[1].append(re_entry)
        for entry in legend_col_3:
            leg[2].append(entry)
    else:
        leg = [legend_col_1, legend_col_2, legend_col_3] 

    fig = polarPlot(data_list,
            airfoil_coords=coords_list,
            legend_entries=leg,
            linecolors=linecolors,
            linestyles=linestyles,
            coordinatecolors=coordinatecolors,
    )
    fig.savefig(figurePath, dpi=300)
    # return fig



def rainbow_plot(path_to_data, figurePath = None, comparison_airfoil = None, color_override = None, tools = None):
    if figurePath is None:
        data_folder = os.path.dirname(path_to_data)
        figurePath = os.path.join(data_folder, 'rainbow_plot.png')

    nafl = 21

    data = json.load(open(path_to_data, 'r'))
    pop = data['population']
    pareto_points = [p for p in pop if p['pareto_index'] == 1]
    pareto_points = sorted(pareto_points, key=lambda p: p['LoD_rough_at_design'])
    rough_LD_vals = np.array([p['LoD_rough_at_design'] for p in pareto_points])
    clean_LD_vals = np.array([p['LoD_clean_at_design'] for p in pareto_points])

    x_norm = float(np.ptp(rough_LD_vals)) or 1.0
    y_norm = float(np.ptp(clean_LD_vals)) or 1.0
    rn = rough_LD_vals / x_norm
    cn = clean_LD_vals / y_norm
    seg_lens = np.sqrt(np.diff(rn) ** 2 + np.diff(cn) ** 2)
    cum_arc = np.concatenate(([0.0], np.cumsum(seg_lens)))
    total_arc = cum_arc[-1]
    if total_arc <= 0 or len(rough_LD_vals) < 2:
        ixs = np.unique(np.linspace(0, max(len(rough_LD_vals) - 1, 0), nafl).astype(int))
    else:
        target_arcs = np.linspace(0.0, total_arc, nafl)
        ixs = np.unique(np.array([int(np.argmin(np.abs(cum_arc - t))) for t in target_arcs]))
    afl_dict_input = {'Rough. Factor=%.2f'%(j / (len(ixs)-1)): pareto_points[k] for j, k in enumerate(ixs)}
    for k in afl_dict_input.keys():
        afl_dict_input[k]['TE_gap'] = data['input_parameters']['TE_gap']

    if comparison_airfoil is not None:
        for ky,vl in comparison_airfoil.items():
            afl_dict_input[ky] = vl

    filename = os.path.basename(path_to_data)
    tokens = filename.split('_')
    for tk in tokens:
        if tk.startswith('e'):
            re = float(tk[1:])*10**5

    reynolds_numbers = [re]
    turb_cases = [[9, 1.0, 1.0], [3, 0.05, 0.05]]
    if tools is None:
        tools = ['neuralfoil']
    
    compare_airfoils(afl_dict_input, reynolds_numbers, turb_cases, tools, figurePath, color_override=color_override)



if __name__ == '__main__':

    import pathlib
    path_to_here = pathlib.Path(__file__).parent.resolve()
    path_to_oso = path_to_here.parent
    path_to_datfiles = path_to_oso / 'historical_airfoils/mhkf1/'

    afl_dict_input = {
        # 'mhkf1-180'  : 'mhkf1-180.dat',
        'DU-96-W-180' : str( path_to_oso / 'historical_airfoils/du/du_96-w-180.dat' ),
        'FFA-W1-182'  : str( path_to_oso / 'historical_airfoils/ffa/fitted/FFA-W1-182_fittedCST10.dat' ),
        'RISO-B-17'   : str( path_to_oso / 'historical_airfoils/riso-b/riso-b-17.dat' ),
        'S831'        : str( path_to_oso / 'historical_airfoils/s/s831.dat' ),
        }

    reynolds_numbers = [3e6, 15e6]
    turb_cases = [[9, 1.0, 1.0], [3, 0.05, 0.05]]
    tools = ['neuralfoil']
    figurePath = 'polar_comparison.png'
    compare_airfoils(afl_dict_input, reynolds_numbers, turb_cases, tools, figurePath)