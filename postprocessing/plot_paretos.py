import json, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import pathlib
path_to_here = pathlib.Path(__file__).parent.resolve()
path_to_oso = path_to_here.parent

colors = ['#0065cc', '#eea800', '#009e73', '#d55e00', '#7860aa', '#ede13f', '#56b4ff', '#fca7c7', '#5d5d5d', '#000000']

plot_dict = {
    'data' : [ 
        # 
        {'data_location':[str(path_to_oso / "postprocessing/cases/cases_101_to_110/case_110/c110_t21_l15_k16_g2000_n752_x3_s3__2026_02_26_09-47"),
                          str(path_to_oso / 'postprocessing/cases/cases_101_to_110/case_109/c109_t21_l15_k16_g2000_n752__2026_01_07_12-47')],                                     
                          'color': colors[0], 'label': r'$\tau=0.21$ Unconstrained', 'linestyle': '-', 'linewidth': 2.0, 'marker': None},
        {'data_location':[str(path_to_oso / 'postprocessing/cases/cases_101_to_110/case_110/c110_t24_l14_k16_g2000_n752_x3_s1__2026_03_03_16-50'),
                        #   str(path_to_oso / 'postprocessing/cases/cases_101_to_110/case_110/c110_t24_l14_k16_g2000_n752_m14_p14__2026_01_22_11-18'
                          ], 
                          'color': colors[1], 'label': r'$\tau=0.24$ Unconstrained', 'linestyle': '-', 'linewidth': 2.0, 'marker': None},
        {'data_location':[str(path_to_oso / 'postprocessing/cases/cases_101_to_110/case_110/c110_t27_l13_k16_g2000_n752_x3_s2__2026_02_26_09-49'),
                          str(path_to_oso / 'postprocessing/cases/cases_101_to_110/case_110/c110_t27_l13_k16_g2000_n752_m14_p14__2026_02_21_12-43')], 
                          'color': colors[2], 'label': r'$\tau=0.27$ Unconstrained', 'linestyle': '-', 'linewidth': 2.0, 'marker': None},
        {'data_location':[str(path_to_oso / 'postprocessing/cases/cases_101_to_110/case_110/c110_t30_l12_k16_g2000_n752_x3_s3__2026_02_26_09-54'),
                          str(path_to_oso / 'postprocessing/cases/cases_101_to_110/case_110/c110_t30_l12_k16_g2000_n752_m14_p14__2026_02_27_13-34')], 
                          'color': colors[3], 'label': r'$\tau=0.30$ Unconstrained', 'linestyle': '-', 'linewidth': 2.0, 'marker': None},
        {'data_location':[str(path_to_oso / 'postprocessing/cases/cases_101_to_110/case_110/c110_t33_l12_k16_g2000_n752_x3_s1__2026_02_26_10-01'),
                          str(path_to_oso / 'postprocessing/cases/cases_101_to_110/case_110/c110_t33_l12_k16_g2000_n752_m14_p14__2026_03_05_10-37')], 
                          'color': colors[4], 'label': r'$\tau=0.33$ Unconstrained', 'linestyle': '-', 'linewidth': 2.0, 'marker': None},
        {'data_location':[str(path_to_oso / 'postprocessing/cases/cases_101_to_110/case_110/c110_t36_l12_k16_g2000_n752_x3_s2__2026_02_26_10-12'),
                          str(path_to_oso / 'postprocessing/cases/cases_111_to_120/case_111/c111_t36_l12_k16_g2000_n752_m14_p14__2026_03_12_23-29')], 
                          'color': colors[5], 'label': r'$\tau=0.36$ Unconstrained', 'linestyle': '-', 'linewidth': 2.0, 'marker': None},
        #
        {'data_location':str(path_to_oso / 'postprocessing/cases/cases_101_to_110/case_110/c110_t21_l15_k16_g2000_n752_m14_p14__2026_01_17_08-45'), 
                             'color': colors[0], 'label': r'$\tau=0.21$ Moment Constrained', 'linestyle': '-.', 'linewidth':2.0, 'marker': None},
        {'data_location':str(path_to_oso / 'postprocessing/cases/cases_101_to_110/case_110/c110_t24_l14_k16_g2000_n752_m14_p14__2026_01_22_11-18'), 
                             'color': colors[1], 'label': r'$\tau=0.24$ Moment Constrained', 'linestyle': '-.', 'linewidth':2.0, 'marker': None},
        {'data_location':str(path_to_oso / 'postprocessing/cases/cases_101_to_110/case_110/c110_t27_l13_k16_g2000_n752_m14_p14__2026_02_21_12-43'), 
                             'color': colors[2], 'label': r'$\tau=0.27$ Moment Constrained', 'linestyle': '-.', 'linewidth':2.0, 'marker': None},
        {'data_location':str(path_to_oso / 'postprocessing/cases/cases_101_to_110/case_110/c110_t30_l12_k16_g2000_n752_m14_p14__2026_02_27_13-34'), 
                             'color': colors[3], 'label': r'$\tau=0.30$ Moment Constrained', 'linestyle': '-.', 'linewidth':2.0, 'marker': None},
        {'data_location':str(path_to_oso / 'postprocessing/cases/cases_101_to_110/case_110/c110_t33_l12_k16_g2000_n752_m14_p14__2026_03_05_10-37'), 
                             'color': colors[4], 'label': r'$\tau=0.33$ Moment Constrained', 'linestyle': '-.', 'linewidth':2.0, 'marker': None},
        {'data_location':str(path_to_oso / 'postprocessing/cases/cases_111_to_120/case_111/c111_t36_l12_k16_g2000_n752_m14_p14__2026_03_12_23-29'), 
                             'color': colors[5], 'label': r'$\tau=0.36$ Moment Constrained', 'linestyle': '-.', 'linewidth':2.0, 'marker': None},
        #
        {'data':(248.2, 120.2), 'color': colors[0], 'label': r'WT2 $\tau=0.21$ RFOIL', 'linestyle': None, 'marker': 'o', 'markersize': 8},
        {'data':(183.9, 118.1), 'color': colors[1], 'label': r'WT2 $\tau=0.24$ RFOIL', 'linestyle': None, 'marker': 'o', 'markersize': 8},
        {'data':(172.5, 112.2), 'color': colors[2], 'label': r'WT2 $\tau=0.27$ RFOIL', 'linestyle': None, 'marker': 'o', 'markersize': 8},
        {'data':(167.4, 97.4),  'color': colors[3], 'label': r'WT2 $\tau=0.30$ RFOIL', 'linestyle': None, 'marker': 'o', 'markersize': 8},
        {'data':(163.1, 87.1),  'color': colors[4], 'label': r'WT2 $\tau=0.33$ RFOIL', 'linestyle': None, 'marker': 'o', 'markersize': 8},
        {'data':(157.1, 75.4),  'color': colors[5], 'label': r'WT2 $\tau=0.36$ RFOIL', 'linestyle': None, 'marker': 'o', 'markersize': 8},
        #
        {'data':(255.8, 122.2), 'color': colors[0], 'label': r'WT2 $\tau=0.21$ XFOIL', 'linestyle': None, 'marker': 'X', 'markersize': 8},
        {'data':(175.2, 118.0), 'color': colors[1], 'label': r'WT2 $\tau=0.24$ XFOIL', 'linestyle': None, 'marker': 'X', 'markersize': 8},
        {'data':(182.0, 109.0), 'color': colors[2], 'label': r'WT2 $\tau=0.27$ XFOIL', 'linestyle': None, 'marker': 'X', 'markersize': 8},
        {'data':(177.8, 98.0),  'color': colors[3], 'label': r'WT2 $\tau=0.30$ XFOIL', 'linestyle': None, 'marker': 'X', 'markersize': 8},
        {'data':(166.6, 90.0),  'color': colors[4], 'label': r'WT2 $\tau=0.33$ XFOIL', 'linestyle': None, 'marker': 'X', 'markersize': 8},
        {'data':(159.7, 78.1),  'color': colors[5], 'label': r'WT2 $\tau=0.36$ XFOIL', 'linestyle': None, 'marker': 'X', 'markersize': 8},
    ],
    'legend_entries': [
        [
        {'text': r'$\tau=0.21$', 'linecolor': colors[0], 'linestyle': '-', 'linewidth':2.0, 'marker': None, 'markersize': None},
        {'text': r'$\tau=0.24$', 'linecolor': colors[1], 'linestyle': '-', 'linewidth':2.0, 'marker': None, 'markersize': None},
        {'text': r'$\tau=0.27$', 'linecolor': colors[2], 'linestyle': '-', 'linewidth':2.0, 'marker': None, 'markersize': None},
        {'text': r'$\tau=0.30$', 'linecolor': colors[3], 'linestyle': '-', 'linewidth':2.0, 'marker': None, 'markersize': None},
        {'text': r'$\tau=0.33$', 'linecolor': colors[4], 'linestyle': '-', 'linewidth':2.0, 'marker': None, 'markersize': None},
        {'text': r'$\tau=0.36$', 'linecolor': colors[5], 'linestyle': '-', 'linewidth':2.0, 'marker': None, 'markersize': None},
        ],
        [
        {'text': r'Unconstrained',      'linecolor': 'k', 'linestyle': '-',  'linewidth':2.0, 'marker': None, 'markersize': None},
        {'text': r'Moment Constrained', 'linecolor': 'k', 'linestyle': '-.', 'linewidth':2.0, 'marker': None, 'markersize': None},
        {'text': r'OSO-2025-WT2 RFOIL', 'linecolor': 'k', 'linestyle': None, 'linewidth':2.0, 'marker': 'o', 'markersize': 8},
        {'text': r'OSO-2025-WT2 XFOIL', 'linecolor': 'k', 'linestyle': None, 'linewidth':2.0, 'marker': 'X', 'markersize': 8},
        ]
    ],
    'plot_metadata': {'title': 'Pareto Plots for an IEA 22MW Airfoil Family', 'xlim': [70, 135], 'ylim': [90,300], 'xlabel': 'Rough L/D', 'ylabel': 'Clean L/D', 'grid': True, 
                        'savefilename': 'pareto_plot.png', 'dpi': 300, 'legend_fontsize': 10, 'legend_loc': 'upper left', 'fontsize': 14, 'figsize': (12,9)},
    'te_gap_lookup' : {
        '15':  0.00196,
        '18':  0.00230,
        '21':  0.00262,
        '24':  0.00751,
        '27':  0.01012,
        '30':  0.01140,
        '33':  0.01140,
        '36':  0.01140,
    },
    'data_index_mapping': {'standard': {'obj1ix': 21, 'obj2ix': 20},
                           'exceptions': {'c109_t21_l15_k16_g2000_n752__2026_01_07_12-47': {'obj1ix': 22, 'obj2ix': 21}}},
}

# write to json
with open('pareto_plot_config.json', 'w') as f:
    json.dump(plot_dict, f, indent=4)



#     'neuralfoil_folders': [
#         "/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_111/c111_t24_l14_k16_g2000_n752__2026_03_07_14-35",
#     ],
#     'neuralfoil_postcompute': [
#         "/Users/codykarcher/Dropbox/research/workstation/wt_airfoil_case_111/population_c111_t24_l14_k16_g2000_n752_g882__2026_03_16_15-47"
#     ],
# }


# ── read JSON and plot ──────────────────────────────────────────────────────

def _gen_num(fname):
    """Sort key: parse the generation number from the trailing g### token."""
    last = fname[:-4].split('_')[-1]
    return int(last[1:]) if last.startswith('g') else -1

def _get_indices(dim_map, loc):
    """Return (obj1ix, obj2ix) for a given directory path.
    Supports both old flat {'obj1ix':...} and new {'standard':..., 'exceptions':{...}} formats."""
    if 'standard' in dim_map:
        m = dim_map.get('exceptions', {}).get(os.path.basename(loc), dim_map['standard'])
    else:
        m = dim_map
    return m['obj1ix'], m['obj2ix']

def _load_pts(loc, dim_map):
    """Load rank-1 Pareto points from one directory.
    Returns Nx2 array of (obj1, obj2), or None if directory/files not found."""
    if not os.path.isdir(loc):
        print(f"  skipping (not found): {os.path.basename(loc)}")
        return None
    txts = sorted([f for f in os.listdir(loc) if f.endswith('.txt')], key=_gen_num)
    if not txts:
        print(f"  skipping (no .txt files): {os.path.basename(loc)}")
        return None
    raw = np.loadtxt(os.path.join(loc, txts[-1]))
    p1  = np.abs(raw[raw[:, -1] == 1])
    o1, o2 = _get_indices(dim_map, loc)
    return p1[:, [o1, o2]]   # Nx2: col0=obj1 (x), col1=obj2 (y)

def _combined_pareto(pts):
    """Return the non-dominated subset of an Nx2 array (maximise both objectives),
    sorted by col0 ascending for plotting."""
    pts_sorted = pts[pts[:, 0].argsort()[::-1]]   # sort by obj1 descending
    front, best_obj2 = [], -np.inf
    for p in pts_sorted:
        if p[1] > best_obj2:
            front.append(p)
            best_obj2 = p[1]
    front = np.array(front)
    return front[front[:, 0].argsort()]            # re-sort ascending for plotting

with open('pareto_plot_config.json', 'r') as fp:
    cfg = json.load(fp)

dim_map = cfg['data_index_mapping']
meta    = cfg.get('plot_metadata', cfg.get('metadata', {}))

fontsize        = meta.get('fontsize', 12)
legend_fontsize = meta.get('legend_fontsize', fontsize)
legend_loc      = meta.get('legend_loc', 'best')
dpi             = meta.get('dpi', 150)

fig, ax = plt.subplots(figsize=meta.get('figsize', (9, 6)))

for entry in cfg['data']:
    color      = entry['color']
    linestyle  = entry.get('linestyle')
    linewidth  = entry.get('linewidth', 1.5)
    marker     = entry.get('marker')
    markersize = entry.get('markersize') or 6

    if 'data_location' in entry:
        loc = entry['data_location']

        if isinstance(loc, list):
            # gather rank-1 points from every directory, then compute combined Pareto front
            all_pts = [_load_pts(l, dim_map) for l in loc]
            all_pts = [p for p in all_pts if p is not None]
            if not all_pts:
                continue
            p1 = _combined_pareto(np.vstack(all_pts))
        else:
            pts = _load_pts(loc, dim_map)
            if pts is None:
                continue
            p1 = pts[pts[:, 0].argsort()]

        ax.plot(p1[:, 0], p1[:, 1],
                color=color, linestyle=linestyle, linewidth=linewidth)

    elif 'data' in entry:
        d = entry['data']   # stored as (clean L/D, rough L/D)
        ax.plot(d[1], d[0],
                marker=marker, color=color,
                markersize=markersize, linestyle='none',
                markeredgewidth=0.5, markeredgecolor='k')

# build one handle-list per group (one group → one legend column)
# matplotlib fills legends column-by-column, so concatenate groups after padding
group_handles = []
for group in cfg.get('legend_entries', []):
    gh = []
    for le in group:
        lc = le.get('linecolor') or '#000000'
        ls = le.get('linestyle') or 'none'
        lw = le.get('linewidth', 1.5)
        mk = le.get('marker')    or 'none'
        ms = le.get('markersize') or 8
        h  = mlines.Line2D([], [], color=lc, linestyle=ls, linewidth=lw,
                           marker=mk, markersize=ms,
                           markeredgewidth=0.5 if mk != 'none' else 0,
                           markeredgecolor='k', label=le['text'])
        gh.append(h)
    group_handles.append(gh)

# pad shorter columns with invisible handles so all columns have equal row count
ncols   = len(group_handles)
max_len = max(len(g) for g in group_handles)
blank   = lambda: mlines.Line2D([], [], linestyle='none', label='')
for g in group_handles:
    while len(g) < max_len:
        g.append(blank())

# concatenate groups: matplotlib column-major fill puts group[0] in col0, group[1] in col1, etc.
flat_handles = []
for g in group_handles:
    flat_handles.extend(g)

ax.legend(handles=flat_handles, ncols=ncols, fontsize=legend_fontsize,
          framealpha=0.9, columnspacing=1.5, loc=legend_loc)
ax.set_xlabel(meta.get('xlabel', 'Rough L/D'), fontsize=fontsize)
ax.set_ylabel(meta.get('ylabel', 'Clean L/D'), fontsize=fontsize)
ax.set_title(meta.get('title', ''), fontsize=fontsize)
if meta.get('xlim'):        ax.set_xlim(meta['xlim'])
if meta.get('ylim'):        ax.set_ylim(meta['ylim'])
if meta.get('grid', True):  ax.grid(True, alpha=0.3)

plt.tight_layout()
save_name = meta.get('savefilename', 'pareto_plot.png')
plt.savefig(save_name, dpi=dpi, bbox_inches='tight')
# plt.show()
print(f"Saved → {save_name}  (dpi={dpi})")


# ── export Pareto-front rows to text files ──────────────────────────────────
import json, os, re
import numpy as np

def _gen_num(fname):
    last = fname[:-4].split('_')[-1]
    return int(last[1:]) if last.startswith('g') else -1

def _get_indices(dim_map, loc):
    if 'standard' in dim_map:
        m = dim_map.get('exceptions', {}).get(os.path.basename(loc), dim_map['standard'])
    else:
        m = dim_map
    return m['obj1ix'], m['obj2ix']

def _load_rank1_full(loc, dim_map):
    """Return (full_rank1_rows, o1ix, o2ix) for the latest .txt in loc, or None.
    Rows are kept with their original signs — abs is NOT applied here."""
    if not os.path.isdir(loc):
        print(f"  skipping (not found): {os.path.basename(loc)}")
        return None
    txts = sorted([f for f in os.listdir(loc) if f.endswith('.txt')], key=_gen_num)
    if not txts:
        print(f"  skipping (no .txt files): {os.path.basename(loc)}")
        return None
    raw  = np.loadtxt(os.path.join(loc, txts[-1]))
    rows = raw[raw[:, -1] == 1]          # rank-1 rows, signs preserved
    o1, o2 = _get_indices(dim_map, loc)
    return rows, o1, o2

def _pareto_rows(data_list):
    """
    data_list : list of (rows_Nx_full, o1ix, o2ix)
    Returns list of (rough_LD, clean_LD, full_row) on the combined Pareto front,
    sorted by rough_LD ascending.  Maximises both objectives.
    abs() applied only to the objective columns (stored negative in some runs).
    """
    pts = []
    for rows, o1, o2 in data_list:
        for row in rows:
            pts.append((abs(row[o1]), abs(row[o2]), row))   # (rough_LD, clean_LD, full_row)
    pts.sort(key=lambda x: -x[0])                           # sort by rough_LD descending
    front, best_clean = [], -np.inf
    for rough, clean, row in pts:
        if clean > best_clean:
            front.append((rough, clean, row))
            best_clean = clean
    front.sort(key=lambda x: x[0])                          # re-sort ascending for output
    return front

def _label_to_filename(label):
    tau_m = re.search(r'tau=0\.(\d+)', label)
    tau   = tau_m.group(1) if tau_m else 'unknown'
    kind  = 'moment_constrained' if 'Moment' in label else 'unconstrained'
    return f"pareto_t{tau}_{kind}.txt"

with open('pareto_plot_config.json', 'r') as fp:
    cfg = json.load(fp)

dim_map   = cfg['data_index_mapping']
out_dir   = str(path_to_oso / 'released_designs/pareto_data')
os.makedirs(out_dir, exist_ok=True)

for entry in cfg['data']:
    if 'data_location' not in entry:
        continue

    label = entry.get('label', '')
    loc   = entry['data_location']
    locs  = loc if isinstance(loc, list) else [loc]

    data_list = [r for r in (_load_rank1_full(l, dim_map) for l in locs) if r is not None]
    if not data_list:
        continue

    front = _pareto_rows(data_list)

    # columns: first 16 design vars (original signs) | clean L/D | rough L/D
    rows_out = np.array([[*row[:16], clean, rough] for rough, clean, row in front])

    fname = os.path.join(out_dir, _label_to_filename(label))
    np.savetxt(fname, rows_out, fmt='%.6f')
    print(f"Wrote {fname}  ({len(rows_out)} Pareto points)")
