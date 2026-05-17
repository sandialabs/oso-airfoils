import matplotlib.pyplot as plt
import numpy as np
import os
import natsort
import json
from PIL import Image
from kulfan import Kulfan

from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

from compare_airfoils import rainbow_plot

import sys
# Usage:
#   python -m mpi4py generate_gif.py <path_to_data> [output_gif_name] [frames_dir_name] [dpi] [frame_duration_ms]
#   All arguments after path_to_data are optional.
# Example:
#   mpirun -n 8 python -m mpi4py generate_gif.py cases/cases_111_to_120/case_114/c114_t21_k16_n752_l13_e15__2026_05_14_02-12-5450 pareto_shapes_evolution.gif _combined_gif_frames 100 100

path_to_data = None
output_gif_name = 'pareto_shapes_evolution.gif'
frames_dir_name = '_combined_gif_frames'
image_dpi = 100
frame_duration = 100  # ms per frame


if len(sys.argv) > 1:
    path_to_data = sys.argv[1]
if len(sys.argv) > 2:
    output_gif_name = sys.argv[2]
if len(sys.argv) > 3:
    frames_dir_name = sys.argv[3]
if len(sys.argv) > 4:
    try:
        image_dpi = int(sys.argv[4])
    except Exception:
        print(f"Warning: Could not parse dpi '{sys.argv[4]}', using default {image_dpi}")
if len(sys.argv) > 5:
    try:
        frame_duration = int(sys.argv[5])
    except Exception:
        print(f"Warning: Could not parse frame duration '{sys.argv[5]}', using default {frame_duration}")

assert path_to_data is not None and os.path.isdir(path_to_data), "First argument must be a valid directory path."

files = natsort.natsorted([f for f in os.listdir(path_to_data) if '.json' in f and 'population' in f], alg=natsort.ns.IGNORECASE)
# print(files[-1])

nafl = 20

# xlim_bottom is always fixed; all other limits are derived from the last population file.
xlim_bottom = (-0.05, 1.05)

# ── Auto-limit tuning factors (multiplicative: 1.0 = tight fit, >1.0 adds padding) ──
xlim_top_factor    = 1.40   # x-axis of K-parameter plot  (rough L/D range)
ylim_top_factor    = 1.20   # y-axis of K-parameter plot  (Kulfan coeff range)
xlim_pareto_factor = 2.00   # x-axis of Pareto plot       (rough L/D range)
ylim_pareto_factor = 1.60   # y-axis of Pareto plot       (clean L/D range)
ylim_bottom_factor = 1.40   # y-axis of airfoil shape plot

# ── Derive limits from the last population file ───────────────────────────────
def _auto_lim(vals, factor):
    lo, hi = float(np.min(vals)), float(np.max(vals))
    ctr = (lo + hi) / 2.0
    half = max((hi - lo) / 2.0, 1e-6) * factor
    return [ctr - half, ctr + half]

_last_data  = json.load(open(os.path.join(path_to_data, files[-1])))
_pareto     = [p for p in _last_data['population'] if p['pareto_index'] == 1]
_rough      = np.array([p['LoD_rough_at_design'] for p in _pareto])
_clean      = np.array([p['LoD_clean_at_design'] for p in _pareto])
_K_all      = np.concatenate([list(p['K_upper']) + list(p['K_lower']) for p in _pareto])

_afl_ys = []
for _p in _pareto:
    _afl = Kulfan(TE_gap=_last_data['input_parameters']['TE_gap'])
    _afl.upperCoefficients = _p['K_upper']
    _afl.lowerCoefficients = _p['K_lower']
    _afl_ys.extend(_afl.ycoordinates)
_afl_ys = np.array(_afl_ys)

xlim_top    = _auto_lim(_rough,  xlim_top_factor)
ylim_top    = _auto_lim(_K_all,  ylim_top_factor)
xlim_pareto = _auto_lim(_rough,  xlim_pareto_factor)
ylim_pareto = _auto_lim(_clean,  ylim_pareto_factor)
ylim_bottom = tuple(_auto_lim(_afl_ys, ylim_bottom_factor))

legend_loc = 'upper right'
fix_colormap_range = None  # (min_val, max_val) to fix; None to auto

# Font size hook: applies to titles, axis labels, and tick labels (NOT the legend)
font_size = 20
legend_font_size = 12  # set independently to leave the legend unaffected

turbo_cmap = plt.get_cmap('turbo', nafl)


frames_dir = os.path.join(path_to_data, frames_dir_name)
os.makedirs(frames_dir, exist_ok=True)
frame_paths = []

gen_indices = range(0, len(files))

for ix, gen_i in enumerate(gen_indices):
    if ix % size != rank:
        continue  # Skip this generation; it's not assigned to this process
    f = files[gen_i]
    data = json.load(open(os.path.join(path_to_data, f)))
    pop = data['population']
    pareto_points = [p for p in pop if p['pareto_index'] == 1]
    pareto_points = sorted(pareto_points, key=lambda p: p['LoD_rough_at_design'])

    seismic_cmap = plt.get_cmap('seismic', int(data['input_parameters']['N_k']))

    rough_LD_vals = np.array([p['LoD_rough_at_design'] for p in pareto_points])
    clean_LD_vals = np.array([p['LoD_clean_at_design'] for p in pareto_points])

    # Sample approximately evenly along arc length of the Pareto front.
    # Normalize each axis by its plotted (or data) range so neither axis
    # dominates the arc-length computation.
    if xlim_pareto is not None:
        x_norm = float(xlim_pareto[1] - xlim_pareto[0])
    else:
        x_norm = float(np.ptp(rough_LD_vals)) or 1.0
    if ylim_pareto is not None:
        y_norm = float(ylim_pareto[1] - ylim_pareto[0])
    else:
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
        ixs = np.array([int(np.argmin(np.abs(cum_arc - t))) for t in target_arcs])

    # Build the airfoil shapes for this generation
    shape_curves = []
    for i in list(reversed(ixs)):
        afl = Kulfan(TE_gap=data['input_parameters']['TE_gap'])
        afl.upperCoefficients = pareto_points[i]['K_upper']
        afl.lowerCoefficients = pareto_points[i]['K_lower']

        rough_ld_value = rough_LD_vals[i]
        if fix_colormap_range is not None:
            cmap_min, cmap_max = fix_colormap_range
            clamped_value = np.clip(rough_ld_value, cmap_min, cmap_max)
            color_idx = (clamped_value - cmap_min) / (cmap_max - cmap_min) * (nafl - 1)
            color_idx = int(np.clip(color_idx, 0, nafl - 1))
        else:
            color_idx = ixs.tolist().index(i)

        shape_curves.append((np.asarray(afl.xcoordinates),
                             np.asarray(afl.ycoordinates),
                             turbo_cmap(color_idx)))

    x_range_bottom = xlim_bottom[1] - xlim_bottom[0]

    if ylim_bottom is None:
        all_y = np.concatenate([yc for _, yc, _ in shape_curves])
        ymin, ymax = float(np.min(all_y)), float(np.max(all_y))
        pad = 0.05 * (ymax - ymin if ymax > ymin else 1.0)
        ylim_b = (ymin - pad, ymax + pad)
    else:
        ylim_b = ylim_bottom

    y_range_bottom = ylim_b[1] - ylim_b[0]

    top_w_in, top_h_in = 20.0, 11.0
    bottom_h_in = top_w_in * (y_range_bottom / x_range_bottom)
    right_col_h_in = top_h_in + bottom_h_in
    fig_h_in = right_col_h_in

    # Left column has same width as right column; total fig width = 2 * top_w_in
    fig = plt.figure(figsize=(2 * top_w_in, fig_h_in))
    gs = fig.add_gridspec(2, 2,
                          width_ratios=[top_w_in, top_w_in],
                          height_ratios=[top_h_in, bottom_h_in],
                          hspace=0.15, wspace=0.12)
    ax_left = fig.add_subplot(gs[:, 0])
    ax_top = fig.add_subplot(gs[0, 1])
    ax_bot = fig.add_subplot(gs[1, 1])

    # --- Left subplot: pareto front (clean vs rough L/D) ---
    # pop = data['population']
    pareto_indices = np.sort(np.unique([p['pareto_index'] for p in pop]))
    seismic_cmap = plt.get_cmap('seismic', int(data['input_parameters']['N_k']))
    for pix in list(reversed(pareto_indices)):
        pareto_points = [p for p in pop if p['pareto_index'] == pix]
        pareto_points = sorted(pareto_points, key=lambda p: p['LoD_rough_at_design'])
        rough_LD_vals = np.array([p['LoD_rough_at_design'] for p in pareto_points])
        clean_LD_vals = np.array([p['LoD_clean_at_design'] for p in pareto_points])
        # if pix == 1:
        #     ax_left.plot(rough_LD_vals, clean_LD_vals, 'ko', ms=5, mew=0)
        # else:
        #     # ax_left.plot(rough_LD_vals, clean_LD_vals, 'o', color='%f'%(max([1-1/pix, 0.3])), ms=5, mew=0)
        clrval = min([1-1/(pix+1), 0.91])
        ax_left.plot(rough_LD_vals, clean_LD_vals, 'o', color=[clrval,clrval,clrval], ms=4, mew=0)
            
    pareto_points = [p for p in pop if p['pareto_index'] == 1]
    pareto_points = sorted(pareto_points, key=lambda p: p['LoD_rough_at_design'])
    rough_LD_vals = np.array([p['LoD_rough_at_design'] for p in pareto_points])
    clean_LD_vals = np.array([p['LoD_clean_at_design'] for p in pareto_points])

    # Overlay the selected airfoils with their turbo colors
    for i in ixs:
        rough_ld_value = rough_LD_vals[i]
        clean_ld_value = clean_LD_vals[i]
        if fix_colormap_range is not None:
            cmap_min, cmap_max = fix_colormap_range
            clamped_value = np.clip(rough_ld_value, cmap_min, cmap_max)
            color_idx = (clamped_value - cmap_min) / (cmap_max - cmap_min) * (nafl - 1)
            color_idx = int(np.clip(color_idx, 0, nafl - 1))
        else:
            color_idx = ixs.tolist().index(i)
        ax_left.plot(rough_ld_value, clean_ld_value, 'o',
                     color=turbo_cmap(color_idx), ms=9)

    ax_left.set_xlabel('Rough L/D', fontsize=font_size)
    ax_left.set_ylabel('Clean L/D', fontsize=font_size)
    ax_left.set_title('Pareto Front, Generation %d' % gen_i, fontsize=font_size)
    ax_left.tick_params(axis='both', labelsize=font_size)
    ax_left.grid(True)
    if xlim_pareto is not None:
        ax_left.set_xlim(xlim_pareto)
    if ylim_pareto is not None:
        ax_left.set_ylim(ylim_pareto)

    # Top: parameters vs rough L/D
    cmap_counter = 0
    for i in range(0, int(data['input_parameters']['N_k']/2)):
        yvals = np.array([p['K_upper'][i] for p in pareto_points])
        ax_top.plot(rough_LD_vals, yvals, 'o-', color=seismic_cmap(cmap_counter), ms=3,
                    label='Upper Surface %d' % (i+1))
        cmap_counter += 1

    cmap_counter = int(data['input_parameters']['N_k']) - 1
    for i in range(0, int(data['input_parameters']['N_k']/2)):
        yvals = np.array([p['K_lower'][i] for p in pareto_points])
        ax_top.plot(rough_LD_vals, yvals, 'o-', color=seismic_cmap(cmap_counter), ms=3,
                    label='Lower Surface %d' % (i+1))
        cmap_counter -= 1

    for i in ixs:
        rough_ld_value = rough_LD_vals[i]
        if fix_colormap_range is not None:
            cmap_min, cmap_max = fix_colormap_range
            clamped_value = np.clip(rough_ld_value, cmap_min, cmap_max)
            color_idx = (clamped_value - cmap_min) / (cmap_max - cmap_min) * (nafl - 1)
            color_idx = int(np.clip(color_idx, 0, nafl - 1))
        else:
            color_idx = ixs.tolist().index(i)
        ax_top.plot(rough_ld_value, 0, 'o', color=turbo_cmap(color_idx), ms=5)

    ax_top.legend(loc=legend_loc, fontsize=legend_font_size)
    ax_top.set_xlabel('Rough L/D', fontsize=font_size)
    ax_top.set_ylabel('Airfoil Shape Parameters', fontsize=font_size)
    ax_top.set_title('Pareto Front Airfoil Shapes, Generation %d' % gen_i, fontsize=font_size)
    ax_top.tick_params(axis='both', labelsize=font_size)
    ax_top.grid(True)
    if xlim_top is not None:
        ax_top.set_xlim(xlim_top)
    if ylim_top is not None:
        ax_top.set_ylim(ylim_top)

    # Bottom: airfoil shapes
    for x, y, c in shape_curves:
        ax_bot.plot(x, y, color=c, alpha=1)

    ax_bot.set_xlim(xlim_bottom)
    ax_bot.set_ylim(ylim_b)
    ax_bot.set_aspect('auto')
    ax_bot.grid(True)
    ax_bot.set_xlabel('x/c', fontsize=font_size)
    ax_bot.set_ylabel('y/c', fontsize=font_size)
    ax_bot.tick_params(axis='both', labelsize=font_size)

    frame_path = os.path.join(frames_dir, f'frame_{ix:04d}.png')
    fig.savefig(frame_path, dpi=image_dpi, bbox_inches='tight')
    frame_paths.append(frame_path)
    plt.close(fig)

frame_paths = comm.gather(frame_paths, root=0)

if rank == 0:
    frame_paths = [fp for sublist in frame_paths for fp in sublist]  # Flatten the list of lists
    frame_paths = natsort.natsorted(frame_paths, alg=natsort.ns.IGNORECASE)

    # Build the GIF from saved frames, closing files after use to avoid 'too many open files'
    gif_path = os.path.join(path_to_data, output_gif_name)
    def frame_generator(paths):
        for p in paths:
            with Image.open(p) as im:
                yield im.copy()

    frames_iter = frame_generator(frame_paths)
    try:
        first_frame = next(frames_iter)
    except StopIteration:
        first_frame = None

    if first_frame is not None:
        frames_list = list(frames_iter)
        first_frame.save(
            gif_path,
            save_all=True,
            append_images=frames_list,
            duration=frame_duration,   # ms per frame
            loop=0,
        )
        print(f'GIF saved to {gif_path}')


if rank == 0:
    path_to_data = path_to_data + os.sep + files[-1]
    rainbow_plot(path_to_data)