"""process_paretos.py
--------------------
Post-processes the pre-extracted Pareto-front txt files in
released_designs/pareto_data/{unconstrained,moment_constrained}/tau_XX/.

Three tasks per tau / constraint-type combination:
  1. Rainbow polar plot  — N_RAINBOW airfoils across the Pareto front + optional
                           historical comparison airfoils
  2. Extreme-end comparison — max-clean-L/D and max-rough-L/D endpoints vs.
                              historical airfoils
  3. Chebyshev + polynomial shape-function fit → shape_functions.json

Run from the pareto_data directory:
    python process_paretos.py
"""

import sys, os, json, shutil
import pathlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

path_to_here = pathlib.Path(__file__).parent.resolve()
path_to_oso  = path_to_here.parent.parent   # released_designs/pareto_data → oso-airfoils/

# Add postprocessing to sys.path so compare_airfoils, neuralfoil/xfoil wrappers,
# and the units-capable kulfan.py are all importable without modifying those files.
# sys.path.insert(0, str(path_to_oso / 'postprocessing'))
import compare_airfoils as _ca_mod
from compare_airfoils import compare_airfoils
from kulfan import Kulfan, units
from neuralfoil_wrapper_noprint import run as _run_neuralfoil
from xfoil_wrapper_noprint import run as _run_xfoil


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

TAUS             = ['21', '24', '27', '30', '33', '36']
CONSTRAINT_TYPES = ['unconstrained', 'moment_constrained']

N_RAINBOW        = 21    # airfoils sampled for the rainbow polar plot
N_PARETO_SAMPLES = 75    # airfoils sampled along the Pareto front for shape fitting
N_DIAG_PLOT      = 21    # airfoils shown in the airfoil-family diagnostic plot
CBY_ORDER        = 20    # Chebyshev basis order for thickness/camber profile fits
POLY_DEG         = 4     # polynomial degree mapping Pareto parameter x → coefficient

# Pareto parameter value used as the "rough" extreme endpoint (0 = max clean, 1 = max rough)
# 1.0 selects the actual maximum-rough-L/D airfoil; 0.9 selects the 90th-percentile position
ROUGH_EXTREME    = 0.85

# Set False to hide the Cp_min panel in all polar plots (cleaner for most comparisons)
SHOW_CPMIN       = False

# Set False to skip all 6 diagnostic plots (speeds up a full re-run significantly)
PRODUCE_DIAG_PLOTS = False

# Apply Cp_min toggle — patches polarPlot in the compare_airfoils module at import time
_orig_polarPlot = _ca_mod.polarPlot
if not SHOW_CPMIN:
    def _no_cpmin(dataList, *args, **kwargs):
        kwargs.setdefault('show_cpmin', False)
        return _orig_polarPlot(dataList, *args, **kwargs)
    _ca_mod.polarPlot = _no_cpmin

# Reynolds number per tau (adjust to match your run conditions)
re_lookup = {
    '15': 12e6,
    '18': 12e6,
    '21': 12e6,
    '24': 13e6,
    '27': 16e6,
    '30': 18e6,
    '33': 16e6,
    '36': 13e6,
}

cl_lookup = {
    '15': 1.5,
    '18': 1.5,
    '21': 1.5,
    '24': 1.4,
    '27': 1.3,
    '30': 1.2,
    '33': 1.2,
    '36': 1.2,
}

# TE gap per tau (from optimisation defaults)
te_gap_lookup = {
    '15': 0.00196,
    '18': 0.00230,
    '21': 0.00262,
    '24': 0.00751,
    '27': 0.01012,
    '30': 0.01140,
    '33': 0.01140,
    '36': 0.01140,
}

# Turbulence cases: [N_crit, xtp_upper, xtp_lower]
#   case 0 → clean  (free transition, N_crit=9)
#   case 1 → rough  (forced at 5%,   N_crit=3)
turb_cases = [[9, 1.0, 1.0], [3, 0.05, 0.05]]


# ══════════════════════════════════════════════════════════════════════════════
# COMPARISON AIRFOILS
# ══════════════════════════════════════════════════════════════════════════════
# Dict keyed by tau string.  Each value is a dict of {display_name: dat_file_path}.
# Comment / uncomment individual entries to control which airfoils appear in
# the rainbow plot and the extreme-end comparison.

_h = path_to_oso / 'historical_airfoils'
_d = path_to_oso / 'released_designs/OSO_2025_WT2/datfiles'


comparison_airfoils = {
    '21': {
        'OSO-2025-WT2-21': str(_d / 'OSO_2025_WT2_T21.dat'),
        # 'FFA-W1-211':    str(_h / 'ffa/fitted/FFA-W1-211_fittedCST10.dat'),
        # 'FFA-W2-210':    str(_h / 'ffa/fitted/FFA-W2-210_fittedCST10.dat'),
        'FFA-W3-211':  str(_h / 'ffa/fitted/FFA-W3-211_fittedCST10.dat'),
        # 'DU-93-W-210':   str(_h / 'du/du_93-w-210.dat'),
        # 'RISO-A-21':   str(_h / 'riso-a/riso-a-21.dat'),
    },
    '24': {
        'OSO-2025-WT2-24': str(_d / 'OSO_2025_WT2_T24.dat'),
        # 'FFA-W1-242':    str(_h / 'ffa/fitted/FFA-W1-242_fittedCST10.dat'),
        'FFA-W3-241':  str(_h / 'ffa/fitted/FFA-W3-241_fittedCST10.dat'),
        # 'MHKF1-240':     str(_h / 'mhkf1/mhkf1-240.dat'),
        # 'RISO-A-24':   str(_h / 'riso-a/riso-a-24.dat'),
    },
    '27': {
        'OSO-2025-WT2-27': str(_d / 'OSO_2025_WT2_T27.dat'),
        # 'FFA-W1-271':    str(_h / 'ffa/fitted/FFA-W1-271_fittedCST10.dat'),
        'FFA-W3-270':  str(_h / 'ffa/fitted/FFA-W3-270_fittedCST10.dat'),
        # 'DU-91-W2-250':  str(_h / 'du/du_91-w2-250.dat'),
        # 'RISO-A-27':   str(_h / 'riso-a/riso-a-27.dat'),
    },
    '30': {
        'OSO-2025-WT2-30': str(_d / 'OSO_2025_WT2_T30.dat'),
        'FFA-W3-301':    str(_h / 'ffa/fitted/FFA-W3-301_fittedCST10.dat'),
        # 'DU-97-W-300':   str(_h / 'du/du_97-w-300.dat'),
        # 'RISO-A-30':   str(_h / 'riso-a/riso-a-30.dat'),
    },
    '33': {
        'OSO-2025-WT2-33': str(_d / 'OSO_2025_WT2_T33.dat'),
        'FFA-W3-332':    str(_h / 'ffa/fitted/FFA-W3-332_fittedCST10.dat'),
        # 'RISO-B-29':   str(_h / 'riso-b/riso-b-29.dat'),
        # 'RISO-B-35':     str(_h / 'riso-b/riso-b-35.dat'),
    },
    '36': {
        'OSO-2025-WT2-36': str(_d / 'OSO_2025_WT2_T36.dat'),
        'FFA-W3-360':    str(_h / 'ffa/fitted/FFA-W3-360_fittedCST10.dat'),
        # 'RISO-B-35':     str(_h / 'riso-b/riso-b-35.dat'),
    },
}

# Colors for comparison airfoils in the rainbow plot (cycling if more than 4)
# _comp_colors = ['#0065cc', '#eea800', '#009e73', '#d55e00', '#7860aa', '#ede13f', '#56b4ff', '#fca7c7', '#5d5d5d', '#000000']
# Only plotting against the oso airfoils
# _comp_colors = ['#fca7c7']
#['#000000', '#7b3f00', '#3d0066', '#006633']


# ══════════════════════════════════════════════════════════════════════════════
# HELPERS — DATA LOADING
# ══════════════════════════════════════════════════════════════════════════════
# File layout written by plot_paretos.py:
#   cols 0–7   : K_upper  (N_k=16, so 8 upper coefficients)
#   cols 8–15  : K_lower  (8 lower coefficients)
#   col  16    : clean L/D
#   col  17    : rough L/D
# Rows are sorted ascending by rough L/D.

def _load_pareto_txt(fpath):
    """Return (K_upper, K_lower, clean_LD, rough_LD) arrays — one row per point."""
    data     = np.loadtxt(fpath)
    K_upper  = data[:, 0:8]
    K_lower  = data[:, 8:16]
    clean_LD = data[:, 16]
    rough_LD = data[:, 17]
    return K_upper, K_lower, clean_LD, rough_LD


def _arc_sample_indices(clean_LD, rough_LD, n):
    """n indices evenly spaced by normalised arc length along the Pareto front."""
    x_norm  = float(np.ptp(rough_LD)) or 1.0
    y_norm  = float(np.ptp(clean_LD)) or 1.0
    cum_arc = np.concatenate([[0.0], np.cumsum(
        np.sqrt(np.diff(rough_LD / x_norm)**2 + np.diff(clean_LD / y_norm)**2))])
    targets = np.linspace(0.0, cum_arc[-1], n)
    return [int(np.argmin(np.abs(cum_arc - t))) for t in targets]


def _x_pareto(clean_LD, rough_LD):
    """Normalised arc-length parameter: 0 = max-clean end, 1 = max-rough end."""
    pts = np.stack([rough_LD, clean_LD], axis=1)
    arc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
    arc /= arc[-1] if arc[-1] > 0 else 1.0
    i_clean = int(np.argmax(clean_LD))
    i_rough = int(np.argmax(rough_LD))
    if arc[i_clean] < arc[i_rough]:
        return (arc - arc[i_clean]) / (arc[i_rough] - arc[i_clean])
    else:
        return 1.0 - (arc - arc[i_rough]) / (arc[i_clean] - arc[i_rough])


def _kulfan_entry(K_upper_row, K_lower_row, te_gap):
    """Build the afl_dict entry accepted by compare_airfoils for a Kulfan airfoil."""
    return {
        'K_upper': K_upper_row.tolist(),
        'K_lower': K_lower_row.tolist(),
        'TE_gap':  te_gap,
    }


def _max_LD_from_shape(x_par_val, thickness_poly, camber_poly, te_gap, psi_vals, re, cl_design):
    """Reconstruct an airfoil from shape-function polynomials at x_par_val, fit 8
    Kulfan coefficients via fit2coordinates, then evaluate clean and rough L/D at
    alpha_design (where cl_clean == cl_design) using xfoil, matching the optimiser."""
    t_params  = [np.polyval(p, x_par_val) for p in thickness_poly]
    c_params  = [np.polyval(p, x_par_val) for p in camber_poly]
    thickness = shape_fcn(psi_vals, t_params)
    camber    = shape_fcn(psi_vals, c_params)
    tfit      = thickness + te_gap * psi_vals
    xcoords   = np.concatenate([psi_vals[::-1], psi_vals[1:]])
    ycoords   = np.concatenate([(camber + tfit / 2)[::-1], (camber - tfit / 2)[1:]])

    afl_fit = Kulfan(TE_gap=te_gap)
    afl_fit.fit2coordinates(xcoords, ycoords, fit_order=8)
    K_u = afl_fit.upperCoefficients.magnitude.tolist()
    K_l = afl_fit.lowerCoefficients.magnitude.tolist()

    # Clean polar — matching optimiser settings (xfoil, N_crit=9)
    res1        = _run_xfoil('alfa', K_u, K_l, [0, 30, 1],
                              Re=re, N_crit=9.0, xtp_u=1.0, xtp_l=1.0, TE_gap=te_gap)
    cl_clean    = np.array(res1['cl'])
    cd_clean    = np.array(res1['cd'])
    alpha_clean = np.array(res1['alpha'])
    LoD_clean   = np.where(cd_clean > 1e-6, cl_clean / cd_clean, 0.0)

    # Find pre-stall CL peak — exactly matching optimiser logic
    # Range extends beyond array length so IndexError signals "no stall found".
    # At i=mid_c the comparison wraps to cl_clean[-1], setting peak=-1 (last
    # element) when xfoil stops before the array end (typical behaviour).
    mid_c  = int(np.argmin(np.abs(alpha_clean)))
    peak_c = mid_c
    try:
        for i in range(mid_c, mid_c + len(alpha_clean)):
            if cl_clean[i] > cl_clean[i - 1]:
                pass
            else:
                peak_c = i - 1
                break
    except IndexError:
        return np.nan, np.nan
    if cl_clean[peak_c] <= cl_design:
        return np.nan, np.nan  # airfoil cannot reach cl_design

    alpha_design        = float(np.interp(cl_design,
                                          cl_clean[mid_c:peak_c],
                                          alpha_clean[mid_c:peak_c]))
    LoD_clean_at_design = float(np.interp(alpha_design,
                                          alpha_clean[mid_c:peak_c],
                                          LoD_clean[mid_c:peak_c]))

    # Rough polar — matching optimiser settings (xfoil, N_crit=3, forced xtp)
    res2        = _run_xfoil('alfa', K_u, K_l, [0, 20, 1],
                              Re=re, N_crit=3.0, xtp_u=0.05, xtp_l=0.05, TE_gap=te_gap)
    cl_rough    = np.array(res2['cl'])
    cd_rough    = np.array(res2['cd'])
    alpha_rough = np.array(res2['alpha'])
    LoD_rough   = np.where(cd_rough > 1e-6, cl_rough / cd_rough, 0.0)

    mid_r  = int(np.argmin(np.abs(alpha_rough)))
    peak_r = mid_r
    try:
        for i in range(mid_r, mid_r + len(alpha_rough)):
            if cl_rough[i] > cl_rough[i - 1]:
                pass
            else:
                peak_r = i - 1
                break
    except IndexError:
        return np.nan, np.nan

    LoD_rough_at_design = float(np.interp(alpha_design,
                                          alpha_rough[mid_r:peak_r],
                                          LoD_rough[mid_r:peak_r]))
    return LoD_clean_at_design, LoD_rough_at_design


# ══════════════════════════════════════════════════════════════════════════════
# HELPERS — CHEBYSHEV + POLYNOMIAL SHAPE FITTING
# (copied from postprocessing/find_shape_functions.py)
# ══════════════════════════════════════════════════════════════════════════════

def cby(IFUN_in, S):
    SWT   = 2.0
    SW    = SWT * S / (1 + (SWT - 1) * S)
    X     = 1.0 - 2.0 * SW
    THETA = np.arccos(np.clip(X, -1.0, 1.0))
    RF    = float(IFUN_in + 1)
    if IFUN_in % 2 == 0:
        return (X - np.cos(RF * THETA)) / RF
    else:
        return (1.0 - np.cos(RF * THETA)) / RF

def build_basis(s_vals, max_order):
    basis = np.zeros((len(s_vals), max_order + 1))
    for i in range(max_order + 1):
        basis[:, i] = [cby(i, s) for s in s_vals]
    return basis

def evaluate_cby(s_vals, coeffs):
    return build_basis(s_vals, len(coeffs) - 1) @ np.asarray(coeffs)

def shape_fcn(x, params):
    r = params[0]
    return r * np.sqrt(1 - (x - 1)**2) - r * x + evaluate_cby(x, params[1:])

def fit_profile(psi, values):
    p0 = [0.0] + [0.0] * CBY_ORDER
    params, _ = curve_fit(
        lambda x, r, *c: shape_fcn(x, [r] + list(c)),
        psi, values, p0=p0,
    )
    return np.asarray(params)

def _poly_coeffs_dict(poly_list):
    labels = ['r'] + [f'cby_{i}' for i in range(len(poly_list) - 1)]
    return {label: poly_list[i].tolist() for i, label in enumerate(labels)}


# ══════════════════════════════════════════════════════════════════════════════
# MAIN LOOP
# ══════════════════════════════════════════════════════════════════════════════

for tau in TAUS:
    te_gap    = te_gap_lookup[tau]
    re        = re_lookup[tau]
    cl_design = cl_lookup[tau]
    comps     = comparison_airfoils.get(tau, {})
    _combined_endpoints = {}   # keyed by ctype; filled inside the inner loop

    for ctype in CONSTRAINT_TYPES:

        # ── locate file ──────────────────────────────────────────────────────
        txt_path = path_to_here / ctype / f'tau_{tau}' / f'pareto_t{tau}_{ctype}.txt'
        if not txt_path.exists():
            print(f'  [skip] {txt_path}')
            continue

        out_dir = txt_path.parent
        print(f'\n── τ={tau}%  {ctype} ──────────────────────────────────')

        K_upper, K_lower, clean_LD, rough_LD = _load_pareto_txt(txt_path)

        if len(clean_LD) < 2:
            print(f'  [skip] only {len(clean_LD)} Pareto point(s)')
            continue

        # ── 1. Rainbow plot ──────────────────────────────────────────────────
        print('  Building rainbow plot ...')
        ixs = _arc_sample_indices(clean_LD, rough_LD, N_RAINBOW)
        n = len(ixs)
        _turbo = plt.get_cmap('turbo', n)
        # Reversed iteration — max-clean entry is last, so it renders on top.
        # color_override preserves the original logical colour per airfoil.
        afl_dict = {}
        rainbow_colors = {}
        for j, ix in enumerate(ixs):
            label = f'Rough. Factor={j / (n - 1):.2f}'
            afl_dict[label] = _kulfan_entry(K_upper[ix], K_lower[ix], te_gap)
            rainbow_colors[label] = _turbo(j)

        compare_airfoils(
            afl_dict, [re], turb_cases, ['neuralfoil'],
            figurePath=str(out_dir / 'rainbow_plot.png'),
            color_override=rainbow_colors,
            cl_design=cl_design,
            reverse_plot_order=True,
        )

        # ── 2. Extreme-end comparison ────────────────────────────────────────
        print('  Building extreme-end comparison ...')
        x_par_full = _x_pareto(clean_LD, rough_LD)
        i_clean = int(np.argmin(np.abs(x_par_full - 0.0)))
        i_rough = int(np.argmin(np.abs(x_par_full - ROUGH_EXTREME)))
        _combined_endpoints[ctype] = {
            'i_clean': i_clean, 'i_rough': i_rough,
            'K_upper': K_upper, 'K_lower': K_lower,
            'out_dir': out_dir,
        }

        ext_dict = {
            r'Rough Bias = 0\%':                          _kulfan_entry(K_upper[i_clean], K_lower[i_clean], te_gap),
            rf'Rough Bias = {int(ROUGH_EXTREME*100)}\%':  _kulfan_entry(K_upper[i_rough], K_lower[i_rough], te_gap),
        }
        ext_dict.update(comps)

        plot2_color_cycle = ['#0065cc', '#eea800', '#009e73', '#d55e00', '#7860aa', '#ede13f', '#56b4ff', '#fca7c7', '#5d5d5d', '#000000']
        ext_colors = {name: plot2_color_cycle[i % len(plot2_color_cycle)]
                      for i, name in enumerate(ext_dict)}

        compare_airfoils(
            ext_dict, [re], turb_cases, ['neuralfoil'],
            figurePath=str(out_dir / 'extreme_comparison.png'),
            color_override=ext_colors,
            cl_design=cl_design,
        )

        # ── 3. Shape-function fitting ────────────────────────────────────────
        print('  Fitting shape functions ...')
        x_par    = x_par_full
        ixs_fit  = _arc_sample_indices(clean_LD, rough_LD, N_PARETO_SAMPLES)
        x_fit    = x_par[ixs_fit]

        # Compute Kulfan chord-wise stations (psi) once — same for all airfoils
        afl0 = Kulfan(TE_gap=te_gap)
        afl0.upperCoefficients = K_upper[0].tolist()
        afl0.lowerCoefficients = K_lower[0].tolist()
        afl0.chord = 1.0 * units.m
        psi_vals = np.array(afl0.psi)

        thickness_bases, camber_bases = [], []
        for ix in ixs_fit:
            afl = Kulfan(TE_gap=te_gap)
            afl.upperCoefficients = K_upper[ix].tolist()
            afl.lowerCoefficients = K_lower[ix].tolist()
            afl.chord = 1.0 * units.m
            thickness_bases.append(np.array([
                afl.computeThickness(1 - p).to('m').magnitude - te_gap * p
                for p in psi_vals
            ]))
            camber_bases.append(np.array(afl.yCamberLine_nondimensional))

        print('    fitting thickness profiles ...')
        thickness_params = np.array([fit_profile(psi_vals, t) for t in thickness_bases])
        print('    fitting camber profiles ...')
        camber_params    = np.array([fit_profile(psi_vals, c) for c in camber_bases])

        thickness_poly = [np.polyfit(x_fit, thickness_params[:, i], POLY_DEG)
                          for i in range(thickness_params.shape[1])]
        camber_poly    = [np.polyfit(x_fit, camber_params[:, i],    POLY_DEG)
                          for i in range(camber_params.shape[1])]

        shape_out = {
            'thickness': _poly_coeffs_dict(thickness_poly),
            'camber':    _poly_coeffs_dict(camber_poly),
            'meta': {
                'tau': tau, 'constraint_type': ctype,
                'TE_gap': te_gap, 'CBY_ORDER': CBY_ORDER, 'POLY_DEG': POLY_DEG,
                'N_PARETO_SAMPLES': N_PARETO_SAMPLES,
            },
        }
        sf_path = out_dir / 'shape_functions.json'
        with open(sf_path, 'w') as f:
            json.dump(shape_out, f, indent=4)
        print(f'    shape functions → {sf_path.relative_to(path_to_oso)}')

        if not PRODUCE_DIAG_PLOTS:
            continue  # skip diagnostic plots; _combined_endpoints already stored

        # ── Diagnostic Plot 1: Thickness profile fit quality ─────────────────
        fit_cmap = plt.get_cmap('turbo', len(thickness_bases))
        fig, ax = plt.subplots(figsize=(14, 5))
        for i, (t, params) in enumerate(zip(thickness_bases, thickness_params)):
            ax.plot(psi_vals, t, '-', color=fit_cmap(i), alpha=0.5, lw=1.0)
            ax.plot(psi_vals, shape_fcn(psi_vals, params), 'k-', alpha=0.2, lw=0.8)
        from matplotlib.lines import Line2D
        ax.legend(handles=[
            Line2D([0], [0], color=fit_cmap(0.5),
                   label=f'Data ({N_PARETO_SAMPLES} sampled, colour = Pareto position)'),
            Line2D([0], [0], color='k', label='shape_fcn fit'),
        ], loc='upper right')
        ax.set_xlabel(r'$\psi$  $(x/c)$')
        ax.set_ylabel('Thickness  (normalised chord)')
        ax.set_title(rf'$\tau$={tau}% {ctype} — thickness profile fit  (CBY_ORDER={CBY_ORDER})')
        ax.grid(True)
        fig.tight_layout()
        fig.savefig(str(out_dir / 'diag_1_thickness_fit.png'), dpi=150)
        plt.close(fig)

        # ── Diagnostic Plot 2: Camber profile fit quality ────────────────────
        fig, ax = plt.subplots(figsize=(14, 5))
        for i, (c, params) in enumerate(zip(camber_bases, camber_params)):
            ax.plot(psi_vals, c, '-', color=fit_cmap(i), alpha=0.5, lw=1.0)
            ax.plot(psi_vals, shape_fcn(psi_vals, params), 'k-', alpha=0.2, lw=0.8)
        ax.legend(handles=[
            Line2D([0], [0], color=fit_cmap(0.5),
                   label=f'Data ({N_PARETO_SAMPLES} sampled, colour = Pareto position)'),
            Line2D([0], [0], color='k', label='shape_fcn fit'),
        ], loc='upper right')
        ax.set_xlabel(r'$\psi$  $(x/c)$')
        ax.set_ylabel('Camber  (normalised chord)')
        ax.set_title(rf'$\tau$={tau}% {ctype} — camber profile fit  (CBY_ORDER={CBY_ORDER})')
        ax.grid(True)
        fig.tight_layout()
        fig.savefig(str(out_dir / 'diag_2_camber_fit.png'), dpi=150)
        plt.close(fig)

        # ── Diagnostic Plot 3: Thickness shape parameters vs x ───────────────
        n_params   = thickness_params.shape[1]
        param_cmap = plt.get_cmap('tab20', 20)
        x_dense    = np.linspace(0, 1, 300)
        fig, ax = plt.subplots(figsize=(12, 12))
        for i in range(n_params):
            label = 'r' if i == 0 else f'cby_{i-1}'
            col   = param_cmap(i % 20)
            ax.plot(x_fit, thickness_params[:, i], 'o', color=col, ms=2.5, alpha=0.6)
            ax.plot(x_dense, np.polyval(thickness_poly[i], x_dense), '-',
                    color=col, lw=1.2, label=label)
        ax.set_xlabel('x  (Pareto parameter:  0 = max clean L/D,  1 = max rough L/D)')
        ax.set_ylabel('Shape parameter value')
        ax.set_title(rf'$\tau$={tau}% {ctype} — thickness params vs x  (degree-{POLY_DEG} fits)')
        ax.legend(loc='upper right', fontsize=7, ncol=3, title='parameter')
        ax.grid(True)
        fig.tight_layout()
        fig.savefig(str(out_dir / 'diag_3_thickness_params.png'), dpi=150)
        plt.close(fig)

        # ── Diagnostic Plot 4: Camber shape parameters vs x ─────────────────
        fig, ax = plt.subplots(figsize=(12, 12))
        for i in range(camber_params.shape[1]):
            label = 'r' if i == 0 else f'cby_{i-1}'
            col   = param_cmap(i % 20)
            ax.plot(x_fit, camber_params[:, i], 'o', color=col, ms=2.5, alpha=0.6)
            ax.plot(x_dense, np.polyval(camber_poly[i], x_dense), '-',
                    color=col, lw=1.2, label=label)
        ax.set_xlabel('x  (Pareto parameter:  0 = max clean L/D,  1 = max rough L/D)')
        ax.set_ylabel('Shape parameter value')
        ax.set_title(rf'$\tau$={tau}% {ctype} — camber params vs x  (degree-{POLY_DEG} fits)')
        ax.legend(loc='upper right', fontsize=7, ncol=3, title='parameter')
        ax.grid(True)
        fig.tight_layout()
        fig.savefig(str(out_dir / 'diag_4_camber_params.png'), dpi=150)
        plt.close(fig)

        # ── Diagnostic Plot 5: Reconstructed airfoil family ──────────────────
        xsamples   = np.linspace(0, 1, N_DIAG_PLOT)
        turbo_cmap = plt.get_cmap('turbo', N_DIAG_PLOT)
        fig, ax = plt.subplots(figsize=(12, 5))
        for j, x in enumerate(xsamples[::-1]):
            orig_j = N_DIAG_PLOT - 1 - j
            t_params = [np.polyval(p, x) for p in thickness_poly]
            c_params = [np.polyval(p, x) for p in camber_poly]
            thickness = shape_fcn(psi_vals, t_params)
            camber    = shape_fcn(psi_vals, c_params)
            tfit      = thickness + te_gap * psi_vals
            xcoords = np.concatenate([psi_vals[::-1], psi_vals[1:]])
            ycoords = np.concatenate([(camber + tfit/2)[::-1], (camber - tfit/2)[1:]])
            ax.plot(xcoords, ycoords, color=turbo_cmap(orig_j))
        ax.set_xlabel('x/c')
        ax.set_ylabel('y/c')
        ax.set_title(rf'$\tau$={tau}% {ctype} — airfoil family'
                     f'  ({N_DIAG_PLOT} samples;  0=max clean L/D → 1=max rough L/D)')
        ax.set_aspect('equal')
        ax.grid(True)
        sm = plt.cm.ScalarMappable(cmap='turbo', norm=plt.Normalize(0, 1))
        fig.colorbar(sm, ax=ax, label='x  (0 = max clean L/D,  1 = max rough L/D)')
        fig.tight_layout()
        fig.savefig(str(out_dir / 'diag_5_airfoil_family.png'), dpi=150)
        plt.close(fig)

        # ── Diagnostic Plot 6: Pareto front vs shape-function reconstruction ─
        print('    computing reconstructed-airfoil performance for diag 6 ...')
        xsamples_d6 = np.linspace(0, 1, N_DIAG_PLOT)
        recon_clean, recon_rough = [], []
        for _x in xsamples_d6:
            _c, _r = _max_LD_from_shape(_x, thickness_poly, camber_poly,
                                        te_gap, psi_vals, re, cl_design)
            recon_clean.append(_c)
            recon_rough.append(_r)
        recon_clean = np.array(recon_clean)
        recon_rough = np.array(recon_rough)

        fig, ax = plt.subplots(figsize=(8, 6))
        sc = ax.scatter(rough_LD, clean_LD, c=x_par, cmap='turbo', s=15,
                        zorder=2, alpha=0.6, label='Pareto data (file)')
        ax.scatter(rough_LD[i_clean], clean_LD[i_clean],
                   marker='*', s=200, color='#0065cc', zorder=5,
                   edgecolors='k', linewidths=0.5, label=r'Rough Bias = 0\%')
        ax.scatter(rough_LD[i_rough], clean_LD[i_rough],
                   marker='*', s=200, color='#d55e00', zorder=5,
                   edgecolors='k', linewidths=0.5,
                   label=rf'Rough Bias = {int(ROUGH_EXTREME*100)}\%')
        ax.scatter(recon_rough, recon_clean, c=xsamples_d6, cmap='turbo',
                   s=60, marker='D', edgecolors='k', linewidths=0.6,
                   zorder=6, label=f'Shape fit ({N_DIAG_PLOT} samples)')
        fig.colorbar(sc, ax=ax, label='x  (Pareto parameter)').solids.set_alpha(1)
        ax.set_xlabel('Rough L/D')
        ax.set_ylabel('Clean L/D')
        ax.set_title(rf'$\tau$={tau}% {ctype} — Pareto front vs shape-fit reconstruction')
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(str(out_dir / 'diag_6_pareto_front.png'), dpi=150)
        plt.close(fig)

        print(f'    diagnostic plots → {out_dir.relative_to(path_to_oso)}')

    # ── Combined extreme-end comparison (unconstrained + moment-constrained) ─
    _ct_abbrev = {'unconstrained': 'UC', 'moment_constrained': 'MC'}
    _comb_color_cycle = ['#0065cc', '#eea800', '#009e73', '#d55e00', '#7860aa',
                         '#ede13f', '#56b4ff', '#fca7c7', '#5d5d5d', '#000000']
    if len(_combined_endpoints) == 2:
        print(f'  Building combined extreme-end comparison ...')
        comb_dict = {}
        for ct, ep in _combined_endpoints.items():
            abbr = _ct_abbrev.get(ct, ct)
            comb_dict[rf'{abbr} Rough Bias = 0\%'] = _kulfan_entry(
                ep['K_upper'][ep['i_clean']], ep['K_lower'][ep['i_clean']], te_gap)
            comb_dict[rf'{abbr} Rough Bias = {int(ROUGH_EXTREME*100)}\%'] = _kulfan_entry(
                ep['K_upper'][ep['i_rough']], ep['K_lower'][ep['i_rough']], te_gap)
        comb_dict.update(comps)
        comb_colors = {name: _comb_color_cycle[i % len(_comb_color_cycle)]
                       for i, name in enumerate(comb_dict)}

        _ep_uc = _combined_endpoints['unconstrained']
        _ep_mc = _combined_endpoints['moment_constrained']
        _comb_fname = 'extreme_comparison_combined.png'
        compare_airfoils(
            comb_dict, [re], turb_cases, ['neuralfoil'],
            figurePath=str(_ep_uc['out_dir'] / _comb_fname),
            color_override=comb_colors,
            cl_design=cl_design,
            legend_ncols=2,
        )
        shutil.copy(_ep_uc['out_dir'] / _comb_fname,
                    _ep_mc['out_dir'] / _comb_fname)
        print(f'    combined comparison → tau_{tau}/*/{_comb_fname}')

print('\nDone.')

