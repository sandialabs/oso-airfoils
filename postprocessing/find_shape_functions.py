"""
find_shape_functions.py

Finds polynomial functions mapping the Pareto blending parameter
    x = 0  ->  maximum clean L/D airfoil
    x = 1  ->  maximum rough L/D airfoil
to the shape coefficients of the thickness and camber profiles,
enabling smooth interpolation across the Pareto front.
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import natsort
import json
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit
from kulfan import Kulfan, units
from wt_objective_nsga2 import core_fitness_function
import sys

# Working Example:
#    python find_shape_functions.py cases/cases_111_to_120/case_112/c112_t24_k16_n752_l14_e130__2026_05_09_16-32-2763

# ==============================================================
# CONFIGURATION  —  all discretization choices live here
# ==============================================================
if len(sys.argv) > 1:
    PATH_TO_DATA = sys.argv[1]
else:
    raise ValueError("Please provide the path to the data directory as a command-line argument.")

# PATH_TO_DATA      = 'c112_t21_k16_n752_l10_e10__2026_05_10_01-02-3638'
# PATH_TO_DATA      = 'c112_t24_k16_n752_l14_e130__2026_05_09_16-32-2763'
N_PARETO_SAMPLES  = 100   # airfoils sampled from the Pareto front for fitting
CBY_ORDER         = 20    # Chebyshev basis order for thickness / camber fits
POLY_DEG          = 4     # polynomial degree for each coefficient(x) fit

# Below here does not affect the fit, just plotting
EVALUATE_FIT      = True   # whether to re-evaluate fitness of the fitted airfoils
N_PLOT            = 21    # airfoils to plot and re-evaluate
OBJ1_IX           = 6     # fitness-output column index (x-axis, rough L/D)
OBJ2_IX           = 5     # fitness-output column index (y-axis, clean L/D)
# ==============================================================


# ---- Chebyshev basis ----

def cby(IFUN_in, S):
    # IFUN  = IFUN_in - 20 if IFUN_in >= 21 else IFUN_in
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


# ---- Profile model: semicircular LE/TE term + Chebyshev corrections ----

def shape_fcn(x, params):
    r = params[0]
    return r * np.sqrt(1 - (x - 1) ** 2) - r * x + evaluate_cby(x, params[1:])

def fit_profile(psi, values):
    p0 = [0.0] + [0.0] * CBY_ORDER
    params, _ = curve_fit(
        lambda x, r, *c: shape_fcn(x, [r] + list(c)),
        psi, values, p0=p0
    )
    return np.asarray(params)


# ---- Load Pareto data ----

files = natsort.natsorted(
    [f for f in os.listdir(PATH_TO_DATA) if '.json' in f and 'population' in f],
    alg=natsort.ns.IGNORECASE
)
data = json.load(open(os.path.join(PATH_TO_DATA, files[-1])))

pop = data['population']
pareto_points = sorted(
    [p for p in pop if p['pareto_index'] == 1],
    key=lambda p: p['LoD_rough_at_design']
)

N_k    = int(data['input_parameters']['N_k'])
TE_gap = data['input_parameters']['TE_gap']

rough_LD = np.array([p['LoD_rough_at_design'] for p in pareto_points])
clean_LD = np.array([p['LoD_clean_at_design']  for p in pareto_points])
all_K    = np.array([list(p['K_upper']) + list(p['K_lower']) for p in pareto_points])


# ---- Compute x_pareto: normalised arc-length parameter along the Pareto front ----

pts = np.stack([rough_LD, clean_LD], axis=1)
arc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
arc /= arc[-1] if arc[-1] > 0 else 1.0

i_clean, i_rough = int(np.argmax(clean_LD)), int(np.argmax(rough_LD))
if arc[i_clean] < arc[i_rough]:
    x_pareto = (arc - arc[i_clean]) / (arc[i_rough] - arc[i_clean])
else:
    x_pareto = 1.0 - (arc - arc[i_rough]) / (arc[i_clean] - arc[i_rough])


# ---- Sample N_PARETO_SAMPLES airfoils evenly along the normalised arc ----

x_norm = float(np.ptp(rough_LD)) or 1.0
y_norm = float(np.ptp(clean_LD)) or 1.0
rn, cn  = rough_LD / x_norm, clean_LD / y_norm
cum_arc = np.concatenate([[0.0], np.cumsum(np.sqrt(np.diff(rn) ** 2 + np.diff(cn) ** 2))])
target_arcs = np.linspace(0.0, cum_arc[-1], N_PARETO_SAMPLES)
ixs = np.array([int(np.argmin(np.abs(cum_arc - t))) for t in target_arcs])


# ---- Extract thickness and camber for each sampled airfoil ----

thickness_bases, camber_bases = [], []
for i in ixs:
    K   = all_K[i]
    afl = Kulfan(TE_gap=TE_gap)
    afl.upperCoefficients = K[:N_k // 2]
    afl.lowerCoefficients = K[N_k // 2:]
    afl.chord = 1.0 * units.m
    psi_vals  = np.array(afl.psi)
    thickness = np.array([afl.computeThickness(1 - p).to('m').magnitude - TE_gap * p
                          for p in psi_vals])
    camber    = np.array(afl.yCamberLine_nondimensional)
    thickness_bases.append(thickness)
    camber_bases.append(camber)


# ---- Fit each profile with shape_fcn ----
# Result shapes: (N_PARETO_SAMPLES, CBY_ORDER + 2)  — params are [r, cby_0, ..., cby_CBY_ORDER]

print("Fitting thickness profiles ...")
thickness_params = np.array([fit_profile(psi_vals, t) for t in thickness_bases])

print("Fitting camber profiles ...")
camber_params = np.array([fit_profile(psi_vals, c) for c in camber_bases])


# ---- Fit each shape parameter as a polynomial in x ----

x_fit = x_pareto[ixs]   # actual Pareto parameter value for each sampled airfoil

thickness_poly = [np.polyfit(x_fit, thickness_params[:, i], POLY_DEG)
                  for i in range(thickness_params.shape[1])]
camber_poly    = [np.polyfit(x_fit, camber_params[:, i], POLY_DEG)
                  for i in range(camber_params.shape[1])]


# ---- Print the coefficient functions ----

def _poly_str(coeffs):
    deg   = len(coeffs) - 1
    parts = []
    for j, v in enumerate(coeffs):
        power = deg - j
        parts.append(f"{v:+.6g}*x^{power}" if power > 0 else f"{v:+.6g}")
    return " ".join(parts)

print(f"\n{'='*65}")
print(f"Thickness shape parameters  (degree-{POLY_DEG} polynomials in x)")
print(f"  CBY_ORDER={CBY_ORDER};  params layout: [r, cby_0, ..., cby_{CBY_ORDER}]")
print(f"{'='*65}")
for i, c in enumerate(thickness_poly):
    label = "r" if i == 0 else f"cby_{i-1}"
    print(f"  t[{i:2d}] ({label:6s})(x) = {_poly_str(c)}")

print(f"\n{'='*65}")
print(f"Camber shape parameters  (degree-{POLY_DEG} polynomials in x)")
print(f"{'='*65}")
for i, c in enumerate(camber_poly):
    label = "r" if i == 0 else f"cby_{i-1}"
    print(f"  c[{i:2d}] ({label:6s})(x) = {_poly_str(c)}")

# ---- Write coefficients to json ----

def _poly_coeffs_dict(poly_list):
    labels = ['r'] + [f'cby_{i}' for i in range(len(poly_list) - 1)]
    return {label: poly_list[i].tolist() for i, label in enumerate(labels)}

shape_function_output = [
    {'thickness': _poly_coeffs_dict(thickness_poly)},
    {'camber':    _poly_coeffs_dict(camber_poly)},
    {'input_parameters': data['input_parameters']},
]

out_path = os.path.join(PATH_TO_DATA, 'shape_functions.json')
with open(out_path, 'w') as f:
    json.dump(shape_function_output, f, indent=4)
print(f"Shape function coefficients written to {out_path}")


# ---- Helper: reconstruct thickness and camber profiles at parameter x ----

def get_profiles_at(x):
    t_params = [np.polyval(p, x) for p in thickness_poly]
    c_params = [np.polyval(p, x) for p in camber_poly]
    return shape_fcn(psi_vals, t_params), shape_fcn(psi_vals, c_params)

if EVALUATE_FIT:
    # ---- Diagnostic Plot 1: Thickness profile fit quality ----

    fit_cmap = plt.get_cmap('turbo', len(thickness_bases))
    fig, ax = plt.subplots(figsize=(14, 5))
    for i, (t, params) in enumerate(zip(thickness_bases, thickness_params)):
        ax.plot(psi_vals, t, '-', color=fit_cmap(i), alpha=0.5, lw=1.0)
        ax.plot(psi_vals, shape_fcn(psi_vals, params), 'k-', alpha=0.2, lw=0.8)
    ax.legend(handles=[
        Line2D([0], [0], color=fit_cmap(0.5), label=f'Data  ({N_PARETO_SAMPLES} sampled airfoils, color = Pareto position)'),
        Line2D([0], [0], color='k', label='shape_fcn fit'),
    ], loc='upper right')
    ax.set_xlabel('ψ  (x/c)')
    ax.set_ylabel('Thickness  (normalised chord)')
    ax.set_title(f'Thickness profile fit quality  (CBY_ORDER={CBY_ORDER})')
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(os.path.join(PATH_TO_DATA, 'diag_1_thickness_fit.png'), dpi=150)


    # ---- Diagnostic Plot 2: Camber profile fit quality ----

    fig, ax = plt.subplots(figsize=(14, 5))
    for i, (c, params) in enumerate(zip(camber_bases, camber_params)):
        ax.plot(psi_vals, c, '-', color=fit_cmap(i), alpha=0.5, lw=1.0)
        ax.plot(psi_vals, shape_fcn(psi_vals, params), 'k-', alpha=0.2, lw=0.8)
    ax.legend(handles=[
        Line2D([0], [0], color=fit_cmap(0.5), label=f'Data  ({N_PARETO_SAMPLES} sampled airfoils, color = Pareto position)'),
        Line2D([0], [0], color='k', label='shape_fcn fit'),
    ], loc='upper right')
    ax.set_xlabel('ψ  (x/c)')
    ax.set_ylabel('Camber  (normalised chord)')
    ax.set_title(f'Camber profile fit quality  (CBY_ORDER={CBY_ORDER})')
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(os.path.join(PATH_TO_DATA, 'diag_2_camber_fit.png'), dpi=150)


    # ---- Diagnostic Plot 3: Thickness shape parameters vs x ----

    n_params   = thickness_params.shape[1]
    param_cmap = plt.get_cmap('tab20', 20)
    x_dense    = np.linspace(0, 1, 300)

    fig, ax = plt.subplots(figsize=(12, 12))
    for i in range(n_params):
        label = 'r' if i == 0 else f'cby_{i-1}'
        col = param_cmap(i % 20)
        ax.plot(x_fit, thickness_params[:, i], 'o', color=col, ms=2.5, alpha=0.6)
        ax.plot(x_dense, np.polyval(thickness_poly[i], x_dense), '-', color=col, lw=1.2, label=label)
    ax.set_xlabel('x  (Pareto parameter:  0 = max clean L/D,  1 = max rough L/D)')
    ax.set_ylabel('Shape parameter value')
    ax.set_title(f'Thickness shape parameters vs Pareto position  (degree-{POLY_DEG} poly fits)')
    ax.legend(loc='upper right', fontsize=7, ncol=3, title='parameter')
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(os.path.join(PATH_TO_DATA, 'diag_3_thickness_params.png'), dpi=150)


    # ---- Diagnostic Plot 4: Camber shape parameters vs x ----

    fig, ax = plt.subplots(figsize=(12, 12))
    for i in range(camber_params.shape[1]):
        label = 'r' if i == 0 else f'cby_{i-1}'
        col = param_cmap(i % 20)
        ax.plot(x_fit, camber_params[:, i], 'o', color=col, ms=2.5, alpha=0.6)
        ax.plot(x_dense, np.polyval(camber_poly[i], x_dense), '-', color=col, lw=1.2, label=label)
    ax.set_xlabel('x  (Pareto parameter:  0 = max clean L/D,  1 = max rough L/D)')
    ax.set_ylabel('Shape parameter value')
    ax.set_title(f'Camber shape parameters vs Pareto position  (degree-{POLY_DEG} poly fits)')
    ax.legend(loc='upper right', fontsize=7, ncol=3, title='parameter')
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(os.path.join(PATH_TO_DATA, 'diag_4_camber_params.png'), dpi=150)


    # ---- Diagnostic Plot 5: Airfoil family + optional fitness re-evaluation ----
    xsamples   = np.linspace(0, 1, N_PLOT)
    turbo_cmap = plt.get_cmap('turbo', N_PLOT)
    fitness_vals = []

    fig, ax = plt.subplots(figsize=(12, 5))
    for i, x in enumerate(xsamples):
        thickness, camber = get_profiles_at(x)
        tfit = thickness + TE_gap * psi_vals

        xcoords = np.concatenate([psi_vals[::-1], psi_vals[1:]])
        ycoords = np.concatenate([(camber + tfit / 2)[::-1], (camber - tfit / 2)[1:]])

        ax.plot(xcoords, ycoords, color=turbo_cmap(i))

        afl = Kulfan(TE_gap=TE_gap)
        afl.fit2coordinates(xcoords, ycoords)
        K_sample = np.concatenate([afl.upperCoefficients, afl.lowerCoefficients])

        ipt = {'params': data['input_parameters'], 'individual': K_sample, 'pid': i}
        fitness = core_fitness_function(ipt)
        fitness_vals.append(fitness)

    ax.set_xlabel('x/c')
    ax.set_ylabel('y/c')
    ax.set_title(f'Airfoil family  ({N_PLOT} samples;  x = 0: max clean L/D  →  x = 1: max rough L/D)')
    ax.set_aspect('equal')
    ax.grid(True)
    sm = plt.cm.ScalarMappable(cmap='turbo', norm=plt.Normalize(0, 1))
    fig.colorbar(sm, ax=ax, label='x  (0 = max clean L/D,  1 = max rough L/D)')
    fig.tight_layout()
    fig.savefig(os.path.join(PATH_TO_DATA, 'diag_5_airfoil_family.png'), dpi=150)


    # ---- Diagnostic Plot 6: Pareto front overlay ----

    fitness_vals = np.array(fitness_vals)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(rough_LD, clean_LD, 'ko-', ms=4, label='Pareto front')
    ax.plot(abs(fitness_vals[:, OBJ1_IX]), abs(fitness_vals[:, OBJ2_IX]),
            'ro', ms=6, label='Shape-function samples')
    ax.set_xlabel('Rough L/D')
    ax.set_ylabel('Clean L/D')
    ax.set_title('Pareto front — shape-function reconstruction')
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(os.path.join(PATH_TO_DATA, 'diag_6_pareto_reconstruction.png'), dpi=150)
    # plt.show()
