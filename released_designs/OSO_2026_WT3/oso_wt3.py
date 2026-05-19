# load the jsons of polynomial data
import json
from pathlib import Path
import pathlib
# path_to_here = pathlib.Path.cwd()
import matplotlib.pyplot as plt
path_to_here = pathlib.Path(__file__).parent.resolve()
path_to_oso  = path_to_here.parent.parent   # released_designs/pareto_data → oso-airfoils/
import numpy as np
np.set_printoptions(linewidth=np.inf)
from kulfan import Kulfan


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

def oso_wt3(tau, rough_factor, N_points=100, moment_constrained = False):
    if rough_factor < 0 or rough_factor > 1:
        raise ValueError("rough_factor must be between 0 and 1")
    if tau not in [ 0.21, 0.24, 0.27, 0.30, 0.33, 0.36]:
        raise ValueError("tau must be one of [0.21, 0.24, 0.27, 0.30, 0.33, 0.36]")
    if rough_factor > 0.85:
        print("Warning: The OSO-2026-WT3 airfoils are known to have unusual leading edge shapes above a rough_factor of 0.85.  Use caution or consider reducing your selection of this value.")
    if rough_factor < 0.15 and moment_constrained:
        print("Warning: The OSO-2026-WT3 airfoils are known to have unusual lower surface shapes for clean-biased moment-constrained airfoils.  Use caution.")

    all_shape_functions = {}
    json_path = path_to_here / f"unconstrained/tau_{int(tau*100)}/shape_functions.json"
    with open(json_path, "r") as f:
        shape_functions = json.load(f)
    all_shape_functions[int(tau*100)] = shape_functions

    all_shape_functions_mc = {}
    json_path = path_to_here / f"moment_constrained/tau_{int(tau*100)}/shape_functions.json"
    with open(json_path, "r") as f:
        shape_functions = json.load(f)
    all_shape_functions_mc[int(tau*100)] = shape_functions
    
    if moment_constrained:
        all_shape_functions = all_shape_functions_mc

    # TE gap values for different tau values)
    te_gap_lookup = {
        '15':  0.00196,
        '18':  0.00230,
        '21':  0.00262,
        '24':  0.00751,
        '27':  0.01012,
        '30':  0.01140,
        '33':  0.01140,
        '36':  0.01140,
    }

    thickness_coeffs_all = all_shape_functions[int(tau*100)]['thickness']
    camber_coeffs_all = all_shape_functions[int(tau*100)]['camber']
    thickness_coeffs = []
    camber_coeffs = []
    ctr = 0
    for ky, vl in thickness_coeffs_all.items():
        evlt = 0
        for i, coeff in enumerate(vl):
            evlt += coeff * rough_factor**(len(vl)-1-i)
        if ctr == 0:
            assert(ky == 'r')
        else:
            assert(ky == 'cby_%d'%(ctr-1))
        thickness_coeffs.append(evlt)
        ctr += 1
    ctr = 0
    for ky, vl in camber_coeffs_all.items():
        evlt = 0
        for i, coeff in enumerate(vl):
            evlt += coeff * rough_factor**(len(vl)-1-i)
        if ctr == 0:
            assert(ky == 'r')
        else:
            assert(ky == 'cby_%d'%(ctr-1))
        camber_coeffs.append(evlt)
        ctr += 1

    ang = np.linspace(0,np.pi/2,N_points)
    psi = -np.cos(ang)+1

    thickness_profile = shape_fcn(psi, thickness_coeffs)
    camber_profile = shape_fcn(psi, camber_coeffs)
    xcoordinates = list(reversed(psi)) + list(psi[1:])
    y_upper = np.array(thickness_profile/2 + camber_profile) + psi*te_gap_lookup[str(int(tau*100))]/2
    y_lower = np.array(camber_profile - thickness_profile/2) - psi*te_gap_lookup[str(int(tau*100))]/2
    ycoordinates = list(reversed(y_upper)) + list(y_lower[1:])
    afl = Kulfan(TE_gap = te_gap_lookup[str(int(tau*100))])
    afl.fit2coordinates(np.array(xcoordinates), np.array(ycoordinates))

    return np.array([xcoordinates, ycoordinates]), afl


if __name__ == "__main__":
    plt.figure(figsize=(14, 5))
    for tau in [0.21, 0.24, 0.27, 0.30, 0.33, 0.36]:
        coords, afl = oso_wt3(tau, 0.8, moment_constrained=False)
        plt.plot(afl.xcoordinates, afl.ycoordinates)
    plt.axis('equal');
    plt.savefig("oso_wt3_unconstrained_080rough.png")