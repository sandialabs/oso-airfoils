import random
from kulfan import Kulfan
import numpy as np
import subprocess
import time
import os
import sys
import copy
import math

def newMember(N_k,tau):
    te_gap_lookup = {
        '15':  0.00196,
        '18':  0.00230,
        '21':  0.00262,
        '24':  0.00751,
        '27':  0.01012,
        '30':  0.01828,
        '33':  0.02644,
        '36':  0.02896,
    }

    design_matrix = [
        # tau,  CL,  spn,     Re
        [0.15, 1.5, 1.00, 10.0e6, ],
        [0.18, 1.5, 1.00, 10.0e6, ],
        [0.21, 1.5, 1.00, 12.0e6, ],
        [0.24, 1.4, 0.85, 13.0e6, ],
        [0.27, 1.3, 0.55, 16.0e6, ],
        [0.30, 1.2, 0.50, 18.0e6, ],
        [0.33, 1.2, 0.35, 16.0e6, ],
        [0.36, 1.2, 0.20, 13.0e6, ],
    ]

    ler_con = [0.007, 0.008, 0.01, 0.025, 0.03, 0.04, 0.06, 0.08]
    lev = ler_con[[dmr[0] for dmr in design_matrix].index(tau)]

    afl = Kulfan(TE_gap = te_gap_lookup[str(int(100*tau))])

    for i in range(0,1000):
        K = [random.uniform(-0.1,0.8) for j in range(0,int(N_k/2))]+[random.uniform(-0.8,0.5) for j in range(0,int(N_k/2))]
        
        Ku = K[0:int(len(K)/2)]
        Kl = K[int(len(K)/2):]
        # Kl[3] = .3
        afl.upperCoefficients = Ku
        afl.lowerCoefficients = Kl
        afl.scaleThickness(tau)
        tau_loc_u = afl.taumax_psi_upper
        tau_loc_l = afl.taumax_psi_lower
        hts = afl.getNormalizedHeight()
        leading_edge_radius_upper, leading_edge_radius_lower = afl.leadingEdgeRadius()

        delta_zeta_upper = afl.zetaUpper[1:] - afl.zetaUpper[0:-1]
        delta_psi = afl.psi[1:] - afl.psi[0:-1]
        first_derivative_approx = (delta_zeta_upper / delta_psi)
        delta_first_derivative = first_derivative_approx[1:] - first_derivative_approx[0:-1]
        delta_delta_psi = delta_psi[1:] - delta_psi[0:-1]
        second_derivative_approx = delta_first_derivative/delta_delta_psi
        positive_curvature = []
        for i in range(0, len(second_derivative_approx)):
            if second_derivative_approx[i] >0:
                positive_curvature.append(second_derivative_approx[i])
        curvature = sum(positive_curvature)

        delta_zeta_lower = afl.zetaLower[1:] - afl.zetaLower[0:-1]
        delta_psi = afl.psi[1:] - afl.psi[0:-1]
        first_derivative_approx_l = (delta_zeta_lower / delta_psi)
        delta_first_derivative_l = first_derivative_approx_l[1:] - first_derivative_approx_l[0:-1]
        delta_delta_psi = delta_psi[1:] - delta_psi[0:-1]
        second_derivative_approx_l = delta_first_derivative_l/delta_delta_psi
        sflips = 0
        sgn = second_derivative_approx_l[0]/abs(second_derivative_approx_l[0])
        for i in range(0, len(second_derivative_approx_l)):
            if second_derivative_approx_l[i]/abs(second_derivative_approx_l[i]) != sgn:
                sgn = second_derivative_approx_l[i]/abs(second_derivative_approx_l[i])
                sflips += 1

        K_candidate = afl.upperCoefficients.magnitude.tolist() + afl.lowerCoefficients.magnitude.tolist()

        conditions = [
            tau_loc_u >= 0.20, 
            tau_loc_u <= 0.40, 
            tau_loc_l >= 0.20, 
            tau_loc_l <= 0.40, 
            not any([math.isnan(rpv) for rpv in K_candidate]),
            max(abs(np.array(Ku)))<2.0, 
            max(abs(np.array(Kl)))<2.0, 
            not any(hts<0),
            leading_edge_radius_upper > .7*lev,
            leading_edge_radius_lower > .7*lev,
            Ku[0] > 0,
            Kl[0] < 0,
            curvature <= 300,
            sflips <= 1,
        ]

        if all(conditions):
            break

    return K_candidate

if __name__ == '__main__':
    N_k = 8
    tau = .30

    for i in range(0,10):
        K = newMember(N_k,tau)
        print(K)    

