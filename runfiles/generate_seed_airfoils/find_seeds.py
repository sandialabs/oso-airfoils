import copy
import random
import numpy as np
import math
import sys

from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

from kulfan import Kulfan

def cprint(x):
    sys.stdout.flush()
    print(x)

N_k = 16
tau = 0.24
afl = Kulfan(TE_gap = 0.01)

# mpirun -n 188 python -m mpi4py postprocess.py

for ii in range(0,10000000):
    if ii % size == rank:
        K = [random.uniform(-0.1,0.8) for j in range(0,int(N_k/2))]+[random.uniform(-0.8,0.5) for j in range(0,int(N_k/2))]
        
        Ku = K[0:int(len(K)/2)]
        Kl = K[int(len(K)/2):]
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
            leading_edge_radius_upper > 0.005,
            leading_edge_radius_lower > 0.005,
            Ku[0] > 0,
            Kl[0] < 0,
            curvature <= 200,
            sflips <= 1,
        ]

        if all(conditions):
            f = open('viable_airfoils/afl_%d.txt'%(ii),'w')
            wst = ''
            for k in afl.upperCoefficients:
                wst += '%f,'%(k.to('dimensionless').magnitude)
            for k in afl.lowerCoefficients:
                wst += '%f,'%(k.to('dimensionless').magnitude)
            wst = wst[:-1] + '\n'
            f.write(wst)
            f.close()

            cprint('%d : %d : %s'%(ii, rank, wst[0:-2]))

# print(counter)
# return K_candidate

# k = 16
# nm = newMember(k)

# afl = Kulfan(TE_gap=0.01)
# afl.upperCoefficients = nm[0:int(k/2)]
# afl.lowerCoefficients = nm[int(k/2):k]
# afl.scaleThickness(0.21)
# import matplotlib.pyplot as plt
# %matplotlib inline

# plt.plot(afl.xcoordinates, afl.ycoordinates, label='Kulfan Airfoil', linewidth=2.0)
# plt.axis('equal')