# from wt_objective_fast import airfoil_fitness
# from wt_objective_updated import airfoil_fitness

# from wt_objective_equal import airfoil_fitness
# from wt_objective_4_1 import airfoil_fitness
# from wt_objective_fullRough import airfoil_fitness
# from wt_objective_equal_12 import airfoil_fitness
# from wt_objective_4_1_12 import airfoil_fitness
# from wt_objective_fullRough_12 import airfoil_fitness
# from wt_objective_equal_14 import airfoil_fitness
# from wt_objective_4_1_14 import airfoil_fitness
# from wt_objective_fullRough_14 import airfoil_fitness
from wt_objective_4_1_mingap import airfoil_fitness

from ga_new_generation_mpi import newGeneration
import random
from kulfan import Kulfan
import numpy as np
import subprocess
import time
import os
import sys

import pathlib
path_to_here = pathlib.Path(__file__).parent.resolve()

from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def cprint(x):
    sys.stdout.flush()
    print(x)

# mpirun -n 12 python -m mpi4py run_ga_mpi.py

# tau = 0.15
# tau = 0.18
# tau = 0.21
tau = 0.24
# tau = 0.27
# tau = 0.30
# tau = 0.33
# tau = 0.36

case_number = 10

N_k = 8
N_pop = 150

if 'dcmania' in str(path_to_here):
    ldr = '/gpfs/dcmania/'
elif 'ahsieh' in str(path_to_here):
    ldr = '/gpfs/ahsieh/'
elif 'karch' in str(path_to_here):
    ldr = './'
else:
    # Default to David's file system
    ldr = '/gpfs/dcmania/'

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

pop = None
if rank == 0:
    pop = [[random.uniform(-0.1,0.8) for j in range(0,int(N_k/2))]+[random.uniform(-0.8,0.1) for j in range(0,int(N_k/2))] for i in range(0,N_pop)]
    for i,p in enumerate(pop):
        afl = Kulfan(TE_gap = te_gap_lookup[str(int(100*tau))])
        # afl = Kulfan(TE_gap = 0)
        K = p
        Ku = K[0:int(len(K)/2)]
        Kl = K[int(len(K)/2):]
        afl.upperCoefficients = Ku
        afl.lowerCoefficients = Kl
        afl.scaleThickness(tau)
        pop[i] = afl.upperCoefficients.magnitude.tolist() + afl.lowerCoefficients.magnitude.tolist()

pop = comm.bcast(pop, root=0)
pop = newGeneration(airfoil_fitness, pop, normalizationVector = [1]*N_k, encodingTypes=[float]*N_k, lowerBounds=[-2.0]*N_k, upperBounds=[2.0]*N_k, tau=tau, initalize=True, comm=comm)
pop = comm.bcast(pop, root=0)

if rank == 0:
    datestr = time.strftime("%d_%b_%Y_%H-%M", time.localtime())
    folderstr = datestr + '__tau_%d'%(100*tau) + '__case_%d'%(case_number)
    # ldr = '/gpfs/cjkarch/'
    if not os.path.isdir(ldr + folderstr):
        os.mkdir(ldr + folderstr) 
    np.savetxt(ldr + folderstr + '/population_t%d_x%d_n%d_g0.txt'%(int(100*tau),N_k,N_pop),np.array(pop))
    # save population

for i in range(0,300):
    if rank == 0:
        cprint('Generation %d'%(i))
    pop = newGeneration(airfoil_fitness, pop, normalizationVector = [1]*N_k, encodingTypes=[float]*N_k, lowerBounds=[-2.0]*N_k, upperBounds=[2.0]*N_k, tau=tau, initalize=False, comm=comm)
    pop = comm.bcast(pop, root=0)
    if rank == 0:
        cprint(pop)
        np.savetxt(ldr + folderstr + '/population_t%d_x%d_n%d_g%d.txt'%(int(100*tau),N_k,N_pop,i),np.array(pop))
        # save population
