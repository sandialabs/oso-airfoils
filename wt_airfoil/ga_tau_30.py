# from wt_objective_fast import airfoil_fitness
from wt_objective_updated import airfoil_fitness
from ga_new_generation_mpi import newGeneration
import random
from kulfan import Kulfan
import numpy as np
import subprocess
import time
import os
import sys

from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def cprint(x):
    sys.stdout.flush()
    print(x)

# mpirun -n 12 python -m mpi4py run_ga_mpi.py

N_k = 8
N_pop = 200
#N_pop = N_k*20
#N_pop = 100

# tau = 0.15
# tau = 0.18
# tau = 0.21
# tau = 0.24
# tau = 0.27
tau = 0.30
# tau = 0.33
# tau = 0.36


pop = None
if rank == 0:
    pop = [[random.uniform(-0.1,0.8) for j in range(0,int(N_k/2))]+[random.uniform(-0.8,0.1) for j in range(0,int(N_k/2))] for i in range(0,N_pop)]
    for i,p in enumerate(pop):
        afl = Kulfan(TE_gap = 0)
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
    folderstr = datestr + '__tau_%d'%(100*tau)
    ldr = '/gpfs/cjkarch/'
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
