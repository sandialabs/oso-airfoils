from wt_objective import airfoil_fitness

from ga_new_generation_mpi import newGeneration
import random
from kulfan import Kulfan
import numpy as np
import subprocess
import time
import os
import sys
import copy

import pathlib
path_to_here = pathlib.Path(__file__).parent.resolve()
filename = str(pathlib.Path(__file__).resolve()).split(os.sep)[-1]

from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def cprint(x):
    sys.stdout.flush()
    print(x)

# mpirun -n 8 python -m mpi4py mpitest.py

task_list = range(200)

for i, task in enumerate(task_list):
    if i%size==rank:
        cprint("Task %d (%d) on processor %d of %d"%(i,task,rank,size))
