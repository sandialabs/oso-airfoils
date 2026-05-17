import subprocess
import warnings
import tempfile
import numpy as np
import pandas as pd
from kulfan import Kulfan
import os
import sys
import math
import shutil
path_to_XFOIL = shutil.which('xfoil')
import pathlib
import random
import string
from datetime import datetime
path_to_here = pathlib.Path(__file__).parent.resolve()
import platform
from collections.abc import Iterable

import neuralfoil as nf  # `pip install neuralfoil`


def run(mode, 
        upperKulfanCoefficients,
        lowerKulfanCoefficients,
        val = 0.0, 
        Re = 1e7,
        M = 0.0,
        xtp_u=1.0,
        xtp_l=1.0,
        N_crit=9.0,
        N_panels = 160,
        flapLocation = None,
        flapDeflection = 0.0,
        polarfile      = None,
        cpDatafile     = None,
        blDatafile     = None,
        defaultDatfile = None,
        executionFile  = None,
        stdoutFile     = None,
        TE_gap = 0.0,
        timelimit = 10,
        max_iter=100,
        file_system = None,
        model = 'xxlarge'):

    if mode not in ['alpha','alfa']:
        raise ValueError('Neuralfoil must be operated in angle of attack mode')


    # check if alpha is a list or numpy array or tuple, etc
    if isinstance(val, Iterable) :
        assert(len(val)==3)
        # or isinstance(val, (str, bytes)):
        alpha = np.linspace(val[0], val[1], int((val[1]-val[0])/val[2])+1)
        # raise ValueError('Alpha must be an iterable (but not string or bytes)')
    # alpha = np.array([5,6,7])
    else:
        alpha = np.array([val])



    afl = Kulfan(TE_gap = TE_gap)
    afl.upperCoefficients = upperKulfanCoefficients
    afl.lowerCoefficients = lowerKulfanCoefficients
    afl.changeOrder(8)


    if model not in ["xsmall", "small", "medium", "large", "xlarge", "xxlarge", "xxxlarge"]:
        raise ValueError('Invalid Mode, must be one of ["xsmall", "small", "medium", "large", "xlarge", "xxlarge", "xxxlarge"]')


    res1 = nf.get_aero_from_kulfan_parameters(
        kulfan_parameters = {
            'upper_weights': afl.upperCoefficients.magnitude,
            'lower_weights': afl.lowerCoefficients.magnitude,
            'TE_thickness': TE_gap,
            'leading_edge_weight':0.0
        },
        alpha=alpha,
        Re=Re,
        n_crit = N_crit,
        xtr_upper = xtp_u,
        xtr_lower= xtp_l,
        model_size=model, 
    )

    stations = np.linspace(0,31,32, dtype=int)
    upper_vel_string = 'upper_bl_ue/vinf_'
    lower_vel_string = 'lower_bl_ue/vinf_'

    max_vel_ratio = 1.0
    cpmin_array = np.zeros(len(alpha))
    for ix, al in enumerate(alpha):
        for st in stations:
            vel_rat_upper = res1[upper_vel_string+str(st)][ix]
            vel_rat_lower = res1[lower_vel_string+str(st)][ix]
            if max([vel_rat_lower,vel_rat_upper]) > max_vel_ratio:
                max_vel_ratio = max([vel_rat_lower,vel_rat_upper])

        cpmin = 1 - max_vel_ratio**2
        cpmin_array[ix] = cpmin


    res = {}
    res['cd'] = res1['CD']
    res['cl'] = res1['CL']
    res['alpha'] = alpha
    res['cm'] = res1['CM']
    res['cpmin'] = np.array(cpmin_array)
    res['xtr_top'] = res1['Top_Xtr']
    res['xtr_bot'] = res1['Bot_Xtr']
    res['xtp_top'] = xtp_u
    res['xtp_bot'] = xtp_l
    res['Re'] = Re
    res['M'] = 0.0
    res['N_crit'] = N_crit
    res['N_panels'] = None

    return res

if __name__ == '__main__':
    res = run('alpha',[0.2,0.2],[-0.2,-0.2],0)
    print('\n')
    print(res)


    res = run('alpha',[0.1,0.1],[-0.1,-0.1],[0,-30,-0.25])
    print('\n')
    print(res)

    res = run('alpha',[0.1,0.1],[-0.1,-0.1],[0,-30,-0.25], file_system=1)
    print('\n')
    print(res)

    res = run('alpha',[0.1,0.1],[-0.1,-0.1],[0,-30,-0.25], file_system=2)
    print('\n')
    print(res)

    res = run('alpha',[0.1,0.1],[-0.1,-0.1],[0,-30,-0.25], file_system=3)
    print('\n')
    print(res)
