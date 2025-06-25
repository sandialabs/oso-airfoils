from kulfan import Kulfan

import subprocess
import warnings
import tempfile
import numpy as np
import pandas as pd
# from kulfan import Kulfan
import os
import sys
import math
import shutil
path_to_XFOIL = shutil.which('xfoil')
import pathlib
import random
import string
path_to_here = pathlib.Path(__file__).parent.resolve()

class FNM(object):
    def __init__(self,ldr,N=5):
        # ldr = '/tmp/t_'
        x = ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(N))
        self.name = ldr + x

def run_xfoil(mode, 
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
        max_iter=100):

    try:
        iter(val)
        is_iterable = True
    except TypeError:
        is_iterable = False

    if 'dcmania' in str(path_to_here):
        tfpre = '/gpfs/dcmania/tempfiles/t_'
    elif 'ahsieh' in str(path_to_here):
        tfpre = '/gpfs/ahsieh/tempfiles/t_'
        # tfpre = '/pscratch/ahsieh/tempfiles/tmp_'
    elif 'karch' in str(path_to_here).lower():
        # tfpre = 't_'
        if sys.platform == "linux" or sys.platform == "linux2":
            # linux
            tfpre = '/home/codykarcher/t_'
        else:
            # OS X
            assert(sys.platform == "darwin")
            tfpre = 't_'
        # elif sys.platform == "win32":
        # Windows...


    else:
        # Default to David's file system
        tfpre = '/gpfs/dcmania/tempfiles/t_'

    tempDatfile    = FNM(tfpre,5)
    tempPolarfile  = FNM(tfpre,5)
    tempStdoutFile = FNM(tfpre,5)
    tempExecFile   = FNM(tfpre,5)
    tempCpDatafile = FNM(tfpre,5)
    tempBlDatafile = FNM(tfpre,5)


    # make sure we dont accidently reuse
    if os.path.exists(tempDatfile.name):
        os.remove(tempDatfile.name)
    if os.path.exists(tempPolarfile.name):
        os.remove(tempPolarfile.name)
    if os.path.exists(tempStdoutFile.name):
        os.remove(tempStdoutFile.name)
    if os.path.exists(tempExecFile.name):
        os.remove(tempExecFile.name)
    if os.path.exists(tempCpDatafile.name):
        os.remove(tempCpDatafile.name)
    if os.path.exists(tempBlDatafile.name):
        os.remove(tempBlDatafile.name)

    numberOfPanels = N_panels

    mode = mode.lower()
    if mode == 'alpha':
        mode = 'alfa'

    if mode not in ['alfa','cl']:
        raise ValueError('Invalid input mode.  Must be one of: alfa, cl ')

    # Removed this to keep to only standard type inputs
    # if isinstance(airfoil, str):        
    #     if ('.dat' in airfoil) or ('.txt' in airfoil):
    #         if os.path.isfile(airfoil):
    #             topline = 'load ' + airfoil + ' \n afl \n'
    #         else:
    #             raise ValueError('Could not find airfoil to be read')
    #     ck1 = 'naca' == airfoil.lower()[0:4]
    #     ck2 = airfoil[-4:].isdigit()
    #     if ck1:
    #         if ck2 and (len(airfoil)!=8):
    #             afl = airfoil.split()
    #             airfoil = afl[0]+afl[1]
    #         if ck2 and (len(airfoil)==8):
    #             topline = airfoil + ' \n'
    #         else:
    #             raise ValueError('Could not parse the NACA 4 digit airfoil')
        
    # elif isinstance(airfoil, Kulfan):
    #     if os.path.isfile(defaultDatfile):
    #         os.remove(defaultDatfile)
    #     airfoil.write2file(defaultDatfile)
    #     topline = 'load ' + defaultDatfile + ' \n' + 'airfoil \n'
        
    # else:
    #     raise ValueError('Could not parse airfoil')

    airfoil = Kulfan(TE_gap=TE_gap)
    airfoil.upperCoefficients = upperKulfanCoefficients
    airfoil.lowerCoefficients = lowerKulfanCoefficients
    airfoil.write2file(tempDatfile.name)
    
    assert(os.path.isfile(tempDatfile.name))
    # print('hello cody\n'*20)
    # print(os.path.isfile(tempDatfile.name))
    # print(tempDatfile.name)
    # ft = open(tempDatfile.name,'r')
    # texttemp = ft.read()
    # ft.close()
    # print(texttemp)
    # shutil.copy(tempDatfile.name, 'temp.txt')

    topline = 'load ' + tempDatfile.name + ' \n' + 'airfoil \n'
    
    estr = ''
    estr += 'plop\n'
    estr += 'g\n'
    estr += '\n'
    estr += topline
    estr += 'ppar\n'
    estr += 'n %d\n'%(numberOfPanels)
    estr += '\n'
    estr += '\n'
    if flapLocation is not None:
        ck1 = flapLocation >= 0.0
        ck2 = flapLocation <= 1.0
        if ck1 and ck2:
            estr += 'gdes \n'
            estr += 'flap \n'
            estr += '%f \n'%(flapLocation)
            estr += '999 \n'
            estr += '0.5 \n'
            estr += '%f \n'%(flapDeflection)
            estr += 'x \n'
            estr += '\n'
        else:
            raise ValueError('Invalid flapLocation.  Must be between 0.0 and 1.0')
    estr += 'oper \n'
    estr += "iter %d\n" %(max_iter)
    #run inviscid first
    if is_iterable:
        if mode == 'alfa':
            estr += "alfa %.2f \n" %(val[0])
        if mode == 'cl':
            estr += "cl %.3f \n" %(val[0])        
    else:
        if mode == 'alfa':
            estr += "alfa %.2f \n" %(val)
        if mode == 'cl':
            estr += "cl %.3f \n" %(val)
    estr += 'visc \n'
    estr += "%.0f \n" %(float(Re))
    estr += "M \n"
    estr += "%.2f \n" %(M)
    if N_crit < 9.0:
        # try to pre-seed rough cases
        if is_iterable:
            if mode == 'alfa':
                estr += "alfa %.2f \n" %(val[0])
            if mode == 'cl':
                estr += "cl %.3f \n" %(val[0])        
        else:
            if mode == 'alfa':
                estr += "alfa %.2f \n" %(val)
            if mode == 'cl':
                estr += "cl %.3f \n" %(val)
    estr += 'vpar \n'
    estr += 'xtr \n'
    estr += '%f \n'%(xtp_u)
    estr += '%f \n'%(xtp_l)
    estr += 'n \n'
    estr += '%f \n'%(N_crit)
    estr += '\n'
    estr += 'pacc \n'
    estr += tempPolarfile.name + ' \n'    #estr += '\n'
    estr += '\n'

    if is_iterable:
        if mode == 'alfa':
            estr += "aseq %.2f %.2f %.2f \n" %(val[0],val[1],val[2])
        if mode == 'cl':
            estr += "cseq %.3f %.3f %.3f \n" %(val[0],val[1],val[2])    

        # estr += 'pwrt \n'
        # estr += tempPolarfile.name + ' \n'    
        estr += '\n'
        estr += 'q \n'

    else:
        # if mode == 'alfa':
        #     estr += "alfa %.2f \n" %(val)
        # if mode == 'cl':
        #     estr += "cl %.3f \n" %(val)
        # # estr += 'pwrt \n'
        # # estr += tempPolarfile.name + ' \n'
        # estr += '\n'
        # estr += 'q \n'
        if mode == 'alfa':
            estr += "alfa %.2f \n" %(val)
        if mode == 'cl':
            estr += "cl %.3f \n" %(val)
        # estr += 'pwrt \n'
        # estr += tempPolarfile.name + ' \n'
        estr += 'cpwr \n'
        estr += tempCpDatafile.name + '\n'
        estr += 'dump \n'
        estr += tempBlDatafile.name + '\n'
        estr += '\n'
        estr += 'q \n'



    exFile = open(tempExecFile.name,'w')
    exFile.write(estr)
    exFile.close()


    cmd = ''
    if sys.platform == "linux" or sys.platform == "linux2":
        # linux
        cmd += 'timeout %d '%(timelimit)
    else:
        # OS X
        assert(sys.platform == "darwin")
        cmd += 'timelimit -t%d '%(timelimit)
    # elif sys.platform == "win32":
    # Windows...


    cmd += path_to_XFOIL
    cmd += ' <' + tempExecFile.name
    if sys.platform == "linux" or sys.platform == "linux2":
        # linux
        cmd += ' >'+tempStdoutFile.name
    else:
        # OS X
        assert(sys.platform == "darwin")
        cmd += ' &>'+tempStdoutFile.name
    # print(estr)

    assert(os.path.isfile(tempDatfile.name))

    try:
        # stdout_val = subprocess.check_output(cmd, shell=True, timeout=5)
        subprocess.run(cmd, shell=True)
    except:
        # process failed or timed out, will be handled below as a normal failure
        # print( upperKulfanCoefficients, lowerKulfanCoefficients, val)
        pass
    
    if os.path.exists(tempDatfile.name):
        os.remove(tempDatfile.name)
    if os.path.exists(tempStdoutFile.name):
        os.remove(tempStdoutFile.name)
    if os.path.exists(tempExecFile.name):
        os.remove(tempExecFile.name)

    # try:
    with warnings.catch_warnings():
        # catch warning for empty file
        warnings.simplefilter('ignore')
        data = np.genfromtxt(tempPolarfile.name, skip_header=12)

    if os.path.exists(tempPolarfile.name):
        os.remove(tempPolarfile.name)

    if not is_iterable:
        try:
            alpha   = data[0]
            cl      = data[1]
            cd      = data[2]
            cdp     = data[3]
            cm      = data[4]
            xtr_top = data[5]
            xtr_bot = data[6]
            Reval   = Re
            Mval    = M
        except:
            if os.path.exists(tempCpDatafile.name):
                os.remove(tempCpDatafile.name)
            if os.path.exists(tempBlDatafile.name):
                os.remove(tempBlDatafile.name)

            return None
    else:
        alpha   = data[:,0]
        cl      = data[:,1]
        cd      = data[:,2]
        cdp     = data[:,3]
        cm      = data[:,4]
        xtr_top = data[:,5]
        xtr_bot = data[:,6]
        Reval   = Re
        Mval    = M

    if not is_iterable:
        cpData = pd.read_csv(tempCpDatafile.name, sep="\\s+",skiprows=1, names = ['x' , 'cp'])
        blData = pd.read_csv(tempBlDatafile.name, sep="\\s+",skiprows=1, names = ['s', 'x', 'y', 'Ue/Vinf', 'Dstar', 'Theta', 'Cf', 'H', 'H*', 'P', 'm', 'K', 'tau', 'Di'])

    if os.path.exists(tempCpDatafile.name):
        os.chmod(tempCpDatafile.name,777)
        # subprocess.call(['rm -f %s'%(tempCpDatafile.name)],shell=True)
        os.remove(tempCpDatafile.name)
    if os.path.exists(tempBlDatafile.name):
        os.chmod(tempBlDatafile.name,777)
        # subprocess.call(['rm -f %s'%(tempBlDatafile.name)],shell=True)
        os.remove(tempBlDatafile.name)

    res = {}
    res['cd'] = cd
    res['cl'] = cl
    res['alpha'] = alpha
    res['cm'] = cm
    res['xtr_top'] = xtr_top
    res['xtr_bot'] = xtr_bot
    res['Re'] = Reval
    res['M'] = Mval
    res['N_crit'] = N_crit
    res['N_panels'] = N_panels
    if not is_iterable:
        res['cp_data'] = cpData.to_dict('list')
        res['bl_data'] = blData.to_dict('list')
    else:
        res['cp_data'] = None
        res['bl_data'] = None

    return res


# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
import os
import json
import numpy as np
import natsort
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})
# colors = ['#0065cc', '#e69f00', '#009e73', '#d55e00', '#56b4ff', '#fca7c7', '#ede13f', '#666666', '#000000']
colors = ['#0065cc', '#e69f00', '#009e73', '#d55e00', '#56b4ff', '#fca7c7', '#ede13f', '#000000']
matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=colors)

# from ada.geometry.airfoils.kulfan import Kulfan
# from ada.analysis.apis.xfoil.run import run as run_xfoil
# PATH_TO_ADA = os.environ['PATH_TO_ADA']


from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def cprint(x):
    sys.stdout.flush()
    print(x)

# ==================================================================================================================================================================
# ==================================================================================================================================================================
def run_xfoil_cl(cl_design, K_upper, K_lower, Re, N_crit, xtp_u, xtp_l, nm, TE_gap):
    res_fastrun = []
    for alpha in [0,5,10,15,20,25]:
        res_q = run_xfoil('alfa', K_upper, K_lower, alpha, Re=Re, N_crit=N_crit, xtp_u=xtp_u, xtp_l=xtp_l, TE_gap = TE_gap)
        if res_q is not None:
            res_fastrun.append(res_q)
            
    cl_quick = [rs['cl'] for rs in res_fastrun]
    alpha_quick = [rs['alpha'] for rs in res_fastrun]

    clean_results = []
    try:
        if cl_design < max(cl_quick):
            slice_idx = [ i for i,vl in enumerate(cl_quick) if vl>cl_design ][0]
            if cl_quick[slice_idx] == cl_quick[slice_idx-1]:
                dist = 0.5
                alpha_design_guess = alpha_quick[slice_idx-1] + (alpha_quick[slice_idx] - alpha_quick[slice_idx-1])*dist
            else:
                dist = ( cl_design - cl_quick[slice_idx-1] )/( cl_quick[slice_idx] - cl_quick[slice_idx-1] )
                alpha_design_guess = alpha_quick[slice_idx-1] + (alpha_quick[slice_idx] - alpha_quick[slice_idx-1])*dist
            res_aguess = run_xfoil('alfa', K_upper, K_lower, round(alpha_design_guess,3), Re=Re, N_crit=N_crit, xtp_u=xtp_u, xtp_l=xtp_l, TE_gap = TE_gap)
            if res_aguess is None:
                res_aguess = run_xfoil('alfa', K_upper, K_lower, round(alpha_design_guess,0), Re=Re, N_crit=N_crit, xtp_u=xtp_u, xtp_l=xtp_l, TE_gap = TE_gap)
                
            if res_aguess is None:
                print(nm + ' failed to find cl, did not converge')
                return None
            else:
                if abs(res_aguess['cl'] - cl_design)/cl_design > 0.05:
                    print(nm + ' failed to find cl, error of %f percent'%(abs(res_aguess['cl'] - cl_design)/cl_design*100))
                return res_aguess
            # else:
                # raise ValueError('Could not find CL')
        else:
            print(nm + ' failed to find cl, could not reach CL target')
            return None    
    except:
        print(nm + ' failed to find cl, could not reach CL target')
        return None    

# ==================================================================================================================================================================

existing_airfoil_dictionary = {
'ffa-w2-152': 
{'upperCoefficients': [0.17046706, 0.33438668, 0.15476145, 0.49129222, 0.17236828, 0.2962321, 0.18079263, 0.0893167],
'lowerCoefficients': [-0.14403581, -0.06194855, -0.19882893, 0.01263156, -0.21452513, -0.02732511, -0.13605042, 0.01492792],
'TE_gap': 0.00188},
'riso-a-15': 
{'upperCoefficients': [0.19542248, 0.4385377, 0.08352634, 0.66159665, 0.02541908, 0.32095343, 0.12465246, 0.11836227],
'lowerCoefficients': [-0.13928135, -0.01753147, -0.16014031, 0.01432329, -0.07359144, -0.0155406, -0.08035804, -0.04604944],
'TE_gap': 0.009644524},
'riso-b-15': 
{'upperCoefficients': [0.24531305, 0.21994651, 0.40331951, 0.08076417, 0.33022497, 0.12332088, 0.30827469, 0.16545387],
'lowerCoefficients': [-0.2356706, -0.06474068, -0.16383772, -0.16611335, -0.05207295, 0.01732643, 0.16465917, 0.10002146],
'TE_gap': 0.004074484},
'riso-p-15': 
{'upperCoefficients': [0.19124502, 0.43690242, 0.01158136, 0.85104171, -0.21667235, 0.40084947, 0.07796756, 0.14215856],
'lowerCoefficients': [-0.14290952, -0.05288197, -0.06074014, -0.13223268, 0.00117351, -0.13872487, -0.0071308, 0.1300919],
'TE_gap': 0.00282367},
's832': 
{'upperCoefficients': [0.1804298, 0.33029002, 0.22285964, 0.37236884, 0.46684779, 0.29243071, 0.43696239, 0.23753289],
'lowerCoefficients': [-0.08649814, 0.05120422, -0.19977357, 0.14689791, -0.24027691, -0.00103672, -0.12835736, 0.09898059],
'TE_gap': 0.0},
's826': 
{'upperCoefficients': [0.16618958, 0.29949937, 0.1402013, 0.42962015, 0.20899093, 0.25881877, 0.33315757, 0.39479959],
'lowerCoefficients': [-0.09194337, -0.08663391, -0.0765494, -0.33676823, 0.28607875, -0.18343956, 0.21115393, 0.22809879],
'TE_gap': 0.0},
'du_96-w-180': 
{'upperCoefficients': [0.17132441, 0.3558899, 0.15175815, 0.54129694, 0.09731448, 0.32825879, 0.23592736, 0.18349312],
'lowerCoefficients': [-0.14156696, -0.19339457, -0.12116041, -0.23327684, -0.12674316, -0.18319159, -0.08597008, 0.01516801],
'TE_gap': 0.00366855},
'ffa-w1-182': 
{'upperCoefficients': [0.20435661, 0.37141468, 0.17705193, 0.55064778, 0.0599962, 0.38180807, 0.14854804, 0.13854278],
'lowerCoefficients': [-0.15779433, -0.08180556, -0.34221352, 0.09207449, -0.45960241, 0.09850959, -0.08437446, 0.04294713],
'TE_gap': 0.0023},
'riso-a-18': 
{'upperCoefficients': [0.21065866, 0.36735624, 0.51013838, 0.05613477, 0.46262089, 0.13142507, 0.14547885, 0.12575698],
'lowerCoefficients': [-0.15689524, -0.05324062, -0.3489733, 0.1677874, -0.29472909, 0.11515652, -0.12257037, 0.00518386],
'TE_gap': 0.00993387},
'riso-b-17': 
{'upperCoefficients': [0.2229299, 0.25426803, 0.43999591, 0.08449051, 0.37073876, 0.15799446, 0.21596679, 0.25953147],
'lowerCoefficients': [-0.20292798, -0.17178474, -0.23256713, -0.18859865, -0.08813462, 0.09519591, 0.08677883, 0.17355953],
'TE_gap': 0.0058782},
'riso-p-18': 
{'upperCoefficients': [0.20270115, 0.44145012, 0.03134601, 0.81277195, -0.12242179, 0.44762556, 0.11394261, 0.21182957],
'lowerCoefficients': [-0.15104553, -0.11876674, -0.25876503, 0.05840517, -0.41413542, 0.21234144, -0.11221318, 0.23652712],
'TE_gap': 0.002267674},
's831': 
{'upperCoefficients': [0.19623838, 0.30673618, 0.33586392, 0.23243316, 0.69550957, 0.21244309, 0.57631864, 0.3466177],
'lowerCoefficients': [-0.07476615, -0.04263616, -0.14854538, -0.00517126, -0.30145584, 0.02069739, -0.03785553, 0.19365549],
'TE_gap': 0.0},
's825': 
{'upperCoefficients': [0.17496738, 0.30223429, 0.19149421, 0.38296943, 0.19527197, 0.28126152, 0.31700931, 0.4130541],
'lowerCoefficients': [-0.13002976, -0.12605804, -0.21698546, -0.46345804, 0.48796559, -0.40427287, 0.32412901, 0.20135574],
'TE_gap': 0.0},
'du_93-w-210': 
{'upperCoefficients': [0.17056056, 0.36019579, 0.14141691, 0.56737107, 0.06090653, 0.38053423, 0.19433952, 0.24414723],
'lowerCoefficients': [-0.16037265, -0.24691292, -0.23247346, -0.43582632, 0.06615759, -0.39116214, 0.12431883, 0.05953002],
'TE_gap': 0.003578184},
'ffa-w3-211': 
{'upperCoefficients': [0.23103006, 0.34720348, 0.2770677, 0.39871009, 0.16817752, 0.3177788, 0.15222807, 0.18010526],
'lowerCoefficients': [-0.17721306, -0.26048573, -0.20480766, -0.34035543, -0.18660585, -0.01224376, 0.0354027, 0.1508669],
'TE_gap': 0.00262},
'riso-a-21': 
{'upperCoefficients': [0.222549, 0.35423949, 0.52515433, 0.05761425, 0.47309112, 0.12196221, 0.16939894, 0.12796911],
'lowerCoefficients': [-0.18412686, -0.1113751, -0.59581571, 0.44462883, -0.78087246, 0.34099662, -0.20115337, 0.10182691],
'TE_gap': 0.010461728},
'riso-b-20': 
{'upperCoefficients': [0.25629888, 0.24992053, 0.44032427, 0.08734482, 0.3986768, 0.12477883, 0.3027098, 0.20727769],
'lowerCoefficients': [-0.27407555, -0.16039775, -0.47561739, -0.10224843, -0.19849772, 0.12422928, 0.0613274, 0.20466287],
'TE_gap': 0.008610528},
'riso-p-20': 
{'upperCoefficients': [0.2118682, 0.43776444, 0.06667393, 0.74530664, -0.03134031, 0.36805299, 0.16690098, 0.1839657],
'lowerCoefficients': [-0.16950715, -0.20653738, -0.41255241, 0.17619447, -0.79786642, 0.47146851, -0.1602198, 0.16963011],
'TE_gap': 0.0017394314},
's830': 
{'upperCoefficients': [0.20622614, 0.34439616, 0.27612088, 0.44804446, 0.39599184, 0.29248272, 0.49611748, 0.42637478],
'lowerCoefficients': [-0.08580882, -0.19106679, 0.06665299, -0.67026441, 0.11214829, -0.11742164, 0.0743748, 0.25363073],
'TE_gap': 0.0},
'du_91-w2-250': 
{'upperCoefficients': [0.22707578, 0.34390523, 0.26285831, 0.51199376, 0.1179782, 0.45998177, 0.18321857, 0.31932369],
'lowerCoefficients': [-0.28848951, -0.31526777, -0.37374971, -0.2570035, -0.38817863, -0.09366499, 0.12409451, 0.20466426],
'TE_gap': 0.0057973},
'ffa-w3-241': 
{'upperCoefficients': [0.26517955, 0.39218777, 0.24667021, 0.43539181, 0.15680384, 0.34962669, 0.16962077, 0.21850923],
'lowerCoefficients': [-0.25331048, -0.30143585, -0.3745882, -0.2015282, -0.45826686, 0.08992511, 0.0137028, 0.16684328],
'TE_gap': 0.00751},
'riso-a-24': 
{'upperCoefficients': [0.22916327, 0.36073498, 0.48983333, 0.13658164, 0.39953185, 0.16104349, 0.18238808, 0.15072615],
'lowerCoefficients': [-0.19453833, -0.20834691, -0.70171489, 0.48293367, -1.11236938, 0.51129168, -0.24122432, 0.22263018],
'TE_gap': 0.011032722},
'riso-b-23': 
{'upperCoefficients': [0.30167158, 0.19998371, 0.51368575, 0.03197528, 0.41708613, 0.16979434, 0.29692835, 0.23479967],
'lowerCoefficients': [-0.33670975, -0.1830624, -0.65570523, -0.10914245, -0.22267713, 0.07159744, 0.09287951, 0.13321409],
'TE_gap': 0.007791394},
'riso-p-23': 
{'upperCoefficients': [0.22105215, 0.38808627, 0.20882467, 0.59056897, 0.10637502, 0.3592841, 0.22848124, 0.25960753],
'lowerCoefficients': [-0.18679469, -0.3338346, -0.40032488, 0.04356117, -0.79412381, 0.41758196, -0.16986101, 0.30845074],
'TE_gap': 0.009597524},
's818': 
{'upperCoefficients': [0.22401092, 0.3191254, 0.35581284, 0.30848478, 0.33209201, 0.24243417, 0.38626227, 0.42118404],
'lowerCoefficients': [-0.18477642, -0.29170421, -0.44008326, -0.36203787, 0.11621318, -0.26443022, 0.08156934, 0.24819197],
'TE_gap': 0.0},
's814': 
{'upperCoefficients': [0.1954169, 0.3438724, 0.20566787, 0.47293529, 0.16285933, 0.32655937, 0.32661484, 0.38999088],
'lowerCoefficients': [-0.22844864, -0.41627951, -0.5132439, -0.08666286, -0.22511942, -0.00627433, -0.03266278, 0.25656376],
'TE_gap': 0.0},
'ffa-w3-270': 
{'upperCoefficients': [0.29571597, 0.43713915, 0.24933277, 0.47923386, 0.16593216, 0.37595725, 0.1806449, 0.23672802],
'lowerCoefficients': [-3.04187807e-01, -3.42727022e-01, -4.67380118e-01, -1.66084204e-01, -5.92587313e-01, 1.22304903e-01, 4.33504563e-04, 1.71333412e-01],
'TE_gap': 0.01012},
'riso-a-27': 
{'upperCoefficients': [0.24226475, 0.37592652, 0.35904575, 0.45305472, 0.07907014, 0.34069977, 0.19374066, 0.24794396],
'lowerCoefficients': [-0.22993112, -0.32120385, -0.51870077, -0.09712376, -0.71019277, 0.21043503, -0.12246814, 0.19302776],
'TE_gap': 0.010406692},
'riso-b-29': 
{'upperCoefficients': [0.3442155, 0.31936724, 0.49944582, 0.0837826, 0.47112241, 0.15946402, 0.29888669, 0.19711854],
'lowerCoefficients': [-0.37603138, -0.24328508, -0.90249249, -0.07959874, -0.39515422, 0.0182801, 0.07614434, 0.19133817],
'TE_gap': 0.011238376},
's815': 
{'upperCoefficients': [0.22296344, 0.30994497, 0.31412238, 0.34752115, 0.3044849, 0.24031627, 0.38500159, 0.38890486],
'lowerCoefficients': [-0.25863791, -0.48562564, -0.52120836, -0.15775011, -0.25292706, -0.02946459, -0.05458873, 0.2667111],
'TE_gap': 0.0},
'du_97-w-300': 
{'upperCoefficients': [0.26292326, 0.38176862, 0.34195592, 0.36427133, 0.260255, 0.34810944, 0.23970753, 0.29341288],
'lowerCoefficients': [-0.30635303, -0.42500691, -0.55371461, -0.43209065, -0.15802568, -0.31049503, 0.03190666, 0.2254738],
'TE_gap': 0.016682226},
'ffa-w3-301': 
{'upperCoefficients': [0.38846653, 0.42325908, 0.40761954, 0.32920532, 0.35262585, 0.29733538, 0.27157092, 0.24179111],
'lowerCoefficients': [-0.33715157, -0.45842177, -0.3639502, -0.38521479, -0.42305812, -0.05659543, 0.0059405, 0.22271602],
'TE_gap': 0.01828},
'riso-a-30': 
{'upperCoefficients': [0.24173236, 0.39296891, 0.32859688, 0.51399612, -0.03118848, 0.39831339, 0.20069154, 0.28974544],
'lowerCoefficients': [-0.21382606, -0.53207268, -0.60086837, 0.00856713, -1.02031982, 0.41492653, -0.18602209, 0.34724866],
'TE_gap': 0.011153848},
'ffa-w3-332': 
{'upperCoefficients': [0.40962967, 0.50921291, 0.34945622, 0.44373682, 0.26843572, 0.36967314, 0.24260932, 0.27658748],
'lowerCoefficients': [-0.44413537, -0.4661261, -0.49785939, -0.30626936, -0.56117073, -0.03612219, -0.09422596, 0.29839339],
'TE_gap': 0.02644},
'ffa-w3-360': 
{'upperCoefficients': [0.49686419, 0.57602135, 0.31116176, 0.63218354, 0.16306294, 0.50488796, 0.23219413, 0.32681488],
'lowerCoefficients': [-0.57619124, -0.35996929, -0.69628566, -0.13195661, -0.75765729, 0.05810408, -0.14168466, 0.32004427],
'TE_gap': 0.02896},
'riso-b-35': 
{'upperCoefficients': [0.41484296, 0.38011185, 0.60810518, 0.08407285, 0.60060771, 0.14423866, 0.40713469, 0.23377322],
'lowerCoefficients': [-0.44143236, -0.32009589, -1.0376858, -0.15126386, -0.42531432, -0.03315246, 0.15509788, 0.13652292],
'TE_gap': 0.01139995}}


def generateAirfoil(aflString):
    aflString = aflString.replace(' ','')
    aflString = aflString.lower()
    if aflString[0:4].lower() == 'naca':
        nacaNumbers = aflString[4:]
        if len(nacaNumbers)!=4:
            return 'Error: invalid naca4'
        else:
            afl1 = Kulfan() #AirfoilGeometry() 
            afl1.naca4_like(int(nacaNumbers[0]), int(nacaNumbers[1]), int(nacaNumbers[2:]))
 
    else:
        if aflString in existing_airfoil_dictionary.keys():
            afl1 = Kulfan(TE_gap = existing_airfoil_dictionary[aflString]['TE_gap'])
            afl1.upperCoefficients = existing_airfoil_dictionary[aflString]['upperCoefficients']
            afl1.lowerCoefficients = existing_airfoil_dictionary[aflString]['lowerCoefficients']
        else:
            raise ValueError('Could not find airfoil')
    return afl1

# ==================================================================================================================================================================
# ==================================================================================================================================================================

if sys.platform == "linux" or sys.platform == "linux2":
    # linux
    path_to_cache = '/home/codykarcher/Dropbox/research/workstation/cached_data/'
else:
    # OS X
    assert(sys.platform == "darwin")
    path_to_cache = '/Users/codykarcher/Dropbox/research/windTurbines/postprocessing/cached_data/'
# elif sys.platform == "win32":
# Windows...

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

comparisonAirfoils = {
    '15' : [  None         , 'ffa-w2-152', 'riso-a-15', 'riso-b-15', 'riso-p-15', 's832', 's826' ],
    '18' : [ 'du_96-w-180' , 'ffa-w1-182', 'riso-a-18', 'riso-b-17', 'riso-p-18', 's831', 's825' ],
    '21' : [ 'du_93-w-210' , 'ffa-w3-211', 'riso-a-21', 'riso-b-20', 'riso-p-20', 's830',  None  ],
    '24' : [ 'du_91-w2-250', 'ffa-w3-241', 'riso-a-24', 'riso-b-23', 'riso-p-23', 's818', 's814' ],
    '27' : [ 'du_91-w2-250', 'ffa-w3-270', 'riso-a-27', 'riso-b-29',  None      ,  None , 's815' ],
    '30' : [ 'du_97-w-300' , 'ffa-w3-301', 'riso-a-30', 'riso-b-29',  None      ,  None ,  None  ],
    '33' : [  None         , 'ffa-w3-332',  None      ,  None      ,  None      ,  None ,  None  ],
    '36' : [  None         , 'ffa-w3-360',  None      , 'riso-b-35',  None      ,  None ,  None  ],
}

files = natsort.natsorted([f for f in os.listdir(path_to_here) if '.txt' in f], alg=natsort.ns.IGNORECASE)

filename = files[0]
filecode = filename.split('.')[0]
filecodes = filecode.split('_')

CL = None
re = None

for fc in filecodes:
    if 'p' not in fc:
        if 'c' in fc:
            case_number = int(fc[1:])
        if 't' in fc:
            # tau = float(fc[1:])/100
            tau = fc[1:]
        if 'k' in fc:
            # N_k = int(fc[1:])
            N_k = int(fc[1:])
        if 'n' in fc:
            N_pop = int(fc[1:])
        if 'l' in fc:
            CL = float(fc[1:])/10.0
        if 'r' in fc:
            rLD = float(fc[1:])
        if 'e' in fc:
            re = float(fc[1:]) * 1e6
        if 'g' in fc:
            generation = float(fc[1:])

if re is None:
    re = design_matrix[ [ dm[0] for dm in design_matrix ].index(int(tau)/100)][3]

if CL is None:
    CL = design_matrix[ [ dm[0] for dm in design_matrix ].index(int(tau)/100)][1]

geo_plots = True
per_plots = True
take_best = False
itermax = int(round(len(files)/100,0)*100) + 100

if len(sys.argv) > 1:
    if sys.argv[1] == '1':
        take_best = True

plts = geo_plots

Nk = int(N_k/2)
Re = re

# mpirun -n 188 python -m mpi4py postprocess.py
# or
# mpirun -n 8 python -m mpi4py postprocess.py 1 
# to take the best airfoil in the run

# ==================================================================================================================================================================
# ==================================================================================================================================================================

if rank == 0:
    bestCandidates = []

    comp_afls = comparisonAirfoils[tau]
    for gen in range(0,len(files)):
        pop = np.loadtxt(str(path_to_here) + os.sep + files[gen])
        for ind in range(0,len(pop)):
            afl = Kulfan(TE_gap=te_gap_lookup[tau])
            afl.upperCoefficients = pop[ind][0:Nk]
            afl.lowerCoefficients = pop[ind][Nk:2*Nk]
            rs = pop[ind][2*Nk:]

            if ind == 0:
                obv = pop[ind][2*Nk]
                bestCandidates.append(pop[ind])
            else: 
                alpha = 0.01

        for iii,ca in enumerate(comp_afls):
            if ca is not None:
                cafl = generateAirfoil(ca)

if rank == 0:
    if per_plots:
        bcs = bestCandidates

        # plot evolution of objective function
        plt.figure(figsize=(12,8))
        plt.plot([b[2*Nk] for b in bcs])
        if take_best:
            best_index = np.argmin([b[2*Nk] for b in bcs])
        else:
            best_index = -1
        plt.title('Tau : %s'%(tau))
        plt.ylabel('Objective Function')
        plt.xlabel('Generation')
        # plt.ylim([20,35])
        plt.xlim([0,itermax])
        plt.grid(1)
        plt.tight_layout()
        plt.savefig(str(path_to_here) + os.sep + 'Objective_Evolution.png', dpi=250)
        plt.close()

        # plot evolution of objective function zoomed
        plt.figure(figsize=(12,8))
        plt.plot([b[2*Nk] for b in bcs])
        plt.title('Tau : %s'%(tau))
        plt.ylabel('Objective Function')
        plt.xlabel('Generation')
        objfcn = [b[2*Nk] for b in bcs]
        zoom_center = round(objfcn[-1]/10,0)*10
        plt.ylim([zoom_center-20,zoom_center+20])
        plt.xlim([0,itermax])
        plt.grid(1)
        plt.tight_layout()
        plt.savefig(str(path_to_here) + os.sep + 'Objective_Evolution_zoomed.png', dpi=250)
        plt.close()

        # plot evolution of design variables
        plt.figure(figsize=(12,8))
        for i in range(0,2*Nk):
            last_v = bcs[-1][i]
            # print(last_v)
            if i <= Nk-1:
                plt.plot([b[i] for b in bcs], color = colors[0], alpha = 0.5)
            else:
                plt.plot([b[i] for b in bcs], color = colors[1], alpha = 0.5)
        plt.title('Tau : %s'%(tau))
        plt.ylabel('Variable Value')
        plt.xlabel('Generation')
        plt.xlim([0,itermax])
        plt.grid(1)
        plt.tight_layout()
        plt.savefig(str(path_to_here) + os.sep + 'Variable_Evolution.png', dpi=250)
        plt.close()

        afl_final = Kulfan(TE_gap=te_gap_lookup[tau])
        afl_final.upperCoefficients = bcs[best_index][0:Nk]
        afl_final.lowerCoefficients = bcs[best_index][Nk:2*Nk]
        afl_final.write2file(str(path_to_here) + os.sep + 'final_airfoil_%s_ix%d.dat'%(tau,best_index))
        afl_final.scaleThickness(float(tau)/100.0)

        # afl_final.upperCoefficients = [.2,.2,.2,.2,.2,.2,.2,.2]
        # afl_final.lowerCoefficients = [-.1,-.1,-.1,-.1,-.1,-.1,-.1,-.1]

        # plot cp comparisons to reference airfoils @ design CL, Clean
        plt.figure(figsize=(12,8))
        for iii,ca in enumerate(comp_afls):
            if ca is not None:
                cafl = generateAirfoil(ca)
                if os.path.isfile(path_to_cache + 'cp_clean_%s.json'%(ca)):
                    with open(path_to_cache + 'cp_clean_%s.json'%(ca), 'r') as f:
                        rs = json.load(f)
                else:
                    rs = run_xfoil_cl(CL, cafl.upperCoefficients, cafl.lowerCoefficients, Re, N_crit=9.0, xtp_u=1.0, xtp_l=1.0, nm=ca, TE_gap=(cafl.zetaUpper[-1]-cafl.zetaLower[-1]))
                    with open(path_to_cache + 'cp_clean_%s.json'%(ca), 'w') as f:
                        json.dump(rs, f)
                if rs is not None:
                    plt.plot(rs['cp_data']['x'], rs['cp_data']['cp'], color=colors[iii], label = ca, alpha=0.5)
                else:
                    print(ca + ': xfoil did not converge')   
            
        rs = run_xfoil('alfa',  afl_final.upperCoefficients, afl_final.lowerCoefficients, bcs[best_index][2*Nk+2],  Re, N_crit=9.0, xtp_u=1.0, xtp_l=1.0, TE_gap = te_gap_lookup[tau])
        if rs is None:
            rs = run_xfoil('alfa',  afl_final.upperCoefficients, afl_final.lowerCoefficients, round(bcs[best_index][2*Nk+2],2),  Re, N_crit=9.0, xtp_u=1.0, xtp_l=1.0, TE_gap = te_gap_lookup[tau])
        if rs is not None:
            plt.plot(rs['cp_data']['x'], rs['cp_data']['cp'], color='k', label = 'Custom Airfoil')
                
        plt.title('Tau : %s , Clean'%(tau))
        plt.ylabel('Cp')
        plt.xlabel('x')
        plt.xlim([-0.02,1.02])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper right')#, fontsize='10')
        plt.gca().invert_yaxis()
        plt.savefig(str(path_to_here) + os.sep + 'CP_comparison_clean.png', dpi=250)
        plt.close()
            
        # plot cp comparisons to reference airfoils @ design CL, Rough
        plt.figure(figsize=(12,8))
        for iii,ca in enumerate(comp_afls):
            if ca is not None:
                cafl = generateAirfoil(ca)
                if os.path.isfile(path_to_cache + 'cp_rough_%s.json'%(ca)):
                    with open(path_to_cache + 'cp_rough_%s.json'%(ca), 'r') as f:
                        rs = json.load(f)
                else:
                    rs = run_xfoil_cl(CL, cafl.upperCoefficients, cafl.lowerCoefficients, Re, N_crit=3.0, xtp_u=0.05, xtp_l=0.05, nm=ca, TE_gap=(cafl.zetaUpper[-1]-cafl.zetaLower[-1]))
                    with open(path_to_cache + 'cp_rough_%s.json'%(ca), 'w') as f:
                        json.dump(rs, f)
                if rs is not None:
                    plt.plot(rs['cp_data']['x'], rs['cp_data']['cp'], color=colors[iii], label = str(ca), alpha=0.5)
                else:
                    print(ca + ': xfoil did not converge')
        
        rs = run_xfoil('alfa',  afl_final.upperCoefficients, afl_final.lowerCoefficients, bcs[best_index][2*Nk+2],  Re, N_crit=3.0, xtp_u=0.05, xtp_l=0.05, TE_gap = te_gap_lookup[tau])
        if rs is None:
            rs = run_xfoil('alfa',  afl_final.upperCoefficients, afl_final.lowerCoefficients, round(bcs[best_index][2*Nk+2],2),  Re, N_crit=3.0, xtp_u=0.05, xtp_l=0.05, TE_gap = te_gap_lookup[tau])
        if rs is not None:
            plt.plot(rs['cp_data']['x'], rs['cp_data']['cp'], color='k', label = 'Custom Airfoil')
            
        plt.title('Tau : %s , Rough'%(tau))
        plt.ylabel('Cp')
        plt.xlabel('x')
        # plt.ylim([0,5])
        plt.xlim([-0.02,1.02])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper right')#, fontsize='10')
        plt.gca().invert_yaxis()
        plt.savefig(str(path_to_here) + os.sep + 'CP_comparison_rough.png', dpi=250)
        plt.close()

        # Sweep AoA clean
        if os.path.isfile(path_to_cache + 'dta_clean_t%s.json'%(tau)):
            with open(path_to_cache + 'dta_clean_t%s.json'%(tau), 'r') as f:
                dta_clean = json.load(f)
        else:
            dta_clean = {}
            for iii,ca in enumerate(comp_afls):
                if ca is not None:
                    dta_clean[ca] = []
                    cafl = generateAirfoil(ca)
                    for alpha in np.linspace(-30,30,31):
                        res = run_xfoil('alfa', cafl.upperCoefficients, cafl.lowerCoefficients, alpha, Re, N_crit=9.0, xtp_u=1.0, xtp_l=1.0, TE_gap = (cafl.zetaUpper[-1]-cafl.zetaLower[-1]))
                        if res is not None:
                            dta_clean[ca].append(res)
                else:
                    dta_clean[iii]=None
            with open(path_to_cache + 'dta_clean_t%s.json'%(tau), 'w') as f:
                json.dump(dta_clean, f)
                
        dta_clean['Custom Airfoil'] = []
        for alpha in np.linspace(-30,30,31):
            res = run_xfoil('alfa', afl_final.upperCoefficients, afl_final.lowerCoefficients, alpha,  Re, N_crit=9.0, xtp_u=1.0, xtp_l=1.0, TE_gap = te_gap_lookup[tau])
            if res is not None:
                dta_clean['Custom Airfoil'].append(res) 
        
        # Sweep AoA rough
        if os.path.isfile(path_to_cache + 'dta_rough_t%s.json'%(tau)):
            with open(path_to_cache + 'dta_rough_t%s.json'%(tau), 'r') as f:
                dta_rough = json.load(f)
        else:
            dta_rough = {}
            for iii,ca in enumerate(comp_afls):
                if ca is not None:
                    dta_rough[ca] = []
                    cafl = generateAirfoil(ca)
                    for alpha in np.linspace(-30,30,31):
                        res = run_xfoil('alfa', cafl.upperCoefficients, cafl.lowerCoefficients, alpha, Re, N_crit=3.0, xtp_u=0.05, xtp_l=0.05, TE_gap = (cafl.zetaUpper[-1]-cafl.zetaLower[-1]))
                        if res is not None:
                            dta_rough[ca].append(res)
                else:
                    dta_rough[iii] = None
            with open(path_to_cache + 'dta_rough_t%s.json'%(tau), 'w') as f:
                json.dump(dta_rough, f)
                
        dta_rough['Custom Airfoil'] = []
        for alpha in np.linspace(-30,30,31):
            res = run_xfoil('alfa', afl_final.upperCoefficients, afl_final.lowerCoefficients, alpha,  Re, N_crit=3.0, xtp_u=0.05, xtp_l=0.05, TE_gap = te_gap_lookup[tau])
            if res is not None:
                dta_rough['Custom Airfoil'].append(res) 
                    
        # Plot CL vs. Alpha Curves
        plt.figure(figsize=(12,8))
        iii = 0
        for ky, val in dta_clean.items():
            if ky == 'Custom Airfoil':
                alpha_plot = 1.0
            else:
                alpha_plot = 0.5
            if val is not None:
                plt.plot([v['alpha'] for v in val], [v['cl'] for v in val], color = colors[iii], alpha = alpha_plot, label=str(ky))
            iii += 1
        plt.title('Tau : %s , Clean'%(tau))
        plt.ylabel('CL')
        plt.xlabel('Alpha')
        plt.xlim([-31,31])
        plt.ylim([-3,3])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper left')
        plt.savefig(str(path_to_here) + os.sep + 'CL_v_Alpha_clean.png', dpi=250)
        plt.close()
        
        plt.figure(figsize=(12,8))
        iii = 0
        for ky, val in dta_rough.items():
            if ky == 'Custom Airfoil':
                alpha_plot = 1.0
            else:
                alpha_plot = 0.5
            if val is not None:
                plt.plot([v['alpha'] for v in val], [v['cl'] for v in val], color = colors[iii], alpha = alpha_plot, label=str(ky))
            iii += 1
        plt.title('Tau : %s , Rough'%(tau))
        plt.ylabel('CL')
        plt.xlabel('Alpha')
        plt.xlim([-31,31])
        plt.ylim([-3,3])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper left')
        plt.savefig(str(path_to_here) + os.sep + 'CL_v_Alpha_rough.png', dpi=250)
        plt.close()

        
        # Plot CL vs. CD curves
        plt.figure(figsize=(12,8))
        iii = 0
        for ky, val in dta_clean.items():
            if ky == 'Custom Airfoil':
                alpha_plot = 1.0
            else:
                alpha_plot = 0.5
            if val is not None:
                plt.plot([v['cd'] for v in val], [v['cl'] for v in val], color = colors[iii], alpha = alpha_plot, label=str(ky))
            iii += 1
        plt.title('Tau : %s , Clean'%(tau))
        plt.ylabel('CL')
        plt.xlabel('CD')
        plt.xlim([0,0.03])
        plt.ylim([-3,3])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper right', framealpha=1)
        plt.savefig(str(path_to_here) + os.sep + 'CL_v_CD_clean.png', dpi=250)
        plt.close()
        
        plt.figure(figsize=(12,8))
        iii = 0
        for ky, val in dta_rough.items():
            if ky == 'Custom Airfoil':
                alpha_plot = 1.0
            else:
                alpha_plot = 0.5
            if val is not None:
                plt.plot([v['cd'] for v in val], [v['cl'] for v in val], color = colors[iii], alpha = alpha_plot, label=str(ky))
            iii += 1
        plt.title('Tau : %s , Rough'%(tau))
        plt.ylabel('CL')
        plt.xlabel('CD')
        plt.xlim([0,0.03])
        plt.ylim([-3,3])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper right', framealpha=1)
        plt.savefig(str(path_to_here) + os.sep + 'CL_v_CD_rough.png', dpi=250)
        plt.close()


        # Plot L/D vs Alpha curves
        plt.figure(figsize=(12,8))
        iii = 0
        for ky, val in dta_clean.items():
            if ky == 'Custom Airfoil':
                alpha_plot = 1.0
            else:
                alpha_plot = 0.5
            if val is not None:
                plt.plot([v['alpha'] for v in val], [v['cl']/v['cd'] for v in val], color = colors[iii], alpha = alpha_plot, label=str(ky))
            iii += 1
        plt.title('Tau : %s , Clean'%(tau))
        plt.ylabel('L/D')
        plt.xlabel('Alpha')
        plt.xlim([-31,31])
        plt.ylim([-300,300])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper left')
        plt.savefig(str(path_to_here) + os.sep + 'LD_v_Alpha_clean.png', dpi=250)
        plt.close()
        
        plt.figure(figsize=(12,8))
        iii = 0
        for ky, val in dta_rough.items():
            if ky == 'Custom Airfoil':
                alpha_plot = 1.0
            else:
                alpha_plot = 0.5
            if val is not None:
                plt.plot([v['alpha'] for v in val], [v['cl']/v['cd'] for v in val], color = colors[iii], alpha = alpha_plot, label=str(ky))
            iii += 1
        plt.title('Tau : %s , Rough'%(tau))
        plt.ylabel('L/D')
        plt.xlabel('Alpha')
        plt.xlim([-31,31])
        plt.ylim([-300,300])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper left')
        plt.savefig(str(path_to_here) + os.sep + 'LD_v_Alpha_rough.png', dpi=250)
        plt.close()


        # Plot L/D vs CL curves
        plt.figure(figsize=(12,8))
        iii = 0
        for ky, val in dta_clean.items():
            if ky == 'Custom Airfoil':
                alpha_plot = 1.0
            else:
                alpha_plot = 0.5
            if val is not None:
                plt.plot([v['cl'] for v in val], [v['cl']/v['cd'] for v in val], color = colors[iii], alpha = alpha_plot, label=str(ky))
            iii += 1
        plt.title('Tau : %s , Clean'%(tau))
        plt.ylabel('L/D')
        plt.xlabel('CL')
        plt.xlim([-3,3])
        plt.ylim([-300,300])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper left')
        plt.plot([CL,CL],[-300,300], color='k', linestyle='--', lw = 0.5)
        plt.savefig(str(path_to_here) + os.sep + 'LD_v_CL_clean.png', dpi=250)
        plt.close()
        
        plt.figure(figsize=(12,8))
        iii = 0
        for ky, val in dta_rough.items():
            if ky == 'Custom Airfoil':
                alpha_plot = 1.0
            else:
                alpha_plot = 0.5
            if val is not None:
                plt.plot([v['cl'] for v in val], [v['cl']/v['cd'] for v in val], color = colors[iii], alpha = alpha_plot, label=str(ky))
            iii += 1
        plt.title('Tau : %s , Rough'%(tau))
        plt.ylabel('L/D')
        plt.xlabel('CL')
        plt.xlim([-3,3])
        plt.ylim([-300,300])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper left')
        plt.plot([CL,CL],[-300,300], color='k', linestyle='--', lw = 0.5)
        plt.savefig(str(path_to_here) + os.sep + 'LD_v_CL_rough.png', dpi=250)
        plt.close()



if rank != 0:
# if rank == 0:
    comp_afls = comparisonAirfoils[tau]
    n_procs = size - 1

    for gen in range(0,len(files)):
        if gen%n_procs + 1 == rank:
            if plts:
                fig,ax=plt.subplots(1, figsize = (12,7)) #open a figure
            pop = np.loadtxt(str(path_to_here) + os.sep + files[gen])
            for ind in range(0,len(pop)):
                afl = Kulfan(TE_gap=te_gap_lookup[tau])
                afl.upperCoefficients = pop[ind][0:Nk]
                afl.lowerCoefficients = pop[ind][Nk:2*Nk]

                if ind == 0:
                    obv = pop[ind][2*Nk]
                else: 
                    alpha = 0.01
                    if plts:
                        plt.plot(afl.xcoordinates, afl.ycoordinates, color='k', alpha = alpha)

            for iii,ca in enumerate(comp_afls):
                if ca is not None:
                    cafl = generateAirfoil(ca)
                    if plts:
                        plt.plot(cafl.xcoordinates, cafl.ycoordinates, color=colors[iii], label = str(ca), alpha=0.5)

            if plts:
                ind = 0
                afl = Kulfan(TE_gap=te_gap_lookup[tau])
                if len(pop) >= 1:
                    afl.upperCoefficients = pop[ind][0:Nk]
                    afl.lowerCoefficients = pop[ind][Nk:2*Nk]
                    rs = pop[ind][2*Nk:]
                    plt.plot(afl.xcoordinates, afl.ycoordinates, color='k', label = 'Best Candidate')

                    plt.title('Thickness: 0.%s , Generation: %d, Objective: %.5f'%(tau,gen,obv))
                    plt.ylim([-0.25,0.25])
                    # plt.axis('equal')
                    ax.set_aspect('equal',adjustable='box')
                    plt.tight_layout()
                    plt.legend(loc='upper right', fontsize="10")
                    plt.savefig(str(path_to_here) + os.sep + 'Generation_%d_airfoils.png'%(gen), dpi=250)
                    plt.close()

if rank == 0:
    print('done')