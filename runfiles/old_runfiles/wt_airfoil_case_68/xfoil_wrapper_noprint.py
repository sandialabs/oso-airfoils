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
path_to_here = pathlib.Path(__file__).parent.resolve()

class FNM(object):
    def __init__(self,ldr,N=5):
        # ldr = '/tmp/t_'
        x = ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(N))
        self.name = ldr + x

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
        tfpre = 't_'
    else:
        # Default to local
        tfpre = 't_'

    tempDatfile    = FNM(tfpre,5)
    tempPolarfile  = FNM(tfpre,5)
    tempStdoutFile = FNM(tfpre,5)
    tempExecFile   = FNM(tfpre,5)

    # make sure we dont accidently reuse
    if os.path.exists(tempDatfile.name):
        os.remove(tempDatfile.name)
    if os.path.exists(tempPolarfile.name):
        os.remove(tempPolarfile.name)
    if os.path.exists(tempStdoutFile.name):
        os.remove(tempStdoutFile.name)
    if os.path.exists(tempExecFile.name):
        os.remove(tempExecFile.name)

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
        if mode == 'alfa':
            estr += "alfa %.2f \n" %(val)
        if mode == 'cl':
            estr += "cl %.3f \n" %(val)
        # estr += 'pwrt \n'
        # estr += tempPolarfile.name + ' \n'
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
    cmd += ' >'+tempStdoutFile.name
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
        alpha   = data[0]
        cl      = data[1]
        cd      = data[2]
        cdp     = data[3]
        cm      = data[4]
        xtr_top = data[5]
        xtr_bot = data[6]
        Reval   = Re
        Mval    = M

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

    return res

if __name__ == '__main__':
    res = run('alpha',[0.2,0.2],[-0.2,-0.2],0)
    print('\n')
    print(res)


    res = run('alpha',[0.1,0.1],[-0.1,-0.1],[0,-30,-0.25])
    print('\n')
    print(res)
