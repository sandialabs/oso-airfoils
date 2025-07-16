import subprocess
import warnings
import tempfile
import numpy as np
import pandas as pd
from kulfan import Kulfan
import os
import math
import shutil
path_to_XFOIL = shutil.which('xfoil')
import pathlib
path_to_here = pathlib.Path(__file__).parent.resolve()

def run(mode, 
        upperKulfanCoefficients,
        lowerKulfanCoefficients,
        val = 0.0, 
        Re = 1e7,
        M = 0.0,
        xtr_u=1.0,
        xtr_l=1.0,
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
        max_iter=200):

    if 'dcmania' in str(path_to_here):
        tfpre = '/gpfs/dcmania/tempfiles/tmp_'
    elif 'ahsieh' in str(path_to_here):
        tfpre = '/gpfs/ahsieh/tempfiles/tmp_'
    elif 'karch' in str(path_to_here):
        tfpre = '.'
    else:
        # Default to David's file system
        tfpre = '/gpfs/dcmania/tempfiles/tmp_'

    tempDatfile    = tempfile.NamedTemporaryFile(mode = 'w', prefix=tfpre, delete=False)
    tempPolarfile  = tempfile.NamedTemporaryFile(mode = 'w', prefix=tfpre, delete=True )
    # tempCpDatafile = tempfile.NamedTemporaryFile(mode = 'w', prefix=tfpre, delete=True )
    # tempBlDatafile = tempfile.NamedTemporaryFile(mode = 'w', prefix=tfpre, delete=True )
    # tempStdoutFile = tempfile.NamedTemporaryFile(mode = 'w', prefix=tfpre, delete=False)
    # tempStdoutFile = tempfile.NamedTemporaryFile(mode = 'w', prefix=tfpre, delete=True)
    tempExecFile   = tempfile.NamedTemporaryFile(mode = 'w', prefix=tfpre, delete=False)
    
    # need names, but not files
    # xfoil will require additional instructions if these files exist
    tempPolarfile.close()
    # tempCpDatafile.close()
    # tempBlDatafile.close()
    # tempStdoutFile.close()

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
    airfoil.write2file(tempDatfile)
    # print(airfoil.tau)
    
    # if not closed, xfoil flags as invalid file tye
    tempDatfile.close()
    
    topline = 'load ' + tempDatfile.name + ' \n' + 'airfoil \n'
    
    # proc = subprocess.Popen([path_to_XFOIL], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
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
    #run inviscid first
    if mode == 'alfa':
        estr += "alfa %.2f \n" %(val)
    if mode == 'cl':
        estr += "cl %.3f \n" %(val)
    estr += "iter %d\n" %(max_iter)
    estr += 'visc \n'
    estr += "%.0f \n" %(float(Re))
    estr += "M \n"
    estr += "%.2f \n" %(M)
    estr += 'vpar \n'
    estr += 'xtr \n'
    estr += '%f \n'%(xtr_u)
    estr += '%f \n'%(xtr_l)
    estr += 'n \n'
    estr += '%f \n'%(N_crit)
    estr += '\n'
    # step to alpha:
    # if mode == 'alfa':
    #     estr += "aseq %.2f %.2f %.2f \n" %(0.0, val, val/abs(val) * max([1.0,val/4]))
    estr += 'pacc \n'
    estr += '\n'
    estr += '\n'
    if mode == 'alfa':
        estr += "alfa %.2f \n" %(val)
    if mode == 'cl':
        estr += "cl %.3f \n" %(val)
    estr += 'pwrt \n'
    estr += tempPolarfile.name + ' \n'
    # estr += 'cpwr \n'
    # estr += tempCpDatafile.name + '\n'
    # estr += 'dump \n'
    # estr += tempBlDatafile.name + '\n'
    estr += '\n'
    estr += 'q \n'

    tempExecFile.write(estr)
    tempExecFile.close()

    if executionFile is not None:
        exFile = open(executionFile,'w')
        exFile.write(estr)
        exFile.close()
    
    # proc = subprocess.Popen([path_to_XFOIL], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    # proc.stdin.write(estr.encode())

    # try:
    #     stdout_val = proc.communicate(timeout=5)[0]
    #     tempStdoutFile.write(stdout_val.decode())
    # except TimeoutExpired:
    #     # process timed out but kept running, will be handled below as a normal failure
    #     proc.kill()
    # except:
    #     # process failed, will be handled below as a normal failure
    #     pass
    # proc.stdin.close()

    cmd = ''
    # cmd += 'timelimit -t2 '
    # cmd += 'timeout 2 '
    if 'dcmania' in str(path_to_here):
        cmd += 'timeout 2 '
    elif 'ahsieh' in str(path_to_here):
        cmd += 'timeout 2 '
    elif 'karch' in str(path_to_here):
        cmd += 'timelimit -t2 '
    else:
        # Default to linux
        cmd += 'timeout 2 '

    cmd += path_to_XFOIL
    cmd += ' <' + tempExecFile.name
    # cmd += ' >'+tempStdoutFile.name
    # print(estr)

    try:
        # stdout_val = subprocess.check_output(cmd, shell=True, timeout=5)
        subprocess.run(cmd, shell=True, timeout=5)
    except:
        # process failed or timed out, will be handled below as a normal failure
        # print( upperKulfanCoefficients, lowerKulfanCoefficients, val)
        pass
    
    try:
        with warnings.catch_warnings():
            # catch warning for empty file
            warnings.simplefilter('ignore')
            data = np.genfromtxt(tempPolarfile.name, skip_header=12)

        alpha   = data[0]
        cl      = data[1]
        cd      = data[2]
        cdp     = data[3]
        cm      = data[4]
        xtr_top = data[5]
        xtr_bot = data[6]
        Reval   = Re
        Mval    = M

        # cpData = pd.read_csv(tempCpDatafile.name, sep="\s+",skiprows=1, names = ['x' , 'cp'])

        # blData = pd.read_csv(tempBlDatafile.name, sep="\s+",skiprows=1, names = ['s', 'x', 'y', 'Ue/Vinf', 'Dstar', 'Theta', 'Cf', 'H', 'H*', 'P', 'm', 'K', 'tau', 'Di'])

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
        # res['cp_data'] = cpData.to_dict('list')
        # res['bl_data'] = blData.to_dict('list')
        # res['stdout'] = stdout_val.decode()

        # std_f = open(tempStdoutFile.name,'r')
        # std_vl = std_f.read()
        # std_f.close()
        # res['stdout'] = std_vl

    except:
        res = None

    if polarfile is not None:
        shutil.copy( tempPolarfile.name , polarfile)
    # if cpDatafile is not None:
    #     shutil.copy( tempCpDatafile.name , cpDatafile)
    # if blDatafile is not None:
    #     shutil.copy( tempBlDatafile.name , blDatafile)
    if defaultDatfile is not None:
        shutil.copy( tempDatfile.name , defaultDatfile)
    # if stdoutFile is not None:
    #     shutil.copy( tempStdoutFile.name , stdoutFile)

    tempDatfile.close()
    # tempStdoutFile.close()
    # tempExecFile.close()


    if os.path.exists(tempDatfile.name):
        os.remove(tempDatfile.name)
    if os.path.exists(tempPolarfile.name):
        os.remove(tempPolarfile.name)
    # if os.path.exists(tempCpDatafile.name):
    #     os.remove(tempCpDatafile.name)
    # if os.path.exists(tempBlDatafile.name):
    #     os.remove(tempBlDatafile.name)
    # if os.path.exists(tempStdoutFile.name):
    #     os.remove(tempStdoutFile.name)
    if os.path.exists(tempExecFile.name):
        os.remove(tempExecFile.name)


    # os.remove(tempDatfile.name)
    # # os.remove(tempPolarfile.name)
    # # os.remove(tempCpDatafile.name)
    # # os.remove(tempBlDatafile.name)
    # # os.remove(tempStdoutFile.name)
    # os.remove(tempExecFile.name)

    return res

if __name__ == '__main__':
    res = run('alpha',[0.2,0.2],[-0.2,-0.2],0)
    print('\n')
    print(res)
