import os
import time
from multiprocessing import Pool
import datetime
import numpy as np
import matplotlib.pyplot as plt
PATH_TO_ADA = os.environ['PATH_TO_ADA']

# from ada.ui.uiManager import AirfoilGeometry
from ada.analysis.apis.xfoil import AnalysisCase
from ada.analysis.apis.xfoil.api_ui import generateXfoilResult
from ada.analysis.apis.xfoil.polarPlot import polarPlot

from ada.analysis.saveLoadAnalysisCase import saveAnalysisCase, loadAnalysisCase

import pandas as pd
from ada.geometry.airfoils.kulfan import Kulfan
from matplotlib.gridspec import GridSpec
# import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})
import matplotlib
colors = ['#0065cc', '#e69f00', '#009e73', '#d55e00', '#56b4ff', '#fca7c7', '#ede13f', '#666666', '#000000']
matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=colors)
from ada.analysis.apis.xfoil.helperFunctions import truncate_colormap, handleZeroDivide, get_colors, computeNormals, get_fractional_color

import natsort

from ada.geometry.airfoils.kulfan import units

from mpi4py import MPI
comm: MPI.Comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# if rank == 0:
#     startTime = time.time()


def generateAirfoil(aflString, useFits=False):
    aflString = aflString.replace(' ','')
    aflString = aflString.lower()
    if aflString[0:4].lower() == 'naca':
        nacaNumbers = aflString[4:]
        if len(nacaNumbers)!=4:
            return 'Error: invalid naca4'
        else:
            afl1 = Kulfan() 
            afl1.naca4_like(int(nacaNumbers[0]), int(nacaNumbers[1]), int(nacaNumbers[2:]))
            # uiManager.geometryItems.append(afl1)

    else:
        # this needs to utilize the large database of airfoils
        afl1 = Kulfan()

        if useFits:
            airfoil_directories = [
                PATH_TO_ADA + 'ada/geometry/airfoil_scraping/raw_airfoil_files/uiuc_airfoils',
                PATH_TO_ADA + 'ada/geometry/airfoil_scraping/raw_airfoil_files/airfoiltools_airfoils',
                PATH_TO_ADA + 'ada/geometry/airfoil_files',
            ]
        else:
            airfoil_directories = [
                '/Users/cjkarch/root/research/windTurbines',
                PATH_TO_ADA + 'ada/geometry/airfoil_scraping/raw_airfoil_files/uiuc_airfoils',
                PATH_TO_ADA + 'ada/geometry/airfoil_scraping/raw_airfoil_files/airfoiltools_airfoils',
                PATH_TO_ADA + 'ada/geometry/airfoil_files',
            ]

        airfoil_directory_data = {}
        for afl_dir in airfoil_directories:
            # file_data_dict = {}
            pth = afl_dir
            combined = os.listdir(pth)
            files   = [f for f in combined if os.path.isfile(pth + os.sep + f) and '.dat' in f]
            for fl in files:
                if fl not in list(airfoil_directory_data.keys()):
                    airfoil_directory_data[fl] = pth + os.sep + fl

        search_file = aflString + '.dat'
        afl1.fit2file(airfoil_directory_data[search_file])

    return afl1


pth = PATH_TO_ADA + 'ada/geometry/airfoil_files'
combined = os.listdir(pth)
files   = [f for f in combined if os.path.isfile(pth + os.sep + f) and '.dat' in f]

combined2 = os.listdir('.')
files2   = natsort.natsorted([f for f in combined2 if os.path.isfile(f) and 'du_' in f and '.dat' in f], alg=natsort.ns.IGNORECASE)

for fl2 in files2:
    if fl2 not in files:
        files.append(fl2)


maxDataSet = 2000
# maxDataSet = 100



airfoilDict = {}
for fl in files:
    airfoilName = fl[0:-4]
    airfoilDict[airfoilName] = {}
    airfoilDict[airfoilName]['filePath'] = pth + os.sep + fl
    airfoilDict[airfoilName]['geometry'] = generateAirfoil(airfoilName)
    
dt_string = datetime.datetime.now().strftime("%m-%d-%Y_%H-%M-%S")
rundataFolder = 'rundata_' + dt_string
try:    
    os.mkdir(rundataFolder)
except:
    pass

airfoilNameList_all = natsort.natsorted(list(airfoilDict.keys()), alg=natsort.ns.IGNORECASE)

# airfoilNameList = []
# for afl in airfoilNameList_all:
#     if afl[0] == 's':
#         if afl in ['s832','s831','s830','s818','s826','s825','s814','s815']:
#             airfoilNameList.append(afl)
#     else:
#         airfoilNameList.append(afl)

# airfoilNameList = [
#  # 'du_00-w2-350',
#  # 'du_00-w2-401',
#  'du_00-w-212',
#  'du_91-w2-250',
#  'du_93-w-210',
#  'du_95-w-180',
#  'du_96-w-180',
#  'du_97-w-300',
#  'ffa-w1-128',
#  'ffa-w1-152',
#  'ffa-w1-182',
#  'ffa-w1-211',
#  'ffa-w1-242',
#  'ffa-w1-271',
#  'ffa-w2-152',
#  'ffa-w2-210',
#  'ffa-w3-211',
#  'ffa-w3-241',
#  'ffa-w3-270',
#  'ffa-w3-301',
#  'ffa-w3-332',
#  'ffa-w3-360',
#  'riso-a-12',
#  'riso-a-15',
#  'riso-a-18',
#  'riso-a-21',
#  'riso-a-24',
#  'riso-a-27',
#  'riso-a-30',
#  'riso-b-15',
#  'riso-b-17',
#  'riso-b-20',
#  'riso-b-23',
#  'riso-b-29',
#  'riso-b-35',
#  'riso-p-15',
#  'riso-p-18',
#  'riso-p-20',
#  'riso-p-23',
#  's814',
#  's815',
#  's818',
#  's825',
#  's826',
#  's830',
#  's831',
#  's832',
#  ]


airfoilNameList = [
 # 's832',
 # 's831',
 # 's830',
 # 's818',
 # 'du_00-w2-350',
 # 'du_00-w2-401',
 # 'du_00-w-212',
 # 'du_91-w2-250',
 # 'du_93-w-210',
 # 'du_95-w-180',
 # 'du_96-w-180',
 # 'du_97-w-300',
 # 'riso-b-15',
 # 'riso-b-17',
 # 'riso-b-20',
 # 'riso-b-23',
 # 'riso-b-29',
 # 'riso-b-35',
 # 'ffa-w1-128',
 # 'ffa-w1-152',
 # 'ffa-w1-182',
 # 'ffa-w1-211',
 # 'ffa-w1-242',
 # 'ffa-w1-271',
 # 'ffa-w2-152',
 # 'ffa-w2-210',
 # 'ffa-w3-211',
 # 'ffa-w3-241',
 # 'ffa-w3-270',
 # 'ffa-w3-301',
 # 'ffa-w3-332',
 # 'ffa-w3-360',
 's826',
 's825',
 's814',
 's815',
 'riso-p-15',
 'riso-p-18',
 'riso-p-20',
 'riso-p-23',
 'riso-a-12',
 'riso-a-15',
 'riso-a-18',
 'riso-a-21',
 'riso-a-24',
 'riso-a-27',
 'riso-a-30',
 ]


if rank == 0:
    # logfileName = 'logfile_' + dt_string + '.txt'
    logfileName = 'logfile.txt'
    if os.path.isfile(logfileName):
        os.remove(logfileName)
    logfile = open(logfileName,'w')
    logfile.close()

if rank == 0:
    startTime = time.time()
    startTime_sub = time.time()

for jj,airfoilName in enumerate(airfoilNameList):
    afl = airfoilDict[airfoilName]['geometry']
    if rank == 0:
        logfile = open(logfileName,'a')
        logfile.write(str(jj+1) + ' of '+ str(len(airfoilNameList)) +  ' : ' + airfoilName + ' ...\n') 
        logfile.write('Total Time : ' + '%.2f\n'%(time.time() - startTime))
        logfile.write('Last Time  : ' + '%.2f\n'%(time.time() - startTime_sub))
        logfile.close()
        startTime_sub = time.time()
        
        # afl = airfoilDict[airfoilName]['geometry']
        
        airfoilDataFolder = rundataFolder + os.sep + airfoilName
        # print(airfoilDataFolder)
        try:
            os.mkdir(airfoilDataFolder)
        except:
            pass

        ac = AnalysisCase()
        # ac.inputs['alpha'] = np.linspace(-10,10,16).tolist()
        # ac.inputs['re']    = [1e5,1e6,1e7]
        # ac.inputs['alpha'] = np.linspace(-20,30,51).tolist()
        # ac.inputs['re']    = (10**np.linspace(np.log10(5e5),np.log10(30e6),5)).tolist()

        lft = np.linspace(-30,-20,11)[0:-1]
        cnr = np.linspace(-20,20,81)
        rgt = np.linspace(20,30,11)[1:]
        alphaList = np.append(lft,np.append(cnr,rgt))
        ac.inputs['alpha']     = alphaList.tolist()
        ac.inputs['re']        = [1e6, 3e6, 5e6, 1e7, 3e7]
        # ac.inputs['mach']      = np.linspace(0,0.5,6).tolist()
        # ac.inputs['n_panels']  = [160,180,200,220]
        ac.inputs['n_crit']    = [3,6,9]
        ac.inputs['xtr_u']     = np.linspace(0,1.0,11).tolist()
        ac.inputs['xtr_l']     = np.linspace(0,1.0,11).tolist()
        
        if 'cl' in list(ac.inputs.keys()):
            xfoil_mode = 'cl'
        else:
            xfoil_mode = 'alfa'

        rcs = ac.generateSingleCases()
        
        rcsPartitions = []
        
        rcsPart = []
        for i,rcc in enumerate(rcs):
            if len(rcsPart) < maxDataSet:
                rcsPart.append(rcc)
            else:
                rcsPartitions.append(rcsPart)
                rcsPart = []
                rcsPart.append(rcc)
                
        if len(rcsPart)>0:
            rcsPartitions.append(rcsPart)


        N_partitions = len(rcsPartitions)
    else:
        N_partitions = None
        rcsPartions = None
        xfoil_mode = None
        airfoilDataFolder = None

    N_partitions = comm.bcast(N_partitions,root=0)
    xfoil_mode = comm.bcast(xfoil_mode,root=0)
    airfoilDataFolder = comm.bcast(airfoilDataFolder,root=0)

    # for ii, rcs in enumerate(rcsPartitions):
    for ii in range(0,N_partitions):
        results = []
        rcDicts = []

        if rank == 0:
            rcs = rcsPartitions[ii]
        else:
            rcs = None
        rcs = comm.bcast(rcs,root=0)

        for i,rc in enumerate(rcs):
            if i%size == rank:
                rcDict = {}
                if xfoil_mode == 'cl':
                    rcDict['cl'] = rc.inputs['cl']
                else:
                    rcDict['alpha'] = rc.inputs['alpha']

                rcDict['fileLeader']              = airfoilDataFolder
                rcDict['upperKulfanCoefficients'] = [c.magnitude for c in afl.upperCoefficients]
                rcDict['lowerKulfanCoefficients'] = [c.magnitude for c in afl.lowerCoefficients]
                rcDict['process_id']              = i
                rcDict['re']                      = rc.inputs['re']
                rcDict['mach']                    = rc.inputs['mach']
                rcDict['n_panels']                = rc.inputs['n_panels']
                rcDict['n_crit']                  = rc.inputs['n_crit']
                rcDict['xtr_u']                   = rc.inputs['xtr_u']
                rcDict['xtr_l']                   = rc.inputs['xtr_l']

                resDict = generateXfoilResult(rcDict)

                if resDict is not None:
                    if xfoil_mode == 'alfa':
                        rc.outputs['cl']      = resDict['cl']
                    else:
                        rc.outputs['alpha']   = resDict['alpha']

                    rc.outputs['bl_data'] = resDict['bl_data']
                    rc.outputs['cd']      = resDict['cd']
                    rc.outputs['cm']      = resDict['cm']
                    rc.outputs['cp_data'] = resDict['cp_data']
                    rc.outputs['xtr_l']   = resDict['xtr_l']
                    rc.outputs['xtr_u']   = resDict['xtr_u']

                    rc.geometry = {}
                    rc.geometry['upperKulfanCoefficients'] = [c.magnitude for c in afl.upperCoefficients]
                    rc.geometry['lowerKulfanCoefficients'] = [c.magnitude for c in afl.lowerCoefficients]

                    rc.mutableLock  = True

                    results.append(rc)

        results_lists = comm.gather(results,root=0)
        results = []
        if rank == 0:
            for rse in results_lists:
                results += rse


            saveAnalysisCase(airfoilDataFolder+os.sep+airfoilName+'_data_%s.json'%(ii),results)







            # rcDicts.append(rcDict)

        # pool = Pool()  
        # resDicts = pool.map(generateXfoilResult, rcDicts)
        # pool.close()

        # for i,rc in enumerate(rcs):
        #     resDict = resDicts[i]
        #     if resDict is not None:
        #         if xfoil_mode == 'alfa':
        #             rc.outputs['cl']      = resDict['cl']
        #         else:
        #             rc.outputs['alpha']   = resDict['alpha']

        #         rc.outputs['bl_data'] = resDict['bl_data']
        #         rc.outputs['cd']      = resDict['cd']
        #         rc.outputs['cm']      = resDict['cm']
        #         rc.outputs['cp_data'] = resDict['cp_data']
        #         rc.outputs['xtr_l']   = resDict['xtr_l']
        #         rc.outputs['xtr_u']   = resDict['xtr_u']

        #         rc.geometry = {}
        #         rc.geometry['upperKulfanCoefficients'] = [c.magnitude for c in afl.upperCoefficients]
        #         rc.geometry['lowerKulfanCoefficients'] = [c.magnitude for c in afl.lowerCoefficients]

        #         rc.mutableLock  = True

        #         results.append(rc)

        # saveAnalysisCase(airfoilDataFolder+os.sep+airfoilName+'_data_%s.json'%(ii),results)

if rank == 0 :
    logfile = open(logfileName,'a')
    logfile.write('done\n')
    logfile.close()
    # print('done')