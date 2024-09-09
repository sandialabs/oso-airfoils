import os
import multiprocessing
import datetime
import numpy as np
import matplotlib.pyplot as plt
# %matplotlib inline
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

# afl = generateAirfoil('du_97-w-300')
# print(afl.upperCoefficients)
# print(afl.lowerCoefficients)