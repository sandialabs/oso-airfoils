from kulfan import Kulfan
import os
import numpy as np
import json

import matplotlib
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})
import matplotlib
# colors = ['#0065cc', '#e69f00', '#009e73', '#d55e00', '#56b4ff', '#fca7c7', '#ede13f', '#666666', '#000000']
colors = ['#0065cc', '#e69f00', '#009e73', '#d55e00', '#56b4ff', '#fca7c7', '#ede13f', '#000000']
matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=colors)

import natsort

# combined = os.listdir('./oso-airfoils/postprocessing')
combined = os.listdir('viable_airfoils/')
folder = 'viable_airfoils/'
files     = files   = natsort.natsorted([f for f in os.listdir(folder) if '.txt' in f], alg=natsort.ns.IGNORECASE)
# folders = list(sorted([f for f in combined if not os.path.isfile(f) and f[0]!='.' and f[0]!='_']))
# print(files)
print(len(files))
plt.figure(figsize=(12, 8))
allstr = ''
for f in files:
    f = open(folder + f, 'r')
    block = f.read()
    lines = block.split('\n')
    f.close()
    pstr = '['
    ents = lines[0].split(',')
    for ent in ents:
        pstr += ent.rjust(9) + ','
    pstr = pstr[0:-1] + '],'
    allstr += pstr + '\n'
    afl = Kulfan(TE_gap=0.01)
    coeffs = np.array([float(x) for x in lines[0].split(',')])
    afl.upperCoefficients = coeffs[:len(coeffs)//2]
    afl.lowerCoefficients = coeffs[len(coeffs)//2:]
    plt.plot(afl.xcoordinates, afl.ycoordinates, color = colors[0], alpha = 0.01)

plt.grid(True)
plt.axis('equal')
plt.savefig('viable_airfoils.png',dpi=200)

f = open('all_viable_airfoils.txt', 'w')
f.write(allstr)
f.close()
print(allstr)
plt.axis('equal')

plt.figure(figsize=(12, 8))
allstr = ''
for f in files:
    f = open(folder + f, 'r')
    block = f.read()
    lines = block.split('\n')
    f.close()
    afl = Kulfan(TE_gap=0.01)
    coeffs = np.array([float(x) for x in lines[0].split(',')])
    for i,c in enumerate(coeffs):
        if i < len(coeffs)//2:
            plt.plot(i,c, '.', alpha = 0.01, color = colors[0])
        else:
            plt.plot(i,c, '.', alpha = 0.01, color = colors[1])
plt.grid(True)
plt.savefig('coefficient_distribution.png',dpi=200)