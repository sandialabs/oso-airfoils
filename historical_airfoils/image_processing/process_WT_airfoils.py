import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as spo
import math
from kulfan import Kulfan
import os
import matplotlib.image as mpimg

def Kulfan_residual(coeff, *args):
        zeta_T = coeff[-1]
        coeff = coeff[0:-1]
        psi,zeta = args
        n = coeff.size - 1
        N1 = 0.5
        N2 = 1.0
        
        C = psi**(N1)*(1-psi)**(N2)
        S = np.zeros([n+1, psi.size])
        S_temp = np.zeros([n+1, psi.size])
        zeta_temp = np.zeros([n+1, psi.size])
        for i in range(0,n+1):
            S[i,:]=math.factorial(n)/(math.factorial(i)*math.factorial(n-i))  *  psi**i  *  (1-psi)**(n-i)
            S_temp[i,:] = coeff[i]*S[i,:]
            zeta_temp[i,:] = C*S_temp[i,:]+psi*zeta_T
        Sf = np.sum(S_temp, axis=0)
        zeta_coeff = psi**(N1)*(1-psi)**(N2)*Sf+psi*zeta_T
        return zeta-zeta_coeff

wd = "png_files"

combined = os.listdir(wd)
files   = list(sorted([f for f in combined if os.path.isfile(wd + os.sep + f)]))
folders = list(sorted([f for f in combined if not os.path.isfile(wd + os.sep + f)]))
filteredFiles = []
for fl in files:
    if fl[0]!='.' and fl[0]!='_' and ('.png' in fl):
        filteredFiles.append(fl)
        

for fl in filteredFiles:
    fig = plt.figure(figsize=(20,12))
    img = mpimg.imread(wd + os.sep + fl)
    imgplot = plt.imshow(img, aspect='equal', alpha=0.2)

    Nrows = len(img)
    Ncols = len(img[0])

    zeroRow    = 0
    zeroColumn = 0
    maxColumn  = 0
    minColumn  = len(img[0])

    min_y_arr = []
    max_y_arr = []

    for j in range(0,len(img[0])):
        min_y = len(img)
        max_y = 0
        for i in range(0,len(img)):
            quad = img[i][j]
            if quad[0] <= 0.1 and quad[1]>=.9 and quad[2]<=0.1:
                zeroRow    = i
                zeroColumn = j
            if quad[0]!=1 or quad[1]!=1 or quad[2]!=1:
                minColumn = min([minColumn,j])
                maxColumn = max([maxColumn,j])

                min_y = min([min_y, i])
                max_y = max([max_y, i])

        min_y_arr.append(min_y)
        max_y_arr.append(max_y)

    plt.plot([zeroColumn, maxColumn],[zeroRow, zeroRow],'b')
    plt.plot(min_y_arr,'g')
    plt.plot(max_y_arr,'r')

    upperCoordinates = []
    lowerCoordinates = []

    sorter = 'u'

    for j in range(minColumn, maxColumn):
        for i in range(0,len(img)):
            quad = img[i][j]

            x =      (j - zeroColumn)/(maxColumn-zeroColumn)
            y = -1 * (i - zeroRow)   /(maxColumn-zeroColumn)
            
            zeroCut = 0.15
            if fl == 'du_93-w-210.png':
                zeroCut = 0.8
            if x<= zeroCut:
                if y>=0:
                    sorter = 'u'
                else:
                    sorter = 'l'
            else:
                if i <= (min_y_arr[j]+max_y_arr[j])/2.0 :
                    sorter = 'u'
                else:
                    sorter = 'l'

            cutoff = 0.8
            if quad[0] <= cutoff and quad[1] <= cutoff and quad[2] <= cutoff and x>=0:
                if sorter=='u':
                    upperCoordinates.append([x,y])
                else:
                    lowerCoordinates.append([x,y])

    lc = np.array(lowerCoordinates)
    uc = np.array(upperCoordinates)

    plt.plot(lc[:,0]*(maxColumn-zeroColumn)+zeroColumn,-1*lc[:,1]*(maxColumn-zeroColumn)+zeroRow,'c.')
    plt.plot(uc[:,0]*(maxColumn-zeroColumn)+zeroColumn,-1*uc[:,1]*(maxColumn-zeroColumn)+zeroRow,'m.')

    psiu = uc[:,0]
    zetau = uc[:,1]

    psil = lc[:,0]
    zetal = lc[:,1]

    fit_order = 8
    coeff_guess = np.ones(fit_order+1)
    resu=spo.leastsq(Kulfan_residual, coeff_guess.tolist(), args=(psiu, zetau))

    coeff_guess = -1*np.ones(fit_order+1)
    resl=spo.leastsq(Kulfan_residual, coeff_guess.tolist(), args=(psil, zetal))

    afl = Kulfan(TE_gap = resu[0][-1] - resl[0][-1])
    afl.upperCoefficients = np.asarray(resu[0][0:-1]).tolist()
    afl.lowerCoefficients = np.asarray(resl[0][0:-1]).tolist()
    plt.plot(afl.xcoordinates, afl.ycoordinates)
    
    plt.plot(afl.xcoordinates*(maxColumn-zeroColumn)+zeroColumn,-1*afl.ycoordinates*(maxColumn-zeroColumn)+zeroRow,'k')
    plt.savefig('verification_plots/' + fl[0:-4] + '.png')
    afl.write2file('fitted_datfiles/' + fl[0:-4] + '.dat')
    plt.close()
    
    
combined = os.listdir(wd)
files   = list(sorted([f for f in combined if os.path.isfile(wd + os.sep + f)]))
folders = list(sorted([f for f in combined if not os.path.isfile(wd + os.sep + f)]))
datFiles = []
for fl in files:
    if ('.dat' in fl):
        datFiles.append(fl)
        
        
correctedTauDict = {
    'afl-du-family1.dat' :  0.18 ,
    'afl-du-family2.dat' :  0.21 ,
    'afl-du-family3.dat' :  0.25 ,
    'afl-du-family4.dat' :  0.30 ,
    'afl-riso-b1.dat' :     0.15 ,
    'afl-riso-b2.dat' :     0.18 ,
    'afl-riso-b3.dat' :     0.21 ,
    'afl-riso-b4.dat' :     0.24 ,
    'afl-riso-b5.dat' :     0.30 ,
    'afl-riso-b6.dat' :     0.36 ,
    'afl-riso-p1.dat' :     0.15 ,
    'afl-riso-p2.dat' :     0.18 ,
    'afl-riso-p3.dat' :     0.21 ,
    'afl-riso-p4.dat' :     0.24 ,
    'du_91-w2-250.dat' :    0.25 ,
    'du_93-w-210.dat' :     0.21 ,
    'du_96-w-180.dat' :     0.18,
    'du_97-w-300.dat' :     0.30 ,
    'riso-a-12.dat' :       0.12 ,
    'riso-a-15.dat' :       0.15 ,
    'riso-a-18.dat' :       0.18 ,
    'riso-a-21.dat' :       0.21 ,
    'riso-a-24.dat' :       0.24 ,
    'riso-a-27.dat' :       0.27 ,
    'riso-a-30.dat' :       0.30 ,
}

du_corrections = {
    'afl-du-family1.dat' : 'du_96-w-180.dat',
    'afl-du-family2.dat' : 'du_93-w-210.dat',
    'afl-du-family3.dat' : 'du_91-w2-250.dat',
    'afl-du-family4.dat' : 'du_97-w-300.dat',
}
        
        
for datfile in datFiles:
    afl = Kulfan()
    afl.readfile(wd + os.sep + datfile)
    afl.scaleThickness(correctedTauDict[datfile])
    if datfile not in list(du_corrections.values()):
        if 'afl-du-family' in datfile:
            if datfile in list(du_corrections.keys()):
                afl.write2file("corrected_datfiles/"+du_corrections[datfile])
        elif 'riso-b' in datfile:
            afl.write2file("corrected_datfiles/riso-b-%d.dat"%(int(100*afl.tau)))
        elif 'riso-p' in datfile:
            afl.write2file("corrected_datfiles/riso-p-%d.dat"%(int(100*afl.tau)))
        else:
            afl.write2file("corrected_datfiles/"+datfile)
