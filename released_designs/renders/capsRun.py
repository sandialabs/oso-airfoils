import math
import pandas as pd
import numpy as np
from mpi4py import MPI


# from kulfan import Kulfan
from ada.geometry.airfoils.kulfan import Kulfan

st = pd.ExcelFile('IEA-22-280-RWT_tabular.xlsx')

# ========================================================
# Load blade geometry parameters
# ========================================================
blade_geometry = pd.read_excel(st, 'Blade Geometry')

airfoils = blade_geometry.iloc[0:16]
airfoil_spans = airfoils['Spanwise position [r/R]'].to_numpy()
airfoil_names = airfoils['Airfoil name'].to_list()

detailed_geometry = blade_geometry.iloc[19:]
normSpans = detailed_geometry['Spanwise position [r/R]'].to_numpy()
chords    = detailed_geometry['Airfoil name'].to_numpy()
twists    = detailed_geometry['Unnamed: 2'].to_numpy()*180.0/np.pi
pas       = detailed_geometry['Unnamed: 3'].to_numpy()
spans     = detailed_geometry['Unnamed: 4'].to_numpy()
prebends  = detailed_geometry['Unnamed: 5'].to_numpy()
sweeps    = detailed_geometry['Unnamed: 6'].to_numpy()

# ========================================================
# Load airfoils
# ========================================================
airfoil_data = {}
for st_nm in st.sheet_names:
    if 'Airfoil ' in st_nm:
        coordinates = []
        nm = st_nm.replace('Airfoil ','')
        # print(nm)
        airfoil_sheet = pd.read_excel(st, st_nm)
        coordinates_end = 0
        for i in range(4,1000):
            coordinates_end = i
            rw = airfoil_sheet.iloc[i].to_numpy()
            if math.isnan(rw[0]):
                break
            else:
                coordinates.append([rw[0], rw[1]])
                
        coordinates = np.array(coordinates)
        coordinates[:,0] -= min(coordinates[:,0])
        airfoil_data[nm] = {}
        airfoil_data[nm]['coordinates'] = coordinates
        # plt.plot(coordinates[:,0], coordinates[:,1])
        # plt.axis('equal')
        
        airfoil_data[nm]['data'] = []
        currentDataStack = []
        for i in range(coordinates_end, len(airfoil_sheet)):
            rw = airfoil_sheet.iloc[i].to_numpy().tolist()
            if i == len(airfoil_sheet)-1:
                currentDataStack.append(rw)
                rw = ['alpha',0,0]
            
            if isinstance(rw[0], (float, int)):
                if not math.isnan(rw[0]):
                    currentDataStack.append(rw)
            elif isinstance(rw[0],str):
                if 'alpha' in rw[0]:
                    if len(currentDataStack)!=0:
                        dta = np.array(currentDataStack)
                        packedData = {'alpha':dta[:,0]*180.0/np.pi, 'cl':dta[:,1], 'cd':dta[:,2], 'cm':dta[:,3]}
                        airfoil_data[nm]['data'].append(packedData)
                        currentDataStack = []
                else:
                    # do nothing
                    assert('alpha' not in rw[0])
                    assert(isinstance(rw[0],str))

            else:
                assert(math.isnan(rw[0]) or rw[0]=='Reynolds')
        
# ========================================================
# Make the Kulfan objects
# ========================================================
for nm, afl in airfoil_data.items():
    kfl = Kulfan()
    kfl.fit2coordinates(afl['coordinates'][:,0], afl['coordinates'][:,1], fit_order=10)
    airfoil_data[nm]['kulfan'] = kfl

# ========================================================
# Stack all the detailed airfoils
# ========================================================
esp_string = 'mark \n'

circle_string = ''
theta = 25 * np.pi/180.0
circle_string += 'skbeg %f %f 0 \n'%(np.cos(theta), np.sin(theta))
circle_string += 'cirarc 0  1 0 -1 0 0 \n'
circle_string += 'cirarc 0 -1 0 %f %f 0 \n'%(np.cos(theta), -1*np.sin(theta))
circle_string += 'cirarc 1 0 0 %f %f 0 \n'%(np.cos(theta), np.sin(theta))
circle_string += 'skend \n'
circle_string += 'scale 0.5 \n'
circle_string += 'translate 0.5 0 0 \n'

# FB80 has a larger TE gap than FB90, which leads to significant blending issues
# this is in a region where flow is slow
airfoil_data['FB90']['kulfan'].constants.TE_gap *= 1.3

for i in range(0,len(normSpans)):
    normSpan = normSpans[i]
    chord    = chords[i]
    twist    = twists[i]
    pa       = pas[i]
    span     = spans[i]
    prebend  = prebends[i]
    sweep    = sweeps[i]
    
    if normSpan <= 0.02:
        esp_string += circle_string

    elif normSpan <= airfoil_spans[2]:
        esp_string += airfoil_data['FB90']['kulfan'].toESP()
        
    else:

        if normSpan != 1:
            right_airfoil_index = np.argmax(airfoil_spans>normSpan)
            left_airfoil_index  = right_airfoil_index - 1
        else:
            right_airfoil_index = -1
            left_airfoil_index  = -1
            
        interval_length = airfoil_spans[right_airfoil_index] - airfoil_spans[left_airfoil_index]
        if interval_length == 0:
            progres_frac = 1
        else:    
            progress_frac = (normSpan - airfoil_spans[left_airfoil_index]) / interval_length
        
        xcoords_left  = airfoil_data[airfoil_names[left_airfoil_index]]['kulfan'].xcoordinates
        xcoords_right = airfoil_data[airfoil_names[right_airfoil_index]]['kulfan'].xcoordinates
        ycoords_left  = airfoil_data[airfoil_names[left_airfoil_index]]['kulfan'].ycoordinates
        ycoords_right = airfoil_data[airfoil_names[right_airfoil_index]]['kulfan'].ycoordinates
        xcoords = (1-progress_frac)*xcoords_left + progress_frac*xcoords_right
        ycoords = (1-progress_frac)*ycoords_left + progress_frac*ycoords_right
        
        kfl = Kulfan()
        kfl.fit2coordinates(xcoords, ycoords)
        esp_string += kfl.toESP()
        
    esp_string += 'translate %f 0 0 \n'%(-1*pa)
    esp_string += 'scale %f \n'%(chord)
    esp_string += 'rotatez %f 0 0 \n'%(twist)
    esp_string += 'translate 0 %f %f \n'%(-prebend,span)
    esp_string += '\n'
    
esp_string += 'rule 1 \n'
f = open('22mw_test_fullDetail.csm','w')
f.write(esp_string)
f.close()

# ========================================================
# Stack a smoother approximation
# ========================================================
esp_smooth = 'mark \n'
for i, airfoil_span in enumerate(airfoil_spans):
    airfoil_name = airfoil_names[i]
    
    spns = [airfoil_span]
    kfls = [ airfoil_data[airfoil_name]['kulfan'] ]
    
    if i in [10,11,12]:
        spns.append(0.5*(airfoil_span + airfoil_spans[i+1]))
        
        left_airfoil_index = i
        right_airfoil_index = i+1
        xcoords_left  = airfoil_data[airfoil_names[left_airfoil_index]]['kulfan'].xcoordinates
        xcoords_right = airfoil_data[airfoil_names[right_airfoil_index]]['kulfan'].xcoordinates
        ycoords_left  = airfoil_data[airfoil_names[left_airfoil_index]]['kulfan'].ycoordinates
        ycoords_right = airfoil_data[airfoil_names[right_airfoil_index]]['kulfan'].ycoordinates
        progress_frac = 0.5
        xcoords = (1-progress_frac)*xcoords_left + progress_frac*xcoords_right
        ycoords = (1-progress_frac)*ycoords_left + progress_frac*ycoords_right
        kfl = Kulfan()
        kfl.fit2coordinates(xcoords, ycoords)
        kfls.append(kfl)
        
    elif i == 13:
        pgfs = [0.25,0.5,0.75]
        left_airfoil_index = i
        right_airfoil_index = i+1
        xcoords_left  = airfoil_data[airfoil_names[left_airfoil_index]]['kulfan'].xcoordinates
        xcoords_right = airfoil_data[airfoil_names[right_airfoil_index]]['kulfan'].xcoordinates
        ycoords_left  = airfoil_data[airfoil_names[left_airfoil_index]]['kulfan'].ycoordinates
        ycoords_right = airfoil_data[airfoil_names[right_airfoil_index]]['kulfan'].ycoordinates
        
        for pgf in pgfs:
            spns.append(airfoil_span + pgf*(airfoil_spans[i+1]-airfoil_span))
            progress_frac = pgf
            xcoords = (1-progress_frac)*xcoords_left + progress_frac*xcoords_right
            ycoords = (1-progress_frac)*ycoords_left + progress_frac*ycoords_right
            kfl = Kulfan()
            kfl.fit2coordinates(xcoords, ycoords)
            kfls.append(kfl)
        
    else:
        pass
    
    for j, spn in enumerate(spns):
        kfl = kfls[j]
        if airfoil_name == 'circular':
            esp_smooth += circle_string
        else:
            esp_smooth += kfl.toESP()
            
        pa_smooth      = np.interp(np.float64(spn),np.array(normSpans, dtype=np.float64),np.array(pas, dtype=np.float64))
        chord_smooth   = np.interp(np.float64(spn),np.array(normSpans, dtype=np.float64),np.array(chords, dtype=np.float64))
        twist_smooth   = np.interp(np.float64(spn),np.array(normSpans, dtype=np.float64),np.array(twists, dtype=np.float64))
        prebend_smooth = np.interp(np.float64(spn),np.array(normSpans, dtype=np.float64),np.array(prebends, dtype=np.float64))
            
        esp_smooth += 'translate %f 0 0 \n'%(-1*pa_smooth)
        esp_smooth += 'scale %f \n'%(chord_smooth)
        esp_smooth += 'rotatez %f 0 0 \n'%(twist_smooth)
        esp_smooth += 'translate 0 %f %f \n'%(-prebend_smooth,spn*max(spans))
        esp_smooth += '\n'
        
esp_smooth += 'blend 1 \n'
esp_smooth += 'dimension tess_params 1 3 0 \n'
esp_smooth += 'set tess_params "1.0; 0.1; 10;" \n'
esp_smooth += 'store blade \n'
esp_smooth += 'restore blade \n'
esp_smooth += 'attribute .tParams tess_params \n'


f = open('22mw_test_smooth.csm','w')
f.write(esp_smooth)
f.close()

# # ========================================================
# # Prep for xfoil
# # ========================================================
# runList = []

# for afl in airfoil_names:
#     if 'FB' not in afl and 'SNL' not in afl and afl != 'circular':
#         runList.append(afl)
        

# # ========================================================
# # Run xfoil
# # ========================================================
# from runSingleCase import runSingleCase
# import time
# import json

# runList = [runList[0]]

# alphaList = np.linspace(-30,30,31)
# # reList    = 10**np.linspace(6,8,7)
# # machList  = np.linspace(0,0.6,7)
# # truList   = np.linspace(0,1.0,11)
# # trlList   = np.linspace(0,1.0,11)
# # nList     = np.array([3,5,7,9,11])

# reList    = [1e7]
# machList  = [0.0]
# truList   = [1.0]
# trlList   = [1.0]
# nList     = [9.0]

# ctr = 0
# case_count = len(alphaList) * len(reList) * len(machList) * len(truList) * len(trlList) * len(nList)
# start = time.time()
# for ky in runList:
#     for alpha in alphaList:
#         for re in reList:
#             for M in machList:
#                 for tru in truList:
#                     for trl in trlList:
#                         for n in nList:
#                             ctr += 1
#                             progress = ctr/case_count
#                             if ctr%100 == 0:
#                                 print('progress: %.3f'%(100*progress), 'time: %.2f'%(time.time()-start))
                                
    
#                             afl = airfoil_data[ky]['kulfan']

#                             try:
#                                 dta = runSingleCase('alpha', afl, alpha, 
#                                                     Re = re, M = M,
#                                                     xtr_u=tru, xtr_l=trl, N_crit=n,
#                                                     max_iter=200
#                                                 )
#                             except:
#                                 pass
                            
#                             with open('data/case_%d.json'%(ctr), 'w') as fp:
#                                 json_string = json.dumps(dta, indent=4)
#                                 json_string = json_string.replace("NaN","null")
#                                 fp.write(json_string)


# # print(case_count)