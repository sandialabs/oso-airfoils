#  These geometry functions serve as the defaults for airfoil geometry parameters
#  They have been roughly fit to existing wind turbine airfoils, but are not meant to be precise
#  They have been modified to that the trends should hold valid for a reasonable domain of airfoil thicknesses
#  For precise design, these values should be provided by the user in the JSON input file

import numpy as np
from scipy.special import logsumexp

# Data used to for the OSO family, which can be compared to
# taus    = [      0.15,       0.18,       0.21,       0.24,       0.27,       0.30,       0.33,       0.36 ]
# te_gaps = [   0.00196,    0.00230,    0.00262,    0.00751,    0.01012,    0.01140,    0.01140,    0.01140 ]
# Ixx_con = [0.00011000, 0.00017438, 0.00027518, 0.00041096, 0.00058321, 0.00079640, 0.00105795, 0.00137822 ]
# Iyy_con = [0.00397999, 0.00436351, 0.00493714, 0.00561409, 0.00633417, 0.00706380, 0.00779600, 0.00855043 ]
# Izz_con = [0.00408809, 0.00454606, 0.00521632, 0.00602287, 0.00691323, 0.00785849, 0.00885328, 0.00991577 ]
# A_con   = [0.08700496, 0.09995900, 0.11477620, 0.13051205, 0.14660942, 0.16289864, 0.17959744, 0.19731100 ]
# ler_con = [     0.007,      0.008,       0.01,      0.025,       0.03,       0.04,       0.06,       0.08 ]

def TE_gap_function(tau):
    L_lo = 0.0049
    k_lo = 10.8885
    tau0_lo = 0.0447
    b_lo = -0.0082
    L_hi = 0.0085
    k_hi = 89.6093
    tau0_hi = 0.2373
    b_hi = 0.0063

    low  = L_lo / (1.0 + np.exp(-k_lo * (tau - tau0_lo))) + b_lo
    high = L_hi / (1.0 + np.exp(-k_hi * (tau - tau0_hi))) + b_hi
    return low + high

def cone_angle_function(tau):
    # if tau <= 0.18:
    #     return 10.0
    # if tau <= 0.27:
    #     return 5.0
    # return 0.0

    k1 = 2368.3255
    tau01 = 0.1800
    k2 = 2322.2349
    tau02 = 0.2700

    drop1 = 5.0 / (1.0 + np.exp(-k1 * (tau - tau01)))
    drop2 = 5.0 / (1.0 + np.exp(-k2 * (tau - tau02)))
    return 10.0 - drop1 - drop2

def Ixx_function(tau):
    # Ixx = 0.15194621  * tau**4 - 0.12744749  * tau**3 + 0.05937931 * tau**2 -0.00976398 * tau + 0.00059178
    a = 0.0288
    b = 2.9791
    return a * tau ** b
    
def Iyy_function(tau):
    # Iyy = 1.36882085  * tau**4 - 1.60347502  * tau**3 + 0.70416555 * tau**2 -0.11306505 * tau + 0.00981479
    a = 0.0237
    return a * tau 
    
def Izz_function(tau):
    # Izz = 1.21138268  * tau**4 - 1.41411698  * tau**3 + 0.64498366 * tau**2 -0.10370946 * tau + 0.00929176
    a = 0.0263
    return a * tau
        
def area_function(tau):
    # area = 19.92943549 * tau**4 - 21.37544951 * tau**3 + 8.68562194 * tau**2 -1.04486644 * tau + 0.1103613
    a = 0.5470
    b = 2.9791
    return a * tau

def ler_function(tau):
    a = 2.3920
    b = 3.3338
    return a * tau ** b

def min_radius_location_upper_function(tau):
    alpha = 10000.0
    a = alpha * np.array([0, (0.0045 / 0.03)])
    b = alpha * np.array([0.011, 0.011 - (0.0045 / 0.03)*0.271])
    return 1/alpha * logsumexp(a*tau + b)

def min_radius_location_lower_function(tau):
    alpha = 10000.0
    a = alpha * np.array([0, (0.035 / 0.03)])
    b = alpha * np.array([0.015, 0.015 - (0.035 / 0.03)*0.30])
    return 1/alpha * logsumexp(a*tau + b)