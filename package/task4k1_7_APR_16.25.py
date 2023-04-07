####
from core import download_IC
from core import download_ATNF
from core import readfiles
import numpy as np
import pandas as pd
import os
import multiprocessing as mul
from multiprocessing import Process
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import scipy.stats as sct
from astropy.coordinates import SkyCoord as scr
from astropy import units as u
from time import sleep
from numba import jit, njit, prange
from tqdm import tqdm
from core.signal_bag import *
from core.stacking_analysis import *
from core.req_arrays import *
import bisect


# UNCOMMENT FOR LINEAR BINS
# all_enu = np.linspace(10**11.001, 10**18.999, 1000)

# UNCOMMENT FOR LOGARITHMIC BINS of WIDTH 0.2 from logE/GeV = 2 to 10 AS GIVEN IN THE ICECUBE DATA
all_enu = e_nu_wall

# UNCOMMENT FOR LOGARITHMIC BINS
# all_enu = np.logspace(11.001, 18.999, 1000)
gamma_arr = [-1, -2.2, -2.5, -3]
phio = np.logspace(-38, -26, 1000) #CHANGING TO LINEAR BINS RESULTS IN STRAIGHT LINES

####
weight_obj =  [weights.weights(gamma).all_weights for gamma in gamma_arr]
sum_weights = [weights.weights(gamma).sum_weights for gamma in gamma_arr]
#Calculating weight/sum(weights) to avoid wrong ns = 0 entries 
normalized_wt = []
for gamma in range(4):
    temp = []
    for season in range(len(weight_obj[gamma])):
        temp.append(weight_obj[gamma][ea_season(season)]/sum_weights[gamma][ea_season(season)])
    normalized_wt.append(temp)

####
#generates ns^ for a single pulsar a single season and a single energy

@njit(fastmath=True)
def ns_singleseason_sing_psr_HAT(dec,enu, gamma, phi0 = 1, season=0):
    '''
    This function returns the number of signal events for a single pulsar as in EQN3 of 2205.15963
    -------------------

    Parameters
    ----------
    dec : float
        The declination of the pulsar in radians
    
    enu : float
        The neutrino energy in eV

    gamma : float
        The spectral index of the neutrino flux

    phi0 : float (optional)
        The normalization constant of the neutrino flux

    season : int (optional)
        The IceCube season number
    
    Returns
    -------
    float
        The number of signal events for the given parameters
    '''


    tt_upt = t_upt[season]
    if enu <= 1e11 or enu >= 1e19:
        return 0.0
    else:
        k=0
        l=0
        for i in range(0, len(e_nu_wall)):
            if e_nu_wall[i] <= enu and e_nu_wall[i+1] > enu:
                for j in range(0, len(dec_nu)):
                    if dec_nu[j] <= dec and dec_nu[j+1] > dec:
                        k=i
                        l=j
                        break
                break

        temp_ea = np.asarray(earea[ea_season(season)])[l*40 + k]
        return tt_upt * temp_ea * dfde(enu, gamma, phi0)     #in s cm2 eV


####
phi0t = 1e-30
season = 0
dec = msdec[0]
print("Decl = ", dec)
tmp = [ns_singleseason_sing_psr_HAT(dec, all_enu[i], gamma_arr[gamma], phi0t ,season=season) for i in range(len(all_enu))]

#UNEQUAL BINWIDTHS IN LINEAR SPACE BUT EQUAL IN LOGARITHMIC SPACE
print("ns^ = ",np.trapz(y=tmp, dx=np.diff(all_enu)))

#EQUAL BINWIDTHS IN BOTH LINEAR 10^9.2 AND LOGARITHMIC SPACE of 0.2
print("ns^ = ",np.trapz(y=tmp, dx = 1e9* (10**0.2)))

####



