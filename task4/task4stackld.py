import numpy as np
import pandas as pd
import os
import multiprocessing as mul
from multiprocessing import Process
import matplotlib.pyplot as plt
#import seaborn as sns
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from IPython.display import clear_output
from scipy.stats import chi2
from scipy.stats import norm
from scipy.stats import gaussian_kde
import scipy.stats as sct
from astropy.coordinates import SkyCoord as scr
from astropy import units as u
from time import sleep
from numba import jit





#@jit(nopython=True)
def psr_wt_quick(nusample_wall, psrno, gamma = 1):

    '''
    Parameters
    ----------
    nusample_wall : float
        The index of catalogue of the given neutrino
    psrno : int
        The index of the pulsar in the ATNF catalogue
    gamma : float, optional
        (default) 1
        The spectral index of the power law
    Returns the weight of psrno^th pulsar for a given neutrino sample {nusample_wall}
    '''
    
    #nu_e = icdata['log10(E/GeV)'].values[nu]
    psr_decl = icdec[psrno]
    #upt = upt_icparts[nusample_wall]
    t_upt = uptdata[nusample_wall]['MJD_stop[days]'].values[-1] - uptdata[nusample_wall]['MJD_start[days]'].values[0]
    t_upt*=86400
    d_ind = 0
    for i in range(len(dec_nu) - 1):
        if dec_nu[i] <= psr_decl and dec_nu[i+1] >= psr_decl:
            d_ind = i
            break
    ea_temp = eadata[nusample_wall][d_ind*40:(d_ind+1)*40]
    earea = np.array(ea_temp['A_Eff[cm^2]'].values) * 1e-4
    weight_kj = t_upt * np.sum(np.multiply(np.multiply(earea, np.power(e_nu, gamma)), de_nu))
    sleep(1e-8)
    return weight_kj                           
            


@jit(nopython=True)
def hvovec(lon1=0.0, lat1=0.0, lon2=0.0, lat2=0.0, rad=False):

    '''
    Parameters
    ----------
    lon1 : float
        Longitude of first point.
    lat1 : float
        Latitude of first point.


    lon2 : float
        Longitude of second point.
    lat2 : float
        Latitude of second point.

    rad : Boolean, optional
        (default) False -> Angle is returned in Degrees
        True -> Angle is returned in radians
    

    Returns
    -------
    float
        Returns the haversine angle between two vectors given their latitude and longitude
    

    Warning
    -------
        This function assumes the input to be in degrees and of equal length\n
        or one of the input pairs to be a single value and the other pair to be an array
    '''

    #Convert decimal degrees to Radians:
    lon1 = np.deg2rad(lon1)
    lat1 = np.deg2rad(lat1)
    lon2 = np.deg2rad(lon2)
    lat2 = np.deg2rad(lat2)

    #Implementing Haversine Formula: 
    dlon = np.subtract(lon2, lon1)
    #dlat = np.subtract(lat2, lat1)

    a = np.add(np.multiply(np.sin(lat1), np.sin(lat2)), np.multiply(np.multiply(np.cos(lat1), np.cos(lat2)), np.cos(dlon)))

    if rad == True:
        return np.arccos(a)
    else:
        return np.abs(np.rad2deg(np.arccos(a)))




@jit(nopython=True)
def S_ij(nu): 

    '''
    Returns the signal PDF for the {psrno}th pulsar and nuind_inp neutrino
    '''
    ang2 = hvovec(msra, msdec, icra[nu], icdec[nu], rad=True) ** 2#, icra[nuind], icdec[nuind], rad=True) ** 2
    sg = np.deg2rad(icang[nu]) ** 2
    return np.divide(np.exp(-1 * np.divide(ang2, 2*sg)), (2 * np.pi * sg))

@jit(nopython=True)
def S_i(nu, all_weights, sum_weights, wall):

    '''
        Parameters
        ----------
        nu : int
            Index of the neutrino in the sample

        Returns
        -------
        Returns the signal PDF for the {nu}th neutrino with all pulsars 
        i.e S_i = wieghted-sum_j S_ij
    '''



    sij = S_ij(nu)
    
    #sleep(1e-8)
    return np.sum(np.multiply(sij, all_weights[wall] / sum_weights[wall]))
    #return np.dot(sij, all_weights[wall] ) / sum_weights[wall]






@jit(nopython=True)
def Bi_stacked(nu, cone):

    '''
    Parameters
    ----------
    nu : int
        Index of the neutrino from IceCube sample
    cone : float
        Cone angle in degrees.
    

    Returns
    -------
    float
        Returns the background PDF for the {nu}th neutrino
    '''

    #count = 0
    #count = np.count_nonzero(np.abs(np.subtract(icdec, icdec[nu])) <= cone)
    #print(count)
    count = 0
    for i in range(len(icdec)):
        if abs(icdec[i] - icdec[nu]) <= cone:
            count+=1
    #print(count)
    binwidth = (np.sin(np.deg2rad(icdec[nu] + cone)) - np.sin(np.deg2rad(icdec[nu] - cone)))*2*np.pi
    #sleep(1e-8)
    return count/(binwidth * lnu)
    

@jit(nopython=True)
def Pr(x, Ns, S, B):

    '''
    Parameters
    ----------
    x : int
        Assumed no.of associated events

    Ns : int
        No.of neutrinos used for analysis

    S : float
        Signal PDF
    B : float
        Background PDF
    
    Returns
    -------
    float 
        Returns the probability of the selected set of neutrinos being associated\n
        with a given pulsar with {Ns} neutrinos, {S} signal and {B} background PDF and {x} assumed associated events
    '''

    nsN = x/Ns
    return np.add(np.multiply(nsN , S), np.multiply(np.subtract(1, nsN), B))





@jit(nopython=True)
def wall_nu(nu):
    
        '''
        Parameters
        ----------
        nu : int
            Index of the neutrino from IceCube sample
        
    
        Returns
        -------
        int
            Returns the index of the wall in which the {nu}th neutrino lies
        '''
        wall = 0
        for i in range(len(vec_uptparts)-1):
            if vec_uptparts[i] <= nu and vec_uptparts[i+1] > nu:
                wall = i
                break        
        return wall




def weights(gamma):
    all_weights= []
    for i in range(len(upt_icparts)-1):
    #print(upt_icparts[i])
        all_weights.append([psr_wt_quick(i, psrno, gamma) for psrno in range(p)])
    return np.asfarray(all_weights)

#nuind, ang = angfinder(psrno, cut)



@jit(nopython=True)
def sigbag_nu(nu):
    '''
        Returns the signal and background PDF for the {nu}th neutrino
    '''

    #sleep(1e-8)

    if nu % 100000 == 0:
        print(nu)
    wall = wall_nu(nu)
    Ns = vec_uptparts[ + 1] - vec_uptparts[wall_nu(nu)]
    return [S_i(nu, all_weights, sum_weights, wall), Bi_stacked(nu, cone),Ns]



@jit(nopython=True)
def TS_st(x, S, B, Ns):

    '''
    Returns
    ----------
        Returns the Test Stastic value at
        $n_s$ = {x} for its parameters S, B, Ns
    '''
    if x >=0:
        return np.sum(np.asfarray(2*np.log(Pr(x,  Ns, S, B)/B)))
    else:
        return -1


def ns_for_TSmax_st(S, B, Ns):
    '''
    Returns the value of $n_s$ for which
    the TS is maximum for {i}^th grb
    '''

    #returns the TSmax for i^th GRB
    nll = lambda x: -TS_st(x, S, B, Ns)
    soln = minimize(nll, 3 , bounds = [(0, None)], tol=1e-12)
    ns = np.round(soln.x, 6)[0]
    #print(soln.success)
    return ns


@jit(nopython=True)
def Ts_arr(x, S_all, B_all, Ns):
    #Ts_arr = lambda x:
    sum = 0.0
    for i in range(lnu):
        sum += TS_st(x, S_all[i], B_all[i], Ns)
    return sum


def drive_stack(gamma):
    all_weights = weights(gamma)
    sum_weights = [np.sum(i) for i in all_weights]
    sum_weights = np.asfarray(sum_weights)
    cut = 5
    cone = 5
    
    pool = mul.Pool(8, maxtasksperchild=100)
    op_async = pool.map_async(sigbag_nu, range(len(icdec)))
    tmp = op_async.get()
    pool.close()
    pool = []
    op_async = []
    tmp = np.asfarray(tmp)
    S_all = tmp[:,0]
    B_all = tmp[:,1]
    Ns = lnu
    nll = lambda x: -Ts_arr(x, S_all, B_all, Ns)
    soln = minimize(nll, 3 , bounds = [(0, None)], tol=1e-12)
    ns = np.round(soln.x, 6)[0]
    allts=[Ts_arr(i,S_all, B_all, Ns) for i in range(0, 1000)]
    return [ns, allts]  


#gaussian = lambda x, mean, var:  amp*norm.pdf(x, mean, var)
def gaussian(x,mean,sd):
#    mean=0.114
    A=p/3.0
    return A*np.exp(-(x - mean) ** 2.0 / sd ** 2.0)
loggaussian = lambda x, mean, var:  np.log(gaussian(x, mean, var))
#gparam, gerr = curve_fit(gaussian, x2, b2,  p0=[0.1, 1.0],sigma=yerr2, absolute_sigma=True, maxfev = 10000, method='dogbox')
#gx = np.linspace(0, wsts[-1], 100)