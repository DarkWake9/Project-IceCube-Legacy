import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize


def hvovec(lon1, lat1, lon2, lat2, rad=False):

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


def angdecfinder(psrno, declcut):

    '''

    Parameters
    ----------
    psrno : int
        Pulsar number from the ATNF catalogue.
    
    declcut : float
        Declination cut-off in degrees.



    Returns
    -------
    float
    Returns the absolute difference between the declinations of pulsar {psrno}\n
    and each of the neutrinos within the specified declcut
    '''
    #aang = hvovec(msra[psrno], msdec[psrno], icra, icdec)
    decdf = list(np.abs(np.subtract(icdec, msdec[psrno])))
    nuind = []
    for i in range(0, len(decdf)):
        if decdf[i] < declcut:
            nuind.append(i)
            #decdf.pop(i)
            #aang[i] = -1
    fdecdf = list(np.abs(np.subtract(icdec[nuind], msdec[psrno])))
    return [nuind, decdf]

def S_ij(psrno, nuind): 

    f'''
    Returns the signal PDF for the {psrno}th pulsar 
    and the all neutrinos within the specified declcut as given in {nuind}
    '''

    ang2 = hvovec(msra[psrno], msdec[psrno], icra[nuind], icdec[nuind], rad=True) ** 2
    sg = np.deg2rad(icang[nuind]) ** 2
    return np.exp(-1 * ang2 / (2 * sg)) / (2 * np.pi * sg)

def bgs(psrno, cone, twopi = False):

    '''
    Parameters
    ----------
    psrno : int
        Pulsar number from the ATNF catalogue.
    
    cone : float
        Solid angle cone in degrees.

    twopi : Boolean, optional
        (default) False -> The background PDF is NOT divided by 2pi
        True -> The background PDF is divided by 2pi
    
    Returns
    -------
    Returns the background PDF for a given pulsar {psrno} within\n 
    a solid angle cone of {cone} degrees for SINGLE STACKED ANALYSIS
    '''

    if twopi == True:
        s_ang = (np.sin(np.deg2rad(msdec[psrno] + cone)) - np.sin(np.deg2rad(msdec[psrno] - cone)))*(2 * np.pi)
        return 1/s_ang
    else:
        s_ang = (np.sin(np.deg2rad(msdec[psrno] + cone)) - np.sin(np.deg2rad(msdec[psrno] - cone)))
        return 1/s_ang

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


def TS(i, x, S, B, Ns):

    '''
    Returns the Test Stastic value fot {i}^th grb
    at $n_s$ = {x} for its parameters S, B, Ns
    '''
    
    return 2*np.sum(np.log(Pr(x,  Ns, S, B)/B))


def ns_for_TSmax(i, S, B, Ns):

    '''
    Parameters
    ----------
    i : int
        Pulsar number from the ATNF catalogue.

    Ns : int
        No.of neutrinos used for analysis

    S : float ->        Signal PDF

    B : float ->       Background PDF

    Returns
    -------
    Returns the value of $n_s$ for which
    the TS is maximum for {i}^th grb
    '''

    #returns the TSmax for i^th GRB
    nll = lambda x: -TS(i,x, S, B, Ns)
    soln = minimize(nll, 3 ,bounds=((0,None),))
    ns = float(soln.x)
    return ns

def plotpsr(psrno = 0):
    
    '''
    Parameters
    ----------
    psrno : int
        Pulsar number from the ATNF catalogue.

    Returns:    None\n
    Plots the signal PDF for the {psrno}th pulsar and\n
    prints the maxns and TSmax for the same

    '''



    cut = 5
    cone = 5
    nuind, decdf = angdecfinder(psrno, cut)
    S = S_ij(psrno, nuind)
    B = bgs(psrno, cone)
    Ns = len(nuind)
    tns=[TS(psrno, i, S, B, Ns) for i in range(0, 1000)]
    print(Ns)
    plt.plot(range(0, 1000), tns)
    plt.show()
    print([ns_for_TSmax(psrno, S, B, Ns), TS(psrno, ns_for_TSmax(psrno, S, B, Ns), S, B, Ns)])