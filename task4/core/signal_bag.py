from numba import jit, njit
import numpy as np
from scipy.optimize import minimize
from core import readfiles, weights
import multiprocessing as mul

all_data = readfiles.Data()
icdata = all_data.icdata
uptdata = all_data.uptdata
eadata = all_data.eadata
mspdata = all_data.mspdata
all_data = []

icwidths = [int(i) for i in "0 36900 107011 93133 136244 112858 122541 127045 129311 123657 145750".split(' ')]
ictimes = [float(i) for i in icdata['MJD[days]']]
icparts = [np.sum(icwidths[:i]) for i in range(1,len(icwidths)+1)]  #paritions of icdata for each season (IC40, IC59, IC79, IC86_I, IC86_II)

upt_icparts = icparts[:5]
upt_icparts.append(icparts[-1])


log_e = np.round(np.arange(2, 10.2, 0.2), 2) #log10(E/GeV) values range as in all 'effectiveArea' files

#dec_nu = mid point of Dec_nu_min and Dec_nu_max as in all 'effectiveArea' files
dec_nu = list(set(eadata[0]['Dec_nu_min[deg]'].values).union(set(eadata[0]['Dec_nu_max[deg]'].values)))

dec_nu.sort()
dec_nu = np.array(dec_nu)

e_nu = ((10**(log_e[:-1])+ 10**(log_e[1:]))/2)*1e9
de_nu = 1e9*(10**log_e[1:] - 10**log_e[:-1])



msra = np.array([float(i) for i in mspdata['RAJD'].values])
msdec = np.array([float(i) for i in mspdata['DECJD'].values])
icra = np.array([float(i) for i in icdata['RA[deg]']])
icdec = np.array([float(i) for i in icdata['Dec[deg]']])
icang = np.array([float(i) for i in icdata['AngErr[deg]']])
iceng = np.array([float(i) for i in icdata['log10(E/GeV)']])
global p, lg, lnu
p = len(msra)
lg = len(icra) // p + 1
lnu = len(icra)
upstop_ttt = np.asfarray([uptdata[i]['MJD_stop[days]'].values[-1] for i in range(len(uptdata))])
upstart_ttt = np.asfarray([uptdata[i]['MJD_start[days]'].values[0] for i in range(len(uptdata))])
earea = np.array([eadata[i]['A_Eff[cm^2]'].values for i in range(len(eadata))]) * 1e-4
vec_uptparts = np.asarray(upt_icparts, dtype=np.int64)
upt_icparts = np.asarray(upt_icparts)




deg2rad_var = np.pi/180
@jit(nopython=True, fastmath=True)
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
    
    a = np.add(np.multiply(np.sin(lat1), np.sin(lat2)), np.multiply(np.multiply(np.cos(lat1), np.cos(lat2)), np.cos(dlon)))

    if rad == True:
        return np.arccos(a)
    else:
        return np.abs(np.arccos(a)/deg2rad_var)



@jit(nopython=True, fastmath=True)
def S_ij(nu): 

    '''
    Parameters
    ----------
    nu : int
        Index of the neutrino in the sample
        
    
    Returns
    -------
        Returns the signal PDF for the {psrno}th pulsar and nuind_inp neutrino
    '''
    ang2 = hvovec(msra, msdec, icra[nu], icdec[nu], rad=True) ** 2#, icra[nuind], icdec[nuind], rad=True) ** 2
    sg = np.deg2rad(icang[nu]) ** 2
    return np.divide(np.exp(-1 * np.divide(ang2, 2*sg)), (2 * np.pi * sg))

@jit(nopython=True, fastmath=True)
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




@jit(nopython=True, fastmath=True)
def Bi_stacked(nu, cone=5):

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




all_sig_bag = []

#gamma = 1
all_weights = weights(-1)

sum_weights = [np.sum(i) for i in all_weights]
sum_weights = np.asfarray(sum_weights)
pool = mul.Pool(8, maxtasksperchild=100)
op_async = pool.map_async(sigbag_nu, tqdm(range(lnu)))
tmp = op_async.get()
pool.close()
pool = []
op_async = []
tmp = np.asfarray(tmp)
all_sig_bag.append(tmp)

#gamma = 2
all_weights = weights(-2)

sum_weights = [np.sum(i) for i in all_weights]
sum_weights = np.asfarray(sum_weights)
pool = mul.Pool(8, maxtasksperchild=100)
op_async = pool.map_async(sigbag_nu, tqdm(range(lnu)))
tmp = op_async.get()
pool.close()
pool = []
op_async = []
tmp = np.asfarray(tmp)
all_sig_bag.append(tmp)

#gamma = -2.5
all_weights = weights(-2.5)

sum_weights = [np.sum(i) for i in all_weights]
sum_weights = np.asfarray(sum_weights)
pool = mul.Pool(8, maxtasksperchild=100)
op_async = pool.map_async(sigbag_nu, tqdm(range(lnu)))
tmp = op_async.get()
pool.close()
pool = []
op_async = []
tmp = np.asfarray(tmp)
all_sig_bag.append(tmp)



#gamma = -3
all_weights = weights(-3)

sum_weights = [np.sum(i) for i in all_weights]
sum_weights = np.asfarray(sum_weights)
pool = mul.Pool(8, maxtasksperchild=100)
op_async = pool.map_async(sigbag_nu, tqdm(range(lnu)))
tmp = op_async.get()
pool.close()
pool = []
op_async = []
tmp = np.asfarray(tmp)
all_sig_bag.append(tmp)