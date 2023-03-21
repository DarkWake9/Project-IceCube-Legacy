from numba import jit, njit
import numpy as np
from scipy.optimize import minimize
from core import readfiles

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

@njit
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
    t_upt = upstop_ttt[nusample_wall] - upstart_ttt[nusample_wall]
    t_upt*=86400
    d_ind = 0

    for i in range(len(dec_nu) - 1):
        if dec_nu[i] <= psr_decl and dec_nu[i+1] >= psr_decl:
            d_ind = i
            break
    ea_temp = earea[nusample_wall][d_ind*40:(d_ind+1)*40]
    
    #earea = np.array(ea_temp['A_Eff[cm^2]'].values) * 1e-4
    
    weight_kj = t_upt * np.sum(np.multiply(np.multiply(ea_temp, np.power(e_nu, gamma)), de_nu))
    
    #sleep(1e-8)
    
    return weight_kj   
                       
            
class weights:
    def __init__(self, gamma):
        self.gamma = gamma
        self.all_weights = []
        for i in range(len(upt_icparts)-1):
            self.all_weights.append([psr_wt_quick(i, psrno, gamma) for psrno in range(p)])
        self.all_weights = np.asfarray(self.all_weights)
        self.sum_weights = [np.sum(i) for i in self.all_weights]
        self.sum_weights = np.asfarray(self.sum_weights)


