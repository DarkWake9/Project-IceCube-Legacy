############################################################################################################################

# This file contains the function to calculate the weights of each pulsar 
# for each neutrino season/sample for a given spectral index gamma
# The function is called psr_wt_quick
# 
# IT FOLLOWS EQN 8 OF 2205.15963
# 
# The class weights is used to store the weights of all pulsars for all neutrino seasons/samples
# For a given gamma an object of class weights is created containing 
# the weights of all pulsars for all neutrino samples and their sum

############################################################################################################################



from numba import njit, prange
import numpy as np
from scipy.optimize import minimize
from core import readfiles

from core.req_arrays import *



@njit
def psr_wt_quick(nusample_wall, psrno, gamma = 1):

    '''
    EQN 8 OF 2205.15963
    -------------------

    Parameters
    ----------
    nusample_wall : float
        The index of catalogue of the given neutrino ===> season
    psrno : int
        The index of the pulsar in the ATNF catalogue
    gamma : float, optional
        (default) 1
        The spectral index of the power law

    Returns
    -------
    
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

    # THIS LINE FOLLOWS EQN 8 OF 2205.15963
    weight_kj = t_upt * np.sum(np.multiply(np.multiply(ea_temp, np.power(e_nu, gamma)), de_nu))
    #sleep(1e-8)
    return weight_kj   
                       
            
class weights:
    def __init__(self, gamma):
        self.gamma = gamma
        all_weights = []
        for i in prange(len(upt_icparts)-1):
            all_weights.append([psr_wt_quick(i, psrno, gamma) for psrno in range(p)])
        self.all_weights = np.asfarray(all_weights)
        sum_weights = [np.sum(i) for i in self.all_weights]
        self.sum_weights = np.asfarray(sum_weights)


