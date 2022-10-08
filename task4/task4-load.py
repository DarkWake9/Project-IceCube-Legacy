import numpy as np
import pandas as pd
import os
import multiprocessing as mul
from multiprocessing import Process
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize

#### [markdown]
# ## IMPORTING AND SORTING THE DATA
# 
# WE HAVE CHOSEN ICECUBE 2008-18 DATA FOR THIS STUDY
# 
# 
# THE ICE CUBE DATA SET HAS 1134450 NEUTRINO EVENTS
# 
# 
# We select neutrino events with Energy >= 100TeV = 10^5 GeV
# 
# i.e log10(E/GeV) > 5
# 
# 
# There are 192107 such neutrino events in this data
# 
# The ms pulsars are taken from the ATNF Catalogue
# 
# There are 441 pulsars (as of May 2022 when the study started)
# 
# All the pulsars lie in the declination range of -87 to +87 degrees

####
####
#### IMPORTING AND SPLITTING ICDATA $$$


def import_icdata():
    path = "/media/darkwake/VIB2/Project-IceCube/icecube_10year_ps/events"
    filenames = ["IC40_exp.csv", "IC59_exp.csv","IC79_exp.csv", "IC86_I_exp.csv", "IC86_II_exp.csv",
    "IC86_III_exp.csv", "IC86_IV_exp.csv", "IC86_V_exp.csv", "IC86_VI_exp.csv", "IC86_VII_exp.csv"]
    file = filenames[0]
    f = open(os.path.join(path, file), 'r')
    lines = f.readlines()
    column=lines[0].split()
    column.pop(0)
    content = []
    for file in filenames:
        f = open(os.path.join(path, file), 'r')
        lines = f.readlines()
        for line in lines[1:]:
            content.append(line.split())
        f.close()
    icdata = pd.DataFrame(content, columns=column)
    icdata['log10(E/GeV)'] = [float(i) for i in icdata['log10(E/GeV)']]

    icdata = icdata.sort_values('log10(E/GeV)')
    icdata = icdata.reset_index()
    icdata = icdata.drop('index', axis=1)
    return icdata

#icdata = icdata.sort_values('log10(E/GeV)')
#icdata = icdata.reset_index()
#icdata = icdata.drop('index', axis=1)
#icdata2 = icdata[icdata['log10(E/GeV)'] > 5]
#icdata2 = icdata2.reset_index()
#icdata2 = icdata2.drop('index', axis=1)
#icdata2
def nu_10gev(icdata):
    icdata2 = icdata[icdata['log10(E/GeV)'] > 5]
    icdata2 = icdata2.reset_index()
    icdata2 = icdata2.drop('index', axis=1)
    return icdata2
#IMPORTING MSPDATA
def import_psrdata(path = "/media/darkwake/VIB2/Project-IceCube/10milsecpsr.txt"):
    f = open(path, 'r')
    lines = f.readlines()

    content=[]
    column=lines[3].split()
    for line in lines[:]:
        content.append(line.split())
        #the INITAL DATABASE IS CLUTTERED SO WE REMOVE THE NULL COLUMNS AND OTHER CLUTTER
    mspdata = pd.DataFrame(content).drop(range(0,6)).dropna().drop([2,6,8,10,11,13,14], axis=1)
    f.close()
    line = []
    lines = []

    mspdata.columns = column
    column = []
    content=[]
    mspdata = mspdata.sort_values('DECJD')
    mspdata.dropna(inplace=True)
    mspdata = mspdata.reset_index()
    mspdata = mspdata.drop(["index", "#"], axis=1)
    return mspdata


def hvovec(lon1, lat1, lon2, lat2):

    #Convert decimal degrees to Radians:
    lon1 = np.deg2rad(lon1)
    lat1 = np.deg2rad(lat1)
    lon2 = np.deg2rad(lon2)
    lat2 = np.deg2rad(lat2)

    #Implementing Haversine Formula: 
    dlon = np.subtract(lon2, lon1)
    #dlat = np.subtract(lat2, lat1)

    a = np.add(np.multiply(np.sin(lat1), np.sin(lat2)), np.multiply(np.multiply(np.cos(lat1), np.cos(lat2)), np.cos(dlon)))

    return np.abs(np.rad2deg(np.arccos(a)))
'''
####
msra = [float(i) for i in mspdata['RAJD']]
msdec = [float(i) for i in mspdata['DECJD']]
icra = [float(i) for i in icdata2['RA[deg]']]
icdec = [float(i) for i in icdata2['Dec[deg]']]
icang = [float(i) for i in icdata2['AngErr[deg]']]
p = len(msra)
lg = len(icra) // p + 1
'''
####
def angfinder(b):
    ang = []
    for a in range(lg):
        
        if a != lg - 1:
        #try:
            ilo = icra[a*p:a*p + p]
            ila = icdec[a*p:a*p + p]
            lo = msra[b] * np.ones(p)
            la = msdec[b] * np.ones(p)
            temp = hvovec(ilo, ila, lo, la)
            for tt in range(len(temp)):
                if temp[tt] > 20:
                    temp[tt] = -1
            ang.extend(temp)
        else:
        #except:
            ilo = icra[a*p:]
            ila = icdec[a*p:]
            ext = len(ilo)
            lo = msra[b] * np.ones(ext)
            la = msdec[b] * np.ones(ext)
            temp = hvovec(ilo, ila, lo, la)
            #ang.extend(hvovec(ilo, ila, lo, la))
            for tt in range(len(temp)):
                if temp[tt] > 20:
                    temp[tt] = -1
            ang.extend(temp)
            
        
    return ang

pool = mul.Pool()
op_async = pool.map_async(angfinder, range(0,p))
aang = op_async.get()
op_async = []
pool = []

#### [markdown]
# we only select those neutrino events, whose angular distance  is within $20^{\circ}$ of a MSP.  For a dataset of $N$ events, if $n_s$ signal events are attributed to a millisecond pulsar, the probability density of an individual event $i$ is given by:
# \begin{equation}
# P_i = \frac{n_s}{N} S_i + (1-\frac{n_s}{N}) B_i,
# %\label{eq1}
# \end{equation}
# where $S_i$ and $B_i$ represent the signal and background pdfs, respectively.
# The likelihood function ($\mathcal{L}$) of the entire dataset, obtained from the product of each individual PDF can be written as 
# \begin{equation}
# \mathcal{L} (n_s) = \prod_{i=1}^N P_i,
# \end{equation}
# where $P_i$ is the same as in eq1 above. The signal PDF is given by:
# \begin{equation}
# \large
# S_i = \frac{1}{2\pi\sigma_i^2}e^{\frac{-(|\theta_i-\theta_s|)^2}{2\sigma_i^2}}
# \end{equation}
# where $|\theta_i-\theta_s|$ is the angular distance between the  MSP and the neutrino, whereas $\sigma_i$ is the angular uncertainty in the neutrino position.
# The background PDF is determined from the average number of events per solid angle  within a declination of $\pm 3^{\circ}$ after averaging over RA. We do not incorporate  the energy information, since the pubic IceCube catalog only contains the energy of the reconstructed muon. The detection statistic used to ascertain the presence of a signal is given by:
# \begin{equation}
# TS (n_s) = 2 \log \frac{\mathcal{L} (n_s)}{\mathcal{L} (0)}
# \end{equation}
#  If the null hypothesis is true, $TS (n_s)$ obeys the $\chi^2$ distribution for one degree of freedom. We calculated $TS (n_s)$ for all MSPs in the  catalog.

#### [markdown]
# ### FOR $\nu_i$ and $p_j$ and  N = len(icdata2) = Total no.of $\nu$ samples
# ### ns[i] = no.of $\nu$ events with angles $<= 20 \deg$ with $psr_i$
# ### $S_i = \frac{1}{2\pi\sigma_i^2}e^{\frac{-(|\theta_i-\theta_s|)^2}{2\sigma_i^2}}$
# ### $|\theta_i-\theta_s|$ = aang[i][s]
# ### $\sigma_i$ = icang[i]           and                  sg =      $\sigma_i^2$
# ###

####
N = len(icdata2)

ns = []
for i in aang:
    count = 0
    for j in i:
        if j != -1:
            count += 1
    ns.append(count)
ns = np.array(ns)
    

#### [markdown]
# Ns = np.sum(ns)//2

####
N

#### [markdown]
# n2 = []
# for i in aang:
#     count = 0
#     for j in i:
#         if j != -1:
#             if j < 5:
#                 count += 1
#     n2.append(count)

#### [markdown]
# for i in aang[0]:
#     if i < 5 and i!= -1:
#         print(aang[0].index(i))

#### [markdown]
# S_ij(10754,0)

#### [markdown]
# [msra[0], msdec[0]]

#### [markdown]
# [icra[10754], icdec[10754], icang[10754]]

#### [markdown]
# [icra[119640], icdec[119640], icang[119640]]
# 

#### [markdown]
# S_ij(119640,0)

#### [markdown]
# for i in aang[0]:
#     if i < 6 and i!= -1 and i >= 5:
#         print(aang[0].index(i))

#### [markdown]
# S_ij(31410, 0)

#### [markdown]
# for i in aang[0]:
#     if i < 7 and i!= -1 and i >= 6:
#         print(aang[0].index(i))

#### [markdown]
# np.exp(-0.05)/150

#### [markdown]
# aang[0][39160] **

#### [markdown]
# S_ij(39160,0)

#### [markdown]
# for i in aang[0]:
#     if i < 8 and i!= -1 and i >= 7:
#         print(aang[0].index(i))

#### [markdown]
# S_ij(99351,0)
# #aang[0][142]

#### [markdown]
# ang = np.deg2rad(aang[0][99351]) ** 2
# 
#     #if ang < 20 and ang != -1:
# sg = np.deg2rad(icang[99351]) ** 2
# np.exp(-1 * ang / (2 * sg)) / (2 * np.pi * sg)

#### [markdown]
# ang = np.deg2rad(aang[0][99351]) ** 2
# err = np.deg2rad(icang[99351]) **2
# np.exp(-ang/(2 * err))/(2 * np.pi * err)

#### [markdown]
# np.exp(-ang/(2 * err))

#### [markdown]
# -ang/(2 * err)

#### [markdown]
# err

#### [markdown]
# a2 = aang[0][99351] ** 2
# a2

#### [markdown]
# icang[99351]

####

#S_ij for ith neutrino and jth pulsar IS WRONG
#EDIT 04102022 - 11.28 - S_ij is for i^th pulsar and j^th neutrino and summed over all NEUTRINOS
def S_ij(i):
    arr = []
    
    for j in range(0,len(icdata2)):
        sg = np.deg2rad(icang[j]) ** 2
        ang = aang[i][j]
        if ang != -1:
            if ang < 20:
                ang = np.deg2rad(aang[i][j])
                ang = ang ** 2
                arr.append(np.exp(-ang / (2 * sg)) / (2 * np.pi * sg))
        #return
        else:
            arr.append(0)
        #return -1
    
    return np.sum(arr)
def sij():
    pool = mul.Pool()
    op_async = pool.map_async(S_ij, range(0,p))
    si = op_async.get()
    pool = []
    op_async = []
    return si
#for i in range(len(si)):

#S = Si()

####
S = np.array(sij()) 

#### [markdown]
# ### S stores $S_i$ values for all pulsars

#### [markdown]
# 2103.12813
# 
# For each sample, given the Dec δi of an event, the background PDF is determined by the relative number of
# events in δi ± $3^◦$ divided by the solid angle.
# 
# so calculate total no of events within delta +/- 3
# and then divide by 2 pi (sin[delta+3]-sin[delta-3])
# 
# you can choose all events
# withiin delta +/- 3
# 

####
def bg(i):
    #Calculating total no of neutrino events within delta +/- 6 of ith PULSAR (not neutrino)
    count = 0
    for j in icdec:
        if abs(msdec[i] - j) <= 6:
            count+=1 
    #calculating solid angle with 3deg
    sang = 2 * np.pi * (np.sin(np.deg2rad(msdec[i] + 6)) - np.sin(np.deg2rad(msdec[i] - 6)))
    return count/sang

def Bi():
    pool = mul.Pool()
    op_async = pool.map_async(bg, range(0,p))
    return op_async.get()

B = np.array(Bi())

#### [markdown]
# ### B stores $B_i$ values for all pulsars

####
B

#### [markdown]
# ### For $a^{th}$ pulsar

#### [markdown]
# ### $P_i[a] = \dfrac{ns[a]}{N} S[a] + \left(1 - \dfrac{ns[a]}{N}\right)B[a]$ if ns = no.of $\nu$ events with angles $<= 20 \deg$ with $psr_a$
# 
#                                                  OR
#                                                  
# ### $P_i[a] = \dfrac{ns}{N} S[a] + \left(1 - \dfrac{ns}{N}\right)B[a]$     where ns is a scalar number specified manually in
# ### $\mathcal{L}(ns) = np.product(P_i(ns, N, S, B))$

####
for i in range(p):
    if B[i] == 0:
        print(i)

#### [markdown]
# sssss = 0
# for i in B:
#     if i == 0:
#         sssss+=1
# sssss

#### [markdown]
# len(B)

####
def Pr(ns, N, S, B):
    nsN = ns/N
    return np.add(np.multiply(nsN , S), np.multiply(np.subtract(1, nsN), B))

####


####
n2 = len(icdata2)
a = ns[0]/n2

####
a * S[0] + (1-a)*B[0]

####
[S[0], B[0]]

####
1 - a

####
0 in l

####
def Ls(ns1, N,S, B):
    return np.prod(Pr(ns1, N,S, B))

#### [markdown]
# def L_n(S, B, N, ns):
#     L = np.ones(len(ns))
#     for i in range(0, len(l)):
#         for j in range(len(l[i])):
#             L[i] *= l[i][j]
#     return L
# L = L_n(S, B, N, ns)
# 

####
def TS(ns, S, B, N):
    ts = []
    for i in ns:
        ts.append(2 * np.log10(Ls(i, N, S, B)/Ls(0, N, S, B)))
    return ts
ts = TS(ns, S, B, N)






####
plt.figure(figsize=(8,6))
plt.yscale('log')
plt.plot(list(ts), list(ns))
plt.show()

####



