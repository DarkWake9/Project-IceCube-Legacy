####
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import os
from astropy import units as u
from astropy.coordinates import SkyCoord as scr
import seaborn as sns
import scipy as sp
import load

#### [markdown]
# # OBJECTIVES

#### [markdown]
# Split the data energy wise 
# ## (0.3-1 TeV, 1-3 TeV, 3-10 TeV)

#### [markdown]
#  Redo both sets of experiments (Task 1 and Task 2)

#### [markdown]
# For method 2 (when counting no of matches within the neutrino error circle) instead of
# counting the background events by looking at the no of neutrino-MSP pairs within the same cos(error region) offset by 5 degrees,
# 
#     Follow the same method as in        https://arxiv.org/pdf/2112.11375v1.pdf      which is briefly summarized below
# 
#         1. Generate a synthetic catalog of MSPs uniformly distributed in the same RA and DEC range as observed population
# 
#         2, Count no of coincidences  within neutrino error circle for the synthetic catalog in (1). 
# 
#         3. Repeat step (1) and (2) 100 times and take the average.

#### [markdown]
# ## IMPORTING DATA
# 

#### [markdown]
# 
# ### SPLITTING THE DATA

####
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

###
icdata

###
icpartitions = (np.log10(300), np.log10(1000), np.log10(3000), np.log10(10000))
lwall = [0]
rwall = []
for i in range(0,len(icpartitions)):
    lwall.append(min(icdata[icdata['log10(E/GeV)'] > icpartitions[i]].index))
    rwall.append(max(icdata[icdata['log10(E/GeV)'] <= icpartitions[i]].index))
lwall.pop(-1)

###
ic = [icdata[lwall[i]:rwall[i] + 1] for i in range(len(lwall))]

####
ic.pop(0)
for i in range(3):
    ic[i] = ic[i].reset_index()
    ic[i] = ic[i].rename(columns={'index': 'oindex'})

extensions = [441 - len(ic[i])%441 for i in range(3)]

icra = [[float(i) for i in ic[j]['RA[deg]']] for j in range(3)]
icdec = [[float(i) for i in ic[j]['Dec[deg]']] for j in range(3)]
icang = [[float(i) for i in ic[j]['AngErr[deg]']] for j in range(3)]

for i in range(3):
    icra[i].extend([0]*extensions[i])
    icdec[i].extend([0]*extensions[i])
    icang[i].extend([0]*extensions[i])

#IMPORTING AND CLEANING MSPDATA
f = open("/media/darkwake/VIB2/Project-IceCube/10milsecpsr.txt", 'r')
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
mspdata.dropna(inplace=True)
mspdata = mspdata.reset_index()
mspdata = mspdata.drop(["index", "#"], axis=1)


msra = [float(i) for i in mspdata['RAJD']]
msdec = [float(i) for i in mspdata['DECJD']]
cos5 = np.cos(np.deg2rad(5))
minra = min(msra)
maxra = max(msra)
mindec = min(msdec)
maxdec = max(msdec)
titles = ['0.3TeV - 1 TeV', '1 TeV - 3 TeV', '3 TeV - 10 TeV']


output3 = []
for i in range(100):
# #### GENERATING SYNTHTIC CATALOGUE
# ####
    msra2, msdec2 = load.gen_cat(minra, maxra, mindec, maxdec)

# #### TASK 2B OUTPUT
# #### COSINE OF SPACE ANGLE DISTRIBUTION
    cosp_ang = load.t2a2b(icra, icdec, msra2, msdec2, extensions, 'b')
    freq, nbins, bground, signal = load.t2bs(cosp_ang)
    cosp_ang=[]

# ## TASK 2C
# ### OUTPUT

    match = load.t2c(icra,icdec, icang, msra2, msdec2, extensions)

    '''for i in range(3):
        print(f"{titles[i]}: {match[i]}")'''

# ## TASK 2D
# ### OUTPUT

    background = load.t2d(icra, icdec, icang, msra2, msdec2, extensions)

    '''for i in range(3):
        print(f"{titles[i]}: {background[i]}")'''
    output3.append([freq, bground, signal, match, background])
