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
import scipy.stats as sct
from astropy.coordinates import SkyCoord as scr
from astropy import units as u
####
#### IMPORTING AND SPLITTING ICDATA $$$


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



#IMPORTING MSPDATA
f = open("/media/darkwake/VIB2/Project-IceCube/allpsr.txt", 'r')
lines = f.readlines()

content=[]
column=lines.pop(0).split()
for line in lines[:]:
    content.append(line.split())
    #the INITAL DATABASE IS CLUTTERED SO WE REMOVE THE NULL COLUMNS AND OTHER CLUTTER
mspdata = pd.DataFrame(content).dropna().drop_duplicates()#.drop(range(0,6)).dropna().drop([2,6,8,10,11,13,14], axis=1)
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
msra = [float(i) for i in mspdata['RAJD']]
msdec = [float(i) for i in mspdata['DECJD']]
icra = [float(i) for i in icdata['RA[deg]']]
icdec = [float(i) for i in icdata['Dec[deg]']]
icang = [float(i) for i in icdata['AngErr[deg]']]
#global p, lg
p = len(msra)
lg = len(icra) // p + 1
#cut = 20
#cone = 3
#psrno = 0


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




def angfinder(b, cut):
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
                if temp[tt] > cut:
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
                if temp[tt] > cut:
                    temp[tt] = -1
            ang.extend(temp)
        
    return ang



def S_ij(i, aang, cut):### $P_i[a] = \dfrac{ns[a]}{N} S[a] + \left(1 - \dfrac{ns[a]}{N}\right)B[a]$ if ns 
        
    arr = []
    for j in range(0,len(icdata)):
        ang = aang[j]
        if ang != -1:
            if ang < cut:
                sg = np.deg2rad(icang[j]) ** 2
                ang = np.deg2rad(ang) ** 2
                arr.append(np.exp(-1 * ang / (2 * sg)) / (2 * np.pi * sg))
    return arr


def nsang(aang, cut):
    count = 0
    for j in aang:
        if j != -1:
            if j < cut:
                count += 1
    return count



def bg(i, cone, ang,Ns):
    icdect = [icdec[i] for i in range(len(ang)) if ang[i] != -1]
    #Calculating total no of neutrino events within delta +/- cone of ith PULSAR (not neutrino)
    count = np.count_nonzero(np.abs(np.subtract(msdec[i],icdect))<cone)
    #calculating solid angle with lt deg conesxipy
    sang = 2 * np.pi * (np.sin(np.deg2rad(msdec[i] + cone)) - np.sin(np.deg2rad(msdec[i] - cone)))
    return count/sang


#B = bg(psrno)

def Pr(x, Ns, S, B):
    nsN = x/Ns
    return np.add(np.multiply(nsN , S), np.multiply(np.subtract(1, nsN), B))

def TS(i, x, S, B, Ns):
    return 2*np.sum(np.log(Pr(x,  Ns, S, B)/B))


def TSmax(i, S, B, Ns):
    #returns the TSmax for i^th GRB
    nll = lambda x: -TS(i,x, S, B, Ns)
    soln = minimize(nll, 3 ,bounds=((0,None),))#313, 314, 349, 370, 374, 378, 381, 394, 396, 397, 423, 434
    ns = float(soln.x)
    return ns

#Driving function
def singpsr2(psrno = 0):
    cut = 20
    cone = 3
    aang = angfinder(psrno, cut)
    Ns = nsang(aang, cut)
    #Ns = len(icra)
    S = S_ij(psrno, aang, cut)
    B = bg(psrno,cone, aang, Ns)
    #TS(psrno, 0)
    maxns = TSmax(psrno, S, B, Ns)
    return [maxns, np.sqrt(TS(psrno, maxns, S, B, Ns))]


pool = mul.Pool(12, maxtasksperchild=800)
op_async = pool.map_async(singpsr2, range(p))
tsmax = op_async.get()
op_async = []
pool = []
maxns = []
sqrtts = []
for i in tsmax:
    maxns.append(i[0])
    sqrtts.append(i[1])


np.savetxt('tsmax.txt', sqrtts)
np.count_nonzero(maxns)
binnedsts, wsts = np.histogram(sqrtts, 15)

ne01 = 4
amp  = np.sum(binnedsts)/ne01
#binnedsts = binnedsts/amp
yerr = np.sqrt(binnedsts)
b2 = binnedsts[(binnedsts>0)]
yerr2 = yerr[(binnedsts>0)]
x=0.5*(wsts[1:]+ wsts[:-1])
x2 = x[(binnedsts>0)]
xerr=x[0]-wsts[0]
gaussian = lambda x, mean, var:  max(b2)*norm.pdf(x, mean, var)

gparam, gerr = curve_fit(gaussian, x2, b2,  p0=[0.1, 1.0],sigma=yerr2, absolute_sigma=True, maxfev = 100000, method='dogbox')
gx = np.linspace(0, max(wsts), 100)

font = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 21,
        }
plt.figure(figsize=(10.8,7.2))
#plt.xticks([0, 1, 2, 3, 4], [0, 1, 2, 3, 4],fontsize=18)
plt.yscale('log')
#plt.yticks([1, 10, 100, 1000, 10000], [1, 10, 100, 1000, 10000], fontsize=18)

#plt.plot(gx,gaussian(gx, gparam[0], gparam[1]),'--', alpha=0.6, label=f'Gaussian Fit', color='darkred', lw=2.2)

plt.errorbar(x, binnedsts, xerr=xerr,yerr=yerr,  aa = True, barsabove = True, fmt='ob',ecolor='blue', elinewidth=1.35, label='$\sqrt{TS_{max}}$')
plt.xlabel("$\mathbf{\sqrt{TS_{max}}}$", fontweight='bold', fontdict=font)
plt.ylabel("No.of Events", fontweight='bold', fontdict=font)
#plt.grid(True, which='major')
legend =plt.legend(prop={'size':17}, framealpha=0.0)
legend.get_frame().set_facecolor('none')
#plt.savefig("outputs/sqrt(TSmax)err.pdf")
plt.show()



plt.hist(sqrtts, 15)
plt.xlabel("$\mathbf{\sqrt{TS_{max}}}$", fontweight='bold', fontdict=font)
plt.ylabel("No.of Events", fontweight='bold', fontdict=font)
plt.show()