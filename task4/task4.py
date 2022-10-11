from json import load
import numpy as np
import matplotlib.pylab as plt
from math import pi
from scipy.optimize import minimize
import task4ld
import os
print(os.getcwd())
data=np.genfromtxt("task4/outputs/nunear4.txt")


nura=data[:,0] #neutrino RA
nudec=data[:,1] #neutrino DEC
angerr=data[:,2] # error in neutrino position
ra=282.5 #pulsar RA
dec= -0.03 # pulsar DEC
cang=np.sin(dec*pi/180.0)*np.sin(nudec*pi/180.0)+ np.cos(dec*pi/180.0)*np.cos(nudec*pi/180.0)*np.cos((ra-nura)*pi/180.0) 
ang=(180.0/pi)*np.arccos(cang) # angle between pulsar and neutrino in degrees
angerrrad=angerr*pi/180.0 # angle between pulsar and neutrino in radins
N=len(nura) # total no of neutrinos  used for analysis
bkgct=np.count_nonzero(np.abs(nudec-dec)<3) # total no of bkgd neutrinos within 3 degrees
Bi=bkgct/(2.0*pi*N*(np.sin((pi/180)*(dec+3.0))- np.sin((pi/180)*(dec-3.0))))
Si=np.exp(-0.5*(ang/angerr)**2)/(2.0*pi*angerrrad*angerrrad)
# the function rtruns TS(ns) in Zhou paper
def L(x):
    TS=np.log(x*Si/N+(1.0-1.0*x/N)*Bi)-np.log(Bi)
    return 2.0*np.sum(TS)



# since scipy.optimize can only do minization calculate negative of L(x) to maximize L(x)
nll = lambda x: -L(x)
soln = minimize(nll, 5.0,bounds=((0,None),)) # to ensure that ns is > 0
ns=soln.x
print(L(ns))
print(ns)
x=np.linspace(0,600,10)
y=[L(p) for p in x]
plt.plot(x,y,'-')
plt.show()