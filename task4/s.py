
import numpy as np
import matplotlib.pylab as plt
from math import pi
from scipy.optimize import minimize
data=np.genfromtxt("task4/nu_near_p57.txt", delimiter=',',names=True)


nura=data['nura'] #neutrino RA
nudec=data['nudec'] #neutrino DEC
angerr=data['angerr'] # error in neutrino position
ra=267.02 #pulsar RA
dec=-24.77917 # pulsar DEC
cang=np.sin(dec*pi/180.0)*np.sin(nudec*pi/180.0)+ np.cos(dec*pi/180.0)*np.cos(nudec*pi/180.0)*np.cos((ra-nura)*pi/180.0) 
ang=(180.0/pi)*np.arccos(cang) # angle between pulsar and neutrino in degrees
angerrrad=angerr*pi/180.0 # angle between pulsar and neutrino in radins
N=len(nura) # total no of neutrinos  used for analysis
bkgct=np.count_nonzero(np.abs(nudec-dec)<3) # total no of bkgd neutrinos within 3 degrees
Bi=bkgct/(2*pi*(np.sin((pi)*(dec+3)/180.0)- np.sin(pi*(dec-3)/180)))
Si=np.exp(-0.5*(ang/angerr)**2)/(2*pi*angerrrad*angerrrad)


# the function rtruns TS(ns) in Zhou paper
def L(x):
    TS=np.log(x*Si/N+(1-x/N)*Bi)-np.log(Bi)
    return 2*np.sum(TS)

# since scipy.optimize can only do minization calculate negative of L(x) to maximize L(x)
nll = lambda x: -L(x)
soln = minimize(nll, 3,bounds=((0,None),)) # to ensure that ns is > 0
ns=soln.x
print(ns)
print(L(ns))
plt.figure()
x = np.array(range(1000))
plt.plot(x, L(x))
plt.show()
