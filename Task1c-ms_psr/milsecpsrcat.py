#PROGRAM TO PLOT THE MILLISECOND PULSARS IN THE GALACTIC COORDINATES
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord as scr

#OPENING THE REQUIRED FILE AND CONSTRUCTING THE DATAFRAME
with open("Task1c-ms_psr/10milsecpsr.txt", 'r') as f:
    lines = f.readlines()

content=[]
column=lines[3].split()

for line in lines[:]:
    content.append(line.split())

#the INITAL DATABASE IS CLUTTERED SO WE REMOVE THE NULL COLUMNS AND OTHER CLUTTER
ms_psr_data = pd.DataFrame(content).drop(range(0,6)).dropna().drop([2,6,8,10,11,13,14], axis=1)
f.close()

ms_psr_data.columns = column


#THE STRING VALUES IN "data" ARE CONVERTED TO FLOAT VALUES AND STORED IN "r_a" AND "decl"
    #WARNING !!!DONOT USE GENERATORS IN ASTROPY.SKYCOORD! THEY ARE NOT SUPPORTED!!!        
r_a = [float(i) for i in ms_psr_data['RAJD']]
decl = [float(i) for i in ms_psr_data['DECJD']]
#!!!DONOT USE GENERATORS IN ASTRPY.SKYCOORD THEY ARE NOT SUPPORTED!!!           
#!!!SKYCOORD ONLY SUPPORTS NORMAL LISTS NOT EVEN NUMPY.NDARRAYS!!!
radec=scr(ra = r_a * u.degree, dec = decl * u.degree, frame = 'icrs')
#CONVERT TO GALACTIC COORDINATES
radec = radec.galactic

#WE ALSO DOWNLOADED THE GALACTIC COORDINATES FROM THE DATABASE
gl = [float(i) for i in ms_psr_data['Gl']]
gb = [float(i) for i in ms_psr_data['Gb']]

#PLOT
fig = plt.figure(figsize=(32,18))
#plt.rcParams.update({'font.size': 48})
ax2 = plt.subplot(121,projection='hammer')
ax2.scatter(radec.l.degree, radec.b.degree,marker='o',color='b', s = 7,alpha=0.5)
ax2.xaxis.set_major_locator(plt.FixedLocator(np.pi / 3 * np.linspace(-2, 2, 5)))
ax2.xaxis.set_minor_locator(plt.FixedLocator(np.pi / 6 * np.linspace(-5, 5, 11)))
ax2.yaxis.set_major_locator(plt.FixedLocator(np.pi / 6 * np.linspace(-2, 2, 5)))
ax2.yaxis.set_minor_locator(plt.FixedLocator(np.pi / 12 * np.linspace(-5, 5, 11)))
ax2.grid(True, which='both', color='r', lw=0.4, ls = '-.')
ax2.set_title('HAMMER Projection: GALACTIC Coordinates\n Calculated\n')

ax1 = plt.subplot(122,projection='hammer')
ax1.scatter(gl, gb, marker='o',color='b', s = 7,alpha=0.5)
ax1.xaxis.set_major_locator(plt.FixedLocator(np.pi / 3 * np.linspace(-2, 2, 5)))
ax1.xaxis.set_minor_locator(plt.FixedLocator(np.pi / 6 * np.linspace(-5, 5, 11)))
ax1.yaxis.set_major_locator(plt.FixedLocator(np.pi / 6 * np.linspace(-2, 2, 5)))
ax1.yaxis.set_minor_locator(plt.FixedLocator(np.pi / 12 * np.linspace(-5, 5, 11)))
ax1.grid(True, which='both', color='r', lw=0.4, ls = '-.')
ax1.set_title('HAMMER Projection: GALACTIC Coordinates\n Downloaded\n')

plt.suptitle("MILLISECOND PULSAR from ATNF psrcat")


plt.show()