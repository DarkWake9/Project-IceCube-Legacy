#PROGRAM TO PLOT THE SKY DISTRIBUTION OF NEUTRINO EVENTS RECORDED BY ICE CUBE BETWEEM 2008-18 IN GALACTIC coordinates
#ONLY "HAMMER" PROJECTIONS ARE PLOTTED!!!!

import pandas as pd
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord as scr
from matplotlib import pyplot as plt
from astroML.plotting import plot_tissot_ellipse
#from pathlib import Path

#folder = Path()
def skydistroof_IC_neu(filenamewithpath):
#read the event csv file and construct a data frame from it
    with open(filenamewithpath, 'r') as f1:
        lines = f1.readlines()

    content=[]
    column=lines[0].split()
    column.pop(0)
    for line in lines[1:]:
        content.append(line.split())

    data = pd.DataFrame(content, columns = column)
    #data = data.drop(0)
    data #ALL THE CONTENTS ARE IN STRING FORMAT NOT FLOAT


    #EDIT 23-05-2022 WITH ASTROPY TO PROJECT IN GALACTIC COORDINATES

    #CONVERTS THE RA, DEC COORDINATES TO GALACTIC COORDINATES

    #THE STRING VALUES IN "data" ARE CONVERTED TO FLOAT VALUES AND STORED IN "r_a" AND "decl"
    #WARNING !!!DONOT USE GENERATORS IN ASTROPY.SKYCOORD! THEY ARE NOT SUPPORTED!!!           Dt: 23-05-2022
    r_a = [float(i) for i in data['RA[deg]']] 
    decl = [float(i) for i in data['Dec[deg]']]     
    radec = scr(ra = r_a * u.degree, dec = decl * u.degree, frame = 'icrs')           
    #WARNING!!!SKYCOORD ONLY SUPPORTS NORMAL LISTS NOT EVEN NUMPY.NDARRAYS!!!    Dt: 23-05-2022
    radec = radec.galactic      #CONVERSION FROM RA,DEC TO GALACTIC COORDINATES


    #PLOT

    fig = plt.figure(figsize=(128,72))
    #plt.rcParams.update({'font.size': 48})

    #ax1 IS FOR HAMMER PROJECTION
    #ax1 = plt.subplot(121,projection='hammer')
    #ax1.scatter(radec.l, data['Dec[deg]'][:],marker='o',color='b', s = 0.7,alpha=0.5)
    #ax1.xaxis.set_major_locator(plt.FixedLocator(np.pi / 3 * np.linspace(-2, 2, 5)))
    #ax1.xaxis.set_minor_locator(plt.FixedLocator(np.pi / 6 * np.linspace(-5, 5, 11)))
    #ax1.yaxis.set_major_locator(plt.FixedLocator(np.pi / 6 * np.linspace(-2, 2, 5)))
    #ax1.yaxis.set_minor_locator(plt.FixedLocator(np.pi / 12 * np.linspace(-5, 5, 11)))
    #ax1.grid(True, which='major', color='r', lw=2)
    #ax1.set_title('HAMMER projection SKY\n\n')

    #ax2 = plt.subplot(122,projection='hammer')     #UNCOMMENT TO COMPARE GALACTIC LAT AND LONG WITH RA,DEC
    
    ax2 = plt.subplot(projection='hammer')
    ax2.scatter(radec.l.degree, radec.b.degree,marker='o',color='b', s = 0.7,alpha=0.5)
    ax2.xaxis.set_major_locator(plt.FixedLocator(np.pi / 3 * np.linspace(-2, 2, 5)))
    ax2.xaxis.set_minor_locator(plt.FixedLocator(np.pi / 6 * np.linspace(-5, 5, 11)))
    ax2.yaxis.set_major_locator(plt.FixedLocator(np.pi / 6 * np.linspace(-2, 2, 5)))
    ax2.yaxis.set_minor_locator(plt.FixedLocator(np.pi / 12 * np.linspace(-5, 5, 11)))
    ax2.grid(True, which='both', color='r', lw=0.8, ls='-.')
    #plot_tissot_ellipse(circ_long[:, None], circ_lat, radius, ax=ax, fc='k', alpha=0.3, linewidth=0)
    ax2.set_title('HAMMER projection \n')
    plt.suptitle(filenamewithpath.replace('icecube_10year_ps/events/','').replace('_exp','') + ' GALACTIC')
    plt.show()

filenames = ["icecube_10year_ps/events/IC40_exp.csv", "icecube_10year_ps/events/IC59_exp.csv",
 "icecube_10year_ps/events/IC79_exp.csv", "icecube_10year_ps/events/IC86_I_exp.csv", "icecube_10year_ps/events/IC86_II_exp.csv",
  "icecube_10year_ps/events/IC86_III_exp.csv", "icecube_10year_ps/events/IC86_IV_exp.csv", "icecube_10year_ps/events/IC86_V_exp.csv",
   "icecube_10year_ps/events/IC86_VI_exp.csv","icecube_10year_ps/events/IC86_VII_exp.csv"]

#To see distribution of a specific event, uncomment the below line and change the index accordingly

#skydistroof_IC_neu(filenames[2])

#To see distributions of all observations uncomment the below lines

for i in filenames:
    skydistroof_IC_neu(i)