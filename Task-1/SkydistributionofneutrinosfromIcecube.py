#PROGRAM TO PLOT THE SKY DISTRIBUTION OF NEUTRINO EVENTS RECORDED BY ICE CUBE BETWEEM 2008-18
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from astroML.plotting import plot_tissot_ellipse


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

# generate a latitude/longitude grid

    circ_long = np.linspace(-np.pi, np.pi, 13)[1:-1]
    circ_lat = np.linspace(-np.pi / 2, np.pi / 2, 7)[1:-1]
    radius = 10 * np.pi / 180.

#plot the sky distributions in Hammer, Aitoff, Mollweide, Lambert projections

    fig = plt.figure(figsize=(64,48))
    plt.rcParams.update({'font.size':15})
    
    for (i, projection) in enumerate(['Hammer', 'Aitoff', 'Mollweide', 'Lambert']):
        ax = plt.subplot(221 + i, projection=projection.lower())
        ax.scatter(data['RA[deg]'][:], data['Dec[deg]'][:],marker='o',color='b', s = 0.2 ,alpha=0.5)
        ax.xaxis.set_major_locator(plt.FixedLocator(np.pi / 3 * np.linspace(-2, 2, 5)))
        ax.xaxis.set_minor_locator(plt.FixedLocator(np.pi / 6 * np.linspace(-5, 5, 11)))
        ax.yaxis.set_major_locator(plt.FixedLocator(np.pi / 6 * np.linspace(-2, 2, 5)))
        ax.yaxis.set_minor_locator(plt.FixedLocator(np.pi / 12 * np.linspace(-5, 5, 11)))
        ax.grid(True, which='minor', color='r', ls='-.', lw=1)

#    plot_tissot_ellipse(circ_long[:, None], circ_lat, radius, ax=ax, fc='k', alpha=0.3, linewidth=0)
    
        ax.set_title('%s projection' % projection)
        plt.suptitle(filenamewithpath.replace('icecube_10year_ps/events/',''))
    plt.show()

#All observations from IceCube from 2008-18

filenames = ["icecube_10year_ps/events/IC40_exp.csv", "icecube_10year_ps/events/IC59_exp.csv",
 "icecube_10year_ps/events/IC79_exp.csv", "icecube_10year_ps/events/IC86_I_exp.csv", "icecube_10year_ps/events/IC86_II_exp.csv",
  "icecube_10year_ps/events/IC86_III_exp.csv", "icecube_10year_ps/events/IC86_IV_exp.csv", "icecube_10year_ps/events/IC86_V_exp.csv",
   "icecube_10year_ps/events/IC86_VI_exp.csv","icecube_10year_ps/events/IC86_VII_exp.csv"]

#To see distribution of a specific event, uncomment the below line and change the index accordingly

#skydistroof_IC_neu(filenames[0])

#To see distributions of all observations uncomment the below lines

for i in filenames:
    skydistroof_IC_neu(i)