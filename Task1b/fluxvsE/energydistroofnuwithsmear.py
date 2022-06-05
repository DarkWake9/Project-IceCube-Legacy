from os import minor
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
    
def energydistroof_IC_neu(filenamewithpath):
    with open("icecube_10year_ps/irfs/IC59_smearing.csv", 'r') as f2:
        lines = f2.readlines()
    content=[]
    column=lines[0].split()
    column.pop(0)
    for line in lines[1:]:
    
        content.append(line.split())

    smear = pd.DataFrame(content, columns = column)
    smear

    de = [pow(10,float(smear['log10(E_nu/GeV)_max'][i])) - pow(10,float((smear['log10(E_nu/GeV)_min'][i]))) for i in range(0,len(smear['log10(E_nu/GeV)_max']))]
    dphi = [float(smear['Dec_nu_max[deg]'][i]) - float((smear['Dec_nu_min[deg]'][i])) for i in range(0,len(smear['Dec_nu_max[deg]']))]
    #print(de)
    #print(de)

    e = [((pow(10,float(smear['log10(E_nu/GeV)_max'][i])) + pow(10,float((smear['log10(E_nu/GeV)_min'][i]))))/2.0) for i in range(0,len(smear['log10(E_nu/GeV)_max']))]
    e2 = [pow((pow(10,float(smear['log10(E_nu/GeV)_max'][i])) + pow(10,float((smear['log10(E_nu/GeV)_min'][i]))))/2.0,2) for i in range(0,len(smear['log10(E_nu/GeV)_max']))]
    e2dphide = [e2[i] * dphi[i]/de[i] for i in range(0,len(smear['log10(E_nu/GeV)_max']))]


    plt.figure(figsize=(26,16))
    plt.xlabel("log(E_nu/GeV)")
    plt.ylabel("Flux E^2dphi/dE")
    plt.title("Energy distribution diagram")
    plt.plot(np.log10(e),e2dphide)
    plt.suptitle(filenamewithpath.replace('icecube_10year_ps/irfs/',''))
    plt.grid(True, which='both')
    plt.show()


filenames = ["icecube_10year_ps/irfs/IC40_smearing.csv", "icecube_10year_ps/irfs/IC59_smearing.csv",
 "icecube_10year_ps/irfs/IC79_smearing.csv", "icecube_10year_ps/irfs/IC86_I_smearing.csv", "icecube_10year_ps/irfs/IC86_II_smearing.csv",
  "icecube_10year_ps/irfs/IC86_III_smearing.csv", "icecube_10year_ps/irfs/IC86_IV_smearing.csv", "icecube_10year_ps/irfs/IC86_V_smearing.csv",
   "icecube_10year_ps/irfs/IC86_VI_smearing.csv","icecube_10year_ps/irfs/IC86_VII_smearing.csv"]


#To see distribution of a specific event, uncomment the below line and change the index accordingly

#energydistroof_IC_neu(filenames[0])

#To see distributions of all observations uncomment the below lines

for i in filenames:
    energydistroof_IC_neu(i)