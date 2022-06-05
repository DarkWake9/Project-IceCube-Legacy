
import pandas as pd
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt

#FUNCTION TO PLOT HISTOGRAMS OF N VS E BINNED IN LOG E

def n_e_hist(file):
    with open(file, 'r') as f2:
        lines = f2.readlines()
    content=[]
    column=lines[0].split()
    column.pop(0)
    for line in lines[1:]:
        
        content.append(line.split())


    smear = pd.DataFrame(content, columns = column)
    smear
    f2.close()
    de = [pow(10,float(smear['log10(E_nu/GeV)_max'][i])) - pow(10,float((smear['log10(E_nu/GeV)_min'][i]))) for i in range(0,len(smear['log10(E_nu/GeV)_max']))]
    dloge = [float(smear['log10(E_nu/GeV)_max'][i]) -float((smear['log10(E_nu/GeV)_min'][i])) for i in range(0,len(smear['log10(E_nu/GeV)_max']))]
    dphi = [float(smear['Dec_nu_max[deg]'][i]) - float((smear['Dec_nu_min[deg]'][i])) for i in range(0,len(smear['Dec_nu_max[deg]']))]
    
    e = [float(de[i])/float(dloge[i]) for i in range(0,len(smear['log10(E_nu/GeV)_max']))]
    
    #eff_a = [float(i) for i in smear['A_Eff[cm^2]']]
    #e2 = [pow((pow(10,float(smear['log10(E_nu/GeV)_max'][i])) + pow(10,float((smear['log10(E_nu/GeV)_min'][i]))))/2.0,2) for i in range(0,len(smear['log10(E_nu/GeV)_max']))]
    #e2dphide = [float(smear['Fractional_Counts'][i]) * e2[i] * dphi[i]/de[i] for i in range(0,len(smear['log10(E_nu/GeV)_max']))]
    

    yax = [e[i] * dphi[i]/dloge[i] for i in range(0, len(dloge))]

    #HISTOGRAMS OF N VS E BINNED IN LOG E

    plt.figure(figsize=(26,16))
    plt.rcParams.update({'font.size': 24})
    plt.hist(yax, bins= np.logspace(np.log10(1),12))
    plt.xscale("log")
    plt.xlabel("log(E_nu/GeV)")
    plt.ylabel("N")
    plt.suptitle(file.replace('icecube_10year_ps/irfs/','').replace('_smearing','').replace('.csv',''))
    plt.grid(True, which="major")
    #plt.show()
    plt.savefig(file.replace('icecube_10year_ps/irfs/','Task1b/histogram/smearing/').replace('.csv','.png'))

filenames = ["icecube_10year_ps/irfs/IC40_smearing.csv", "icecube_10year_ps/irfs/IC59_smearing.csv",
 "icecube_10year_ps/irfs/IC79_smearing.csv", "icecube_10year_ps/irfs/IC86_I_smearing.csv", "icecube_10year_ps/irfs/IC86_II_smearing.csv"]

for filename in filenames:
    n_e_hist(filename)