#FINAL VERSION ON 06-06-2022
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def n_e_histf(file):
    with open(file, 'r') as f2:
        lines = f2.readlines()
    content=[]
    column=lines[0].split()
    column.pop(0)
    for line in lines[1:]:            
        
        content.append(line.split())

    data = pd.DataFrame(content, columns = column)
    f2.close()

    e = [float(data['log10(E/GeV)'][i]) for i in range(0,len(data['log10(E/GeV)']))]

    #HISTOGRAMS OF N VS E BINNED IN LOG E
    plt.figure(figsize=(28,16))
    plt.rcParams.update({'font.size': 24})
    plt.hist(e, bins= np.logspace(np.log10(2.6),np.log10(7)))
    plt.xscale("log")
    plt.xlabel("E_nu (GeV)")
    plt.ylabel("N")
    plt.suptitle(file.replace('icecube_10year_ps/events/','').replace('_exp','').replace('.csv',''))
    plt.grid(True, which="both")
    #plt.show()
    plt.savefig(file.replace('icecube_10year_ps/events/','Task1b/output/').replace('_exp','').replace('.csv','.png'))
    

filenames = ["icecube_10year_ps/events/IC40_exp.csv", "icecube_10year_ps/events/IC59_exp.csv",
 "icecube_10year_ps/events/IC79_exp.csv", "icecube_10year_ps/events/IC86_I_exp.csv", "icecube_10year_ps/events/IC86_II_exp.csv",
  "icecube_10year_ps/events/IC86_III_exp.csv", "icecube_10year_ps/events/IC86_IV_exp.csv", "icecube_10year_ps/events/IC86_V_exp.csv",
   "icecube_10year_ps/events/IC86_VI_exp.csv","icecube_10year_ps/events/IC86_VII_exp.csv"]

for filename in filenames:
    n_e_histf(filename)