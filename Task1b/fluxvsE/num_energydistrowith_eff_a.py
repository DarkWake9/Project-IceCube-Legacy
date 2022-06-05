
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
with open("icecube_10year_ps/irfs/IC59_effectiveArea.csv", 'r') as f2:
    lines = f2.readlines()
content=[]
column=lines[0].split()
column.pop(0)
for line in lines[1:]:
    
    content.append(line.split())
data = pd.DataFrame(content, columns = column)
#data = pd.DataFrame(content).drop(0)

#ASSUMING number of particles = flux * effective area

de = [pow(10,float(data['log10(E_nu/GeV)_max'][i])) - pow(10,float((data['log10(E_nu/GeV)_min'][i]))) for i in range(0,len(data['log10(E_nu/GeV)_max']))]
dphi = [float(data['Dec_nu_max[deg]'][i]) - float((data['Dec_nu_min[deg]'][i])) for i in range(0,len(data['Dec_nu_max[deg]']))]
#print(de)
#print(de)
eff_a = [float(data['A_Eff[cm^2]'][i]) for i in range(0,len(data['A_Eff[cm^2]']))]
e = [((pow(10,float(data['log10(E_nu/GeV)_max'][i])) + pow(10,float((data['log10(E_nu/GeV)_min'][i]))))/2.0) for i in range(0,len(data['log10(E_nu/GeV)_max']))]
e2 = [pow((pow(10,float(data['log10(E_nu/GeV)_max'][i])) + pow(10,float((data['log10(E_nu/GeV)_min'][i]))))/2.0,2) for i in range(0,len(data['log10(E_nu/GeV)_max']))]
#flux = [e2[i] * dphi[i]/de[i] for i in range(0,len(data['log10(E_nu/GeV)_max']))]
n_e = [eff_a[i]*(e2[i] * dphi[i]/de[i]) for i in range(0,len(data['log10(E_nu/GeV)_max']))]


plt.figure(figsize=(32,18))
plt.rcParams.update({'font.size': 48})
plt.xlabel("log(E_nu/GeV)")
plt.ylabel("n(E) = A_eff * E2dphi/dE")
plt.title("Energy distribution diagram")
plt.plot(np.log10(e),n_e)
#plt.suptitle(filenamewithpath.replace('icecube_10year_ps/irfs/',''))
plt.grid(True, which='both')
#plt.show()
plt.savefig('test.jpg')




