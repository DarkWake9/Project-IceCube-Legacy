import numpy as np
import random as rand

rounds = []
nsims = 10000000
p = 0.5
for i in range(1, nsims):
    r = 0
    nloss = 0
    while nloss != 2:
        r+=1
        if rand.random() < p:
            nloss = 0
        else: 
            nloss += 1
    rounds.append(r)

print(np.average(rounds))