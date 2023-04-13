import numpy as np
import time

arr=np.random.rand(1000000)

sum=0
start=time.time()
for j in range(100):
    for i in arr:
        sum+=i
end=time.time()
print(end-start)

sum2=0
start=time.time()
for j in range(100):
    sum2+=np.sum(arr)+j
end=time.time()
print(end-start)

