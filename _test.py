import numpy as np
from scipy.constants import *
import time

a = np.array([[1, 2], [3, 4]])
b = np.array([[5, 5], [5, 5]])

c = np.kron(a, b)
print(c)

t1 = time.time()

for i in range(2**20): #1 million (for e.g. 20 atoms)
    c = a * b

t2 = time.time()
print(t2 - t1)

"""
matmul 1 million times in
Python: 
Time on HPC ~0.75 seconds
Time on desktop ~0.53 seconds ...
C++:
Time on HPC ~
Time on desktop ~0.04 seconds
"""