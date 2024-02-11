import numpy as np
from scipy.constants import *
import time
import matplotlib.pyplot as plt

def coherence_operators(i, j, N):
    
    if i == 0:
        space = np.array([[0, 0], [1, 0]])
    elif j == 0:
        space = np.array([[0, 1], [0, 0]])
    else:
        space = np.identity(2)
    for k in range(1, N):
        if k == i:
            space = np.kron(space, np.array([[0, 0], [1, 0]]))
        elif k == j:
            space = np.kron(space, np.array([[0, 1], [0, 0]]))
        else:
            space = np.kron(space, np.identity(2))
        #print(k)

    return space

N = 10

def doSum(N):
    sum = np.zeros((2**N, 2**N), dtype=float)
    for i in range(N):
        for j in range(N):
            sum += coherence_operators(i, j, N)
    return sum

def getNonZeroEntries(mat):
    return len(np.where(mat == 1)[0])

x = range(1, N+1)
y = [getNonZeroEntries(doSum(k)) for k in x]

plt.plot(x, y, 'x')
plt.xscale('linear')
plt.yscale('linear')
plt.show()

"""
Unfortunately, the scaling for the number of non-zero entries in the "sum of coherence operators" is also exponential.

In conclusion, I need to figure out how the "block diagonalization" works.
"""
