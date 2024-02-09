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

x = np.linspace(1, N, 2)
y = [getNonZeroEntries(doSum(k)) for k in x]

plt.plot(x, y)
plt.xscale('linear')
plt.yscale('linear')
plt.show()
for i in range(1,N):
    print("dimension: " + str(2**i * 2**i))
    print("the non-zero entries: " + str(len(np.where(doSum(i) == 1)[0])))

print((y[1]-y[0])/(x[1]-x[0]))
print((y[2]-y[1])/(x[2]-x[1]))
print((y[3]-y[2])/(x[3]-x[2]))
print((y[4]-y[3])/(x[4]-x[3]))