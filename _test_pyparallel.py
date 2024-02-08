import numpy as np
from scipy.constants import *
import time
from multiprocessing import Pool

a = np.array([[1, 2], [3, 4]])
b = np.array([[5, 5], [5, 5]])

pool = Pool()

def matmul(a, b, N):
    for i in range(N):
        return a @ b

N = 2**20

t1 = time.time()

results = pool.apply_async(matmul, args=(a, b, N))
answer = results.get(timeout=10)

t2 = time.time()
print(t2 - t1)