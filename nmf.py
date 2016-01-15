# MULTIPLICATIVE UPDATES FOR NMF ---------------------------
# nmf.py
#
# Thu Jan 14 2016
# From Lee and Seung, "Algorithms for NMF", NIPS 2001
#

import random
import scipy as sp
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

### Helper function to format matrices as tables
def printmat(mat):
    print(pd.DataFrame(mat))

### Set up initial conditions
pp = 100
nn = 100
rr = 50
W0 = np.random.exponential(10, [pp, rr])
H0 = np.random.exponential(10, [rr, nn])
V = W0@H0

### Generate random initial guesses (exponential => positive entries)
W = np.random.exponential(10, [pp, rr])
H = np.random.exponential(10, [rr, nn])

### Data frame for storing errors after each iteration
res = np.full([iters, 3], np.nan)
res = pd.DataFrame(res)
res.columns = ["L2", "L1", "max"]
# printmat(res)

### Run the main algorithm
iters = 1000
updated = W @ H
for i in range(iters):
    L2err = np.sqrt(np.sum(np.square(V - updated)))
    L1err = np.sum(np.absolute(V - updated))
    maxerr = np.max(np.absolute(V-updated))
    if maxerr > L1err:
        print("Whoa!")

    res.iloc[i] = [L2err, L1err, maxerr]
    H = H * (W.T @ V) / (W.T @ W @ H)       # Update for H
    W = W * (V @ H.T) / (W @ H @ H.T)       # Update for W
    updated = W @ H

# printmat(V)
# printmat(updated)
# printmat(V - updated)
printmat(res)

### Plot the results
plt.plot(res['L2'])
plt.yscale('log')
plt.show()
