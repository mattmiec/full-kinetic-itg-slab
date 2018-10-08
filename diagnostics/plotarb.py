#!/usr/bin/python3

import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage

#call as plotarb.py filename [sigma] [tmax]

filename = sys.argv[1]
if len(sys.argv)>2:
    sigma = int(sys.argv[2])
else:
    sigma = 0

data = np.genfromtxt(filename, skip_header = 1)

if len(sys.argv)>3:
    tmax = int(sys.argv[3])
else:
    tmax = data.shape[0]

with open(filename) as f:
    headers = f.readline().split()

t = data[:tmax,0]

for i in range(1, data.shape[1]):

    col = data[:tmax, i]

    if not all(val == 0.0 for val in col):

        if sigma != 0:
          col = scipy.ndimage.filters.gaussian_filter(col, sigma)

        plt.figure(0)
        plt.plot(t, col, label=headers[i])

        plt.figure(i)

        plt.suptitle(headers[i])

        plt.subplot(211)
        plt.plot(t, col)

        plt.subplot(212)
        plt.plot(t, abs(col))
        plt.yscale('log')
        plt.xlabel('t')

        plt.savefig(headers[i] + '.png')

plt.figure(0)
plt.title(filename[:-4])
plt.legend()
plt.savefig(filename[:-4] + '.png')
plt.show()
