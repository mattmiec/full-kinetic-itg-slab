#!/usr/bin/python3

import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage

#call as plotmodes.py filename [sigma] [tmax]

filename = sys.argv[1]

data = np.genfromtxt(filename, skip_header = 1)

if len(sys.argv)>2:
    sigma = int(sys.argv[2])
else:
    sigma = 0

if len(sys.argv)>3:
    tmax = int(sys.argv[3])
else:
    tmax = data.shape[0]

nmodes = (data.shape[1]-1)//2

with open(filename) as f:
    headers = f.readline().split()

t = data[:tmax,0]

for i in range(nmodes):

    re = data[:tmax, 2 * i + 1]
    im = data[:tmax, 2 * i + 2]

    if not (all(val == 0.0 for val in re) and all(val == 0.0 for val in im)):

        if sigma != 0:
          re = scipy.ndimage.filters.gaussian_filter(re, sigma)
          im = scipy.ndimage.filters.gaussian_filter(im, sigma)

        am = np.sqrt(re**2+im**2)

        plt.figure(i)

        plt.subplot(211)
        plt.plot(t,re,t,im,t,am)
        plt.legend([headers[2*i + 1], headers[2*i + 2],'(Re + Im)^2'],loc=0) #best

        plt.subplot(212)
        plt.plot(t,am)
        plt.yscale('log')
        plt.xlabel('t')
        plt.legend(['(Re + Im)^2'])

        plt.savefig('mode'+str(i)+'.png')

plt.show()
