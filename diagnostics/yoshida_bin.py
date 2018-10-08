#!/usr/bin/python3

import sys
import numpy as np
import matplotlib.pyplot as plt

#call yoshida_bin.py filename mode [imag] [minfreq] [maxfreq] [plot]

def main():
    filename = sys.argv[1]
    mode = int(sys.argv[2])

    if len(sys.argv) > 3:
      imag = int(sys.argv[3])
    else:
      imag = 0

    data = np.genfromtxt(filename, skip_header=1)

    t = data[:, 0]
    re = data[:, 2 * mode + 1]
    im = data[:, 2 * mode + 2]

    dt = t[1] - t[0]

    if len(sys.argv) > 4:
      minfreq = float(sys.argv[4])
    else:
      minfreq = 0.

    if len(sys.argv) > 5:
      maxfreq = float(sys.argv[5])
    else:
      maxfreq = np.pi/dt

    if len(sys.argv) > 6:
      plot = int(sys.argv[6])
    else:
      plot = 1


    if imag == 1:
      (w, d, a, p) = fyoshida_bin(im, dt, minfreq, maxfreq, plot)
    else:
      (w, d, a, p) = fyoshida_bin(re, dt, minfreq, maxfreq, plot)

    print('frequency is ', w)
    print('growth rate is ', d)


def fyoshida_bin(x, dt, minf, maxf, plot = False):
#return [W, D, A, p]

  N = x.size
  nn = np.array(range(N))
  n2 = np.array(range(N//2))
  Xw = np.fft.fft(x)

  df = 2.*np.pi/N/dt
  n2bin = n2[int(minf/df): int(maxf/df+1)]

  Xabs = np.amax(np.absolute(Xw[n2bin]))
  ind = np.argmax(np.absolute(Xw[n2bin])) + int(minf/df)

  #plot frequency spectrum
  if (plot == True):
    plt.plot(2*np.pi*n2/dt/N,np.absolute(Xw[n2])**2)
    plt.yscale('log')
    plt.xlabel('$\omega$')
    plt.savefig('freq.png')
    plt.show()

  #set four indices to be used
  k = np.array([ind-1, ind, ind+1])
  if (ind-2)>0:
    if abs(Xw[ind+2]) > abs(Xw[ind-2]):
      k = np.append(k,ind+2)
    else:
      k = np.insert(k,0,ind-2)
  else:
    k = np.append(k,ind+2)
  R = (Xw[k[0]]-2*Xw[k[1]]+Xw[k[2]])/(Xw[k[1]]-2*Xw[k[2]]+Xw[k[3]])
  #growth rate
  D = (2*np.pi)/N*(-3/(R-1)-1).imag
  #frequency
  if (k[3]-ind == 1):
    W = (2*np.pi/N)*(ind-3/(R-1)-2).real
  else:
    W = (2*np.pi/N)*(ind-3/(R-1)-1).real
  exarg = -1j*W*nn
  wsum = np.sum(np.exp(-1*D*nn))
  v = x*np.exp(exarg)
  A = 2*abs(np.sum(v)/wsum)
  p = -np.arctan(np.sum(v).imag/np.sum(v).real)

  W = W/dt
  D = -D/dt

  return (W,D,A,p)


if __name__ == "__main__":
    main()
