#!/usr/bin/python3

import sys
import numpy as np
import matplotlib.pyplot as plt
import time

#call as plotgrid.py filename nx ny dt [halfx] [sleep] [skip]

filename = sys.argv[1]
nx = int(sys.argv[2])
ny = int(sys.argv[3])
dt = float(sys.argv[4])

xlim = nx
if len(sys.argv) > 5:
    if int(sys.argv[5]) == 1:
        xlim = nx/2 + 1

if len(sys.argv) > 6:
    sleep = float(sys.argv[6])
else:
    sleep = 0.01

if len(sys.argv) > 7:
    skip = int(sys.argv[7])
else:
    skip = 1

raw = np.genfromtxt(filename)
sraw = raw.size
nt = raw.size//nx//ny

grid = np.zeros((nx,ny,nt))

n=0
for t in range(nt):
  for j in range(ny):
    for i in range(nx):
      grid[i,j,t]=raw[n]
      n=n+1


f = plt.figure()
a = f.gca()
f.show()

qmesh = a.pcolormesh(grid[:,:,0].T)
a.set_xlabel('x')
a.set_ylabel('y')
a.set_title('T = ' + str(0.0))
a.set_xlim(0,xlim)
a.set_ylim(0,ny)
cb = f.colorbar(qmesh, ax = a)
time.sleep(sleep)

for t in range(1, nt-1,skip):
  a.clear()
  a.set_xlabel('x')
  a.set_ylabel('y')
  a.set_title('T = ' + str(0.0))
  a.set_xlim(0,xlim)
  a.set_ylim(0,ny)
  qmesh = a.pcolormesh(grid[:,:,t].T)
  cb.update_bruteforce(qmesh)
  a.set_title('T = ' + str(dt * t))
  f.canvas.draw()
  time.sleep(sleep)
