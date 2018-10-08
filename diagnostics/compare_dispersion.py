import sys
sys.path.insert(0, '/home/mtmiec/fk-itg-slab/diagnostics/')
import yoshida_bin
import numpy as np
import matplotlib.pyplot as plt
import dielectric
from scipy.optimize import fsolve

# GET FREQUENCIES FROM PHI_HIST

omega_array_full = []
itg_array_full = []
filename = '63' + '/phi_hist.out'
data = np.genfromtxt(filename, skip_header=1)

for mode in range(0,5):
  t = data[:, 0]
  re = data[:, 2 * mode + 1]
  im = data[:, 2 * mode + 2]

  dt = t[1] - t[0]

  my_omega_array = []
  for j in range(4):
    minfreq = float(j) + 0.05
    maxfreq = float(j+1)

    (w, d, a, p) = yoshida_bin.fyoshida_bin(im, dt, minfreq, maxfreq, False)
    my_omega_array.append(w)

  omega_array_full.append(my_omega_array)

omega_array_full = np.array(omega_array_full)
print(omega_array_full)

# CALCULATE FREQUENCIES FROM DISPERSION RELATION

th = 1.
kx = np.linspace(0.1,0.5,5)
alpha = 0.0

guess = np.array([[1.99,0.0],[2.999,0.0],[3.99999,0.0]])
om = np.zeros((guess.shape[0],len(kx)))
gm = np.zeros((guess.shape[0],len(kx)))

for i in range(guess.shape[0]):
    for j in range(len(kx)):
        om[i,j], gm[i,j] = fsolve(dielectric.dielectricFunc_kpar0_ky0, guess[i,:], args = (th, kx[j], alpha))
        guess[i,:] = np.array([om[i,j],gm[i,j]])

# PLOT

plt.figure()
plt.title("Drift-kinetic-electron Theory vs Drift-fluid-electron Simulation\n($k_\parallel = k_y = 0.0$, $\kappa_T = \kappa_n = 0.0$)", fontsize=14)

lines = ['r-','b-','g-','y-']
markers = ['ro', 'bo', 'go', 'yo']

for i in range(4):
  plt.plot(np.linspace(.1,.5,5), omega_array_full[:,i], markers[i])

for i in range(0, guess.shape[0]):
  plt.plot(kx,om[i,:],lines[i+1])

plt.ylim([0.0,4.5])
plt.xlim([0.0,0.6])
plt.xlabel("$k_x$")
plt.ylabel("Re[$\omega$]")
plt.legend(["EIC?","PIB 1","PIB 2","PIB 3"],loc="upper right")

plt.gcf().subplots_adjust(top=.8, left=.15, wspace = 0.4)
plt.savefig("disp.png")
plt.show()
