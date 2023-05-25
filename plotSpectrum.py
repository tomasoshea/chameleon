# Tom O'Shea 2023

# plot prob against E for solar chameleons

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
#ax2.set(xlim=(0,1))#, ylim=(1e0, 1.2e3))
ax2.set(xlim=(0,3e3), ylim=(-0.01,1.01))

dat = loadtxt("data/cham-spectrum.dat")
dat[:,1] = dat[:,1] / np.nanmax(dat[:,1])
ax2.plot(dat[:,0],dat[:,1], color='magenta', label='Bm = 1e2, n=1')

dat = loadtxt("data/chamn4-spectrum.dat")
dat[:,1] = dat[:,1] / np.nanmax(dat[:,1])
ax2.plot(dat[:,0],dat[:,1], color='green', ls='--', label='Bm = 1e2, n=4')

dat = loadtxt("data/chamB6-spectrum.dat")
dat[:,1] = dat[:,1] / np.nanmax(dat[:,1])
ax2.plot(dat[:,0],dat[:,1], color='red', ls=':', label='Bm = 1e6, n=1')

dat = loadtxt("data/chamB5-spectrum.dat")
dat[:,1] = dat[:,1] / np.nanmax(dat[:,1])
ax2.plot(dat[:,0],dat[:,1], color='cyan', ls=':', label='Bm = 1e5, n=1')

# axes
ax2.set_xlabel("Chameleon energy [eV]")
ax2.set_ylabel("Tachocline chameleon flux [m-2 s-1 eV-1]")
#ax2.set_xscale('log')
#ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/cham-spectrum3.jpg')
plt.show()
