# Tom O'Shea 2023

# plot scalar primakoff spectrum

from numpy import loadtxt
from matplotlib import pyplot as plt
import numpy as np


plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
#ax2.set(xlim=(0,1e1), ylim=(0, 1.01))
ax2.set(xlim=(1e-0,1e5), ylim=(0, 1.6e22))


dat = loadtxt("data/primakoff_spectrum_1e-3.dat")
#dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
#ax2.plot(dat[:,0],dat[:,1], ls='-', label='1 meV')

dat = loadtxt("data/primakoff_spectrum_1e0_newer.dat")
#dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
ax2.plot(dat[:,0],dat[:,1], ls='--', label='1 eV')

dat = loadtxt("data/primakoff_spectrum_1e3.dat")
#dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
#ax2.plot(dat[:,0],dat[:,1], ls=':', label='1 keV')

dat = loadtxt("data/primakoff_spectrum_1e4.dat")
#dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
#ax2.plot(dat[:,0],dat[:,1], ls='-.', label='10 keV')


# axes
ax2.set_xlabel("Scalar energy [eV]")
ax2.set_ylabel("Lambda^2 dPhi/dw [eV4]")	#[m-2 s-1 eV-1]")
ax2.set_xscale('log')
#ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/primakoff_spectrum_new-log.jpg')
plt.show()
