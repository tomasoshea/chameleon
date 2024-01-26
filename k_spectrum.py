# Tom O'Shea 2023

# plot scalar primakoff spectrum

from numpy import loadtxt
from matplotlib import pyplot as plt
import numpy as np


plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
#ax2.set(xlim=(1e-3,20), ylim=(1e12, 1e22))
#ax2.set(xlim=(2e2,2e4), ylim=(0, 1.01))


# Primakoff

#dat = loadtxt("data/L_primakoff_k-spectrum_fixed_1e-3--test.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-', label='m = 1 meV (L)')
#
#dat = loadtxt("data/L_primakoff_k-spectrum--test2.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', label='m = 1 eV (L)')
#
#dat = loadtxt("data/L_primakoff_diff_k-spectrum_fixed_1e-3.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', label='m = 1 meV, w = 20 keV')
#

dat = loadtxt("data/L_primakoff_diff_k-spectrum_fixed_1e-3.dat")
ax2.plot(dat[:,0],dat[:,1], ls='-', label='m = 1 meV, w = 20 keV')

dat = loadtxt("data/L_primakoff_diff_k-spectrum_fixed_1e-3_1e3.dat")
ax2.plot(dat[:,0],dat[:,1], ls='-', label='m = 1 meV, w = 1 keV')

dat = loadtxt("data/L_primakoff_diff_k-spectrum_fixed_1e-3_1e2.dat")
ax2.plot(dat[:,0],dat[:,1], ls='-', label='m = 1 meV, w = 0.1 keV')

#dat = loadtxt("data/L_primakoff_diff_k-spectrum_fixed_1e0.dat")
#ax2.plot(dat[:,0],dat[:,1], ls='-.', label='m = 1 eV, w = 20 keV')
#
#dat = loadtxt("data/L_primakoff_diff_k-spectrum_fixed_1e3.dat")
#ax2.plot(dat[:,0],dat[:,1], ls='--', label='m = 1 keV, w = 20 keV')
#
#dat = loadtxt("data/L_primakoff_diff_k-spectrum_fixed_1e3.dat")
#ax2.plot(dat[:,0],dat[:,1], ls=':', label='m = 19 keV, w = 20 keV')

#dat = loadtxt("data/L_primakoff_diff_k-spectrum_fixed_1e-3.dat")
#dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
#ax2.plot(np.sqrt((20e3)**2 - (dat[:,0])**2),dat[:,1], ls='-', label='m = 1 meV, w = 20 keV')
#
#dat = loadtxt("data/L_primakoff_diff_k-spectrum_fixed_1e0.dat")
#dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
#ax2.plot(np.sqrt((20e3)**2 - (dat[:,0])**2),dat[:,1], ls='--', label='m = 1 eV, w = 20 keV')
#
#dat = loadtxt("data/L_primakoff_diff_k-spectrum_fixed_1e3.dat")
#dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
#ax2.plot(np.sqrt((20e3)**2 - (dat[:,0])**2),dat[:,1], ls=':', label='m = 1 keV, w = 20 keV')



# axes
ax2.set_xlabel("k_l [eV]")
ax2.set_ylabel("dPhi/dk [eV2] * Lambda2")	#[m-2 s-1 eV-1]")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/k-diff_k-spectrum_fixed_(w-k)-log--test.jpg')
plt.show()
