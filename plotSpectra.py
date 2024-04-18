# Tom O'Shea 2023

# plot scalar primakoff spectrum

from numpy import loadtxt
from matplotlib import pyplot as plt
import numpy as np


plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set(xlim=(1e-3,2e1), ylim=(1e-4, 1e6))
#ax2.set(xlim=(0,2e1), ylim=(0, 1.5e3))		# T linear limits


Mpl = 2e27		# Planck mass [eV]
s2eV = (6.582119569e-16)#	// Hz to eV
J2eV = (1. / 1.602176634e-19)#	// Joules to eV (1 / e)
m2eV = (1.973269804e-7)#	// m-1 to eV
K2eV = (8.617333262e-5)#	// Kelvin to eV
kg2eV = 5.609588604e35#	// from hbar/c2
T2eV = 2e-16 * 1e18#		// Tesla to eV2 conversion [eV2/T]


# Electron/Ion
dat = loadtxt("data/T_spectrum_1e2.dat")
ax2.plot(dat[:,0]/1e3,dat[:,1], ls='--', label="T")

dat = loadtxt("data/L_spectrum_1e2.dat")
ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-', color='r', label="L")
plt.vlines(dat[-1,0]/1e3,1e-5,dat[-1,1],color='r')

# B-field
dat = loadtxt("data/B_spectrum_1e2.dat")
ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-.', label='B')

# coalescence
dat = loadtxt("data/coalescence_ll_spectrum_1e2.dat")
ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', label="L-L coalescence")

dat = loadtxt("data/coalescence_lt_spectrum_1e2.dat")
ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', label="L-T coalescence")


# axes
ax2.set_xlabel("Scalar energy [keV]")
ax2.set_ylabel(r'$\beta_\gamma^{-2}\; \frac{\mathrm{d}N}{\mathrm{d}\omega}$')	#[m-2 s-1 eV-1]")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/spectrum_LTB_ll_lt--test.jpg')
plt.show()
