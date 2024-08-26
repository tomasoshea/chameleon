# Tom O'Shea 2023

# plot Debye screening scale

from numpy import loadtxt
from matplotlib import pyplot as plt
import numpy as np


plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
#ax2 = fig2.add_axes((.15,.15,.8,.8))
ax2 = fig2.subplots()
ax2.set(xlim=(0,1), ylim=(0, 10))


Mpl = 2e27		# Planck mass [eV]
s2eV = (6.582119569e-16)#	// Hz to eV
J2eV = (1. / 1.602176634e-19)#	// Joules to eV (1 / e)
m2eV = (1.973269804e-7)#	// m-1 to eV
K2eV = (8.617333262e-5)#	// Kelvin to eV
kg2eV = 5.609588604e35#	// from hbar/c2
T2eV = 2e-16 * 1e18#		// Tesla to eV2 conversion [eV2/T]
rSolar = 6.957e8/m2eV#				// solar radius [eV-1]



# Electron/Ion T
r = loadtxt("data/r.dat")
kappa = loadtxt("data/kappa.dat")
ax2.plot(r/rSolar, kappa/1e3, ls='-', label="kappa")



# axes
ax2.set_xlabel(r'$r/R_\odot$')
ax2.set_ylabel(r'$\kappa$ [keV]')	#[m-2 s-1 eV-1]")
#ax2.set_xscale('log')
#ax2.set_yscale('log')
#ax2.legend()

plt.tight_layout()
plt.show()
