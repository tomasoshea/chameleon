# Tom O'Shea 2024

# plot scalar primakoff spectrum

from numpy import loadtxt
from matplotlib import pyplot as plt
import numpy as np


plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.subplots()
#ax2.set(xlim=(1e-3,2e1), ylim=(1e-4, 9e5))
ax2.set(xlim=(1e-2,2e1), ylim=(1e14, 1e22))
#ax2.set(xlim=(0,2e1), ylim=(0, 1.5e3))		# T linear limits


Mpl = 2e27		# Planck mass [eV]
s2eV = (6.582119569e-16)#	// Hz to eV
J2eV = (1. / 1.602176634e-19)#	// Joules to eV (1 / e)
m2eV = (1.973269804e-7)#	// m-1 to eV
K2eV = (8.617333262e-5)#	// Kelvin to eV
kg2eV = 5.609588604e35#	// from hbar/c2
T2eV = 2e-16 * 1e18#		// Tesla to eV2 conversion [eV2/T]


# Electron/Ion T
dat = loadtxt("data/T_spectrum_1e2.dat")
dat[:,1] = dat[:,1]*1e3/s2eV	# convert [eV/eV] to [s-1 keV-1]
ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-',color='k')
ax2.text(3e0,4e20,"T")

# B-field
dat = loadtxt("data/B_spectrum_1e2.dat")
dat[:,1] = dat[:,1]*1e3/s2eV	# convert [eV/eV] to [s-1 keV-1]
ax2.plot(dat[:,0]/1e3,dat[:,1],color='k', ls='--')
ax2.text(2e-2,6e16,"B")
dat2 = loadtxt("data/B_spectrum_1e2--lowB.dat")
dat2[:,1] = dat2[:,1]*1e3/s2eV	# convert [eV/eV] to [s-1 keV-1]
ax2.plot(dat2[:,0]/1e3,dat2[:,1],color='k', ls='--')
plt.fill_between(dat2[:,0]/1e3,dat2[:,1],dat[:,1], facecolor='k', alpha=0.2)#, hatch='xx'

# axes
ax2.set_xlabel("Scalar energy [keV]")
ax2.set_ylabel(r'$\beta_\gamma^{-2}\; \frac{\mathrm{d}\dot{N}}{\mathrm{d}\omega}$ [s$^{-1}$ keV$^{-1}$]')	#[m-2 s-1 eV-1]")
ax2.set_xscale('log')
ax2.set_yscale('log')
#ax2.legend(loc='lower right')

plt.tight_layout()
name = "spectrum_TB--black"
plt.savefig('plots/{}.jpg'.format(name))
plt.savefig('plots/pdfs/{}.pdf'.format(name))
plt.show()
