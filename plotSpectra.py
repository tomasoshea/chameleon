# Tom O'Shea 2023

# plot scalar primakoff spectrum

from numpy import loadtxt
from matplotlib import pyplot as plt
import numpy as np


plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
#ax2 = fig2.add_axes((.15,.15,.8,.8))
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
ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-')
ax2.text(3e0,4e20,"T")

# B-field
dat = loadtxt("data/B_spectrum_1e2.dat")
dat[:,1] = dat[:,1]*1e3/s2eV	# convert [eV/eV] to [s-1 keV-1]
ax2.plot(dat[:,0]/1e3,dat[:,1],color='orange', ls='--')
ax2.text(2e-2,6e16,"B")
dat2 = loadtxt("data/B_spectrum_1e2--lowB.dat")
dat2[:,1] = dat2[:,1]*1e3/s2eV	# convert [eV/eV] to [s-1 keV-1]
ax2.plot(dat2[:,0]/1e3,dat2[:,1],color='orange', ls='--')
plt.fill_between(dat2[:,0]/1e3,dat2[:,1],dat[:,1], facecolor='orange', alpha=0.3)

## LL
##dat = loadtxt("data/LLspectrum_tot.dat")
#dat = loadtxt("data/coalescence_ll_spectrum_1e2--test.dat")
#dat[:,1] = dat[:,1]*1e3/s2eV/4/np.pi/np.pi	# convert [eV/eV] to [s-1 keV-1]
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-')

"""
# LL coalescence & L-primakoff 
datL = loadtxt("data/L_spectrum_1e2--test.dat")
wpmax = np.max(datL[:,0])

dat = loadtxt("data/coalescence_ll_spectrum_1e2--test.dat")
for i in range(len(dat[:,0])):
	if dat[i,0] < wpmax:
		dat[i,1] = np.nan
	else:
		dat[i,1] = dat[i,1]*1e3/s2eV	# convert [eV/eV] to [s-1 keV-1]
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-', color='r', label="L-L coalescence")
#plt.vlines(dat[-1,0]/1e3,1e-5,dat[-1,1],ls='-',color='r')

dat = datL
dat[:,1] = dat[:,1]*1e3/s2eV	# convert [eV/eV] to [s-1 keV-1]
for i in range(len(dat[:,0])):
	if dat[i,0] > 288: dat[i,0] = np.nan
ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-', color='r', label="L")
#plt.vlines(dat[-1,0]/1e3,1e-5,dat[-1,1],color='r')
ax2.text(2e-2,1e22,"L")
"""

#dat = loadtxt("data/L_spectrum_1e2.dat")
#dat[:,1] = dat[:,1]*1e3/s2eV	# convert [eV/eV] to [s-1 keV-1]
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-')
#ax2.text(2e-2,1e22,"L")

#dat = loadtxt("data/L_spectrum_1e2--omega_unbounded.dat")
#dat[:,1] = dat[:,1]*1e3/s2eV	# convert [eV/eV] to [s-1 keV-1]
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-.', label='L unbounded')

#dat = loadtxt("data/L_spectrum_1e2--omega_bounded1.dat")
#dat[:,1] = dat[:,1]*1e3/s2eV	# convert [eV/eV] to [s-1 keV-1]
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-.', label='L bounded (kD)')

#dat = loadtxt("data/L_spectrum_1e2--omega_kD.dat")
#dat[:,1] = dat[:,1]*1e3/s2eV	# convert [eV/eV] to [s-1 keV-1]
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-.', label='L bounded (kD)')


#datL = loadtxt("data/L_spectrum_1e2--omega_kD.dat")
#wpmax = np.max(datL[:,0])
#wpmax = 275
#dat = loadtxt("data/coalescence_ll_spectrum_1e2--test.dat")
#for i in range(len(dat[:,0])):
#	if dat[i,0] < wpmax:
#		dat[i,1] = np.nan
#	else:
#		dat[i,1] = dat[i,1]*1e3/s2eV	# convert [eV/eV] to [s-1 keV-1]
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', color='r', label="L-L coalescence")
#plt.vlines(dat[-1,0]/1e3,1e-5,dat[-1,1],ls=':',color='r')
#
#dat = datL
#dat[:,1] = dat[:,1]*1e3/s2eV	# convert [eV/eV] to [s-1 keV-1]
#for i in range(len(dat[:,0])):
#	if dat[i,0] > wpmax: dat[i,0] = np.nan
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', color='r', label="L")
#ax2.text(2e-2,4e18,"L")
#plt.vlines(dat[-1,0]/1e3,1e-5,dat[-1,1],color='r')

#dat = loadtxt("data/coalescence_lt_spectrum_1e2--uncapped.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', label="L-T coalescence")


# axes
ax2.set_xlabel("Scalar energy [keV]")
ax2.set_ylabel(r'$\beta_\gamma^{-2}\; \frac{\mathrm{d}\dot{N}}{\mathrm{d}\omega}$ [s$^{-1}$ keV$^{-1}$]')	#[m-2 s-1 eV-1]")
ax2.set_xscale('log')
ax2.set_yscale('log')
#ax2.legend(loc='lower right')

plt.tight_layout()
name = "spectrum_TB"
plt.savefig('plots/{}.jpg'.format(name))
plt.savefig('plots/pdfs/{}.pdf'.format(name))
plt.show()
