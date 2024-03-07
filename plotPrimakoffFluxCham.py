# Tom O'Shea 2023

# plot scalar primakoff spectrum

from numpy import loadtxt
from matplotlib import pyplot as plt
import numpy as np


plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set(xlim=(1e-3,2e1), ylim=(1e-6, 5e3))

Mpl = 2e27		# Planck mass [eV]
s2eV = (6.582119569e-16)#	// Hz to eV
J2eV = (1. / 1.602176634e-19)#	// Joules to eV (1 / e)
m2eV = (1.973269804e-7)#	// m-1 to eV
K2eV = (8.617333262e-5)#	// Kelvin to eV
kg2eV = 5.609588604e35#	// from hbar/c2
T2eV = 2e-16 * 1e18#		// Tesla to eV2 conversion [eV2/T]


# Primakoff

#dat = loadtxt("data/primakoff_spectrum_1e-3_newer.dat")
##dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
##ax2.plot(dat[:,0]/1e3,dat[:,1], ls='--', label='m = 1 meV')
#
#dat = loadtxt("data/primakoff_spectrum_1e0_newer.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-', label='m = 1 eV (T)')
#
#dat = loadtxt("data/primakoffV3_spectrum_fixed_1e0.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1]/(16*Mpl**2), ls='-', label='T-plasmon')

dat = loadtxt("data/primakoffV3_spectrum_1e2.dat")
ax2.plot(dat[:,0]/1e3,dat[:,1], ls='--', label="T")

dat = loadtxt("data/primakoffV3_L_spectrum_1e2.dat")
ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-', color='r', label="L")
plt.vlines(dat[-1,0]/1e3,1e-99,dat[-1,1],color='r')

# total
#dat = loadtxt("data/primakoffV3_total-spectrum_cham_1e3.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-.',lw='2',color='k', label="combined")

# B-field
dat = loadtxt("data/scalarB_spectrum_1e2.dat")
ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-.', label='B')

# coalescence
dat = loadtxt("data/coalescence_ll_spectrum_1e2.dat")
ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', label="L-L coalescence")

#dat = loadtxt("data/coalescence_lt_spectrum_1e4.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', label="L-T coalescence")


# luca-linda
#dat = loadtxt("luca-linda/data/plotres_1e3.dat", delimiter=',')	# w [keV], dPhi/dw [cm-2 s-1 keV-1]
#dat[:,1] = dat[:,1]*1e4*m2eV*m2eV*s2eV/1e3
#ax2.plot(dat[:,0],dat[:,1], ls=':', label="Luca/Linda 1")

# luca-linda
#dat = loadtxt("luca-linda/data/plotnonres_1e3.dat", delimiter=',')	# w [keV], dPhi/dw [cm-2 s-1 keV-1]
#dat[:,1] = dat[:,1]*1e4*m2eV*m2eV*s2eV/1e3
#ax2.plot(dat[:,0],dat[:,1], ls=':', label="Luca/Linda 2")

#dat = loadtxt("data/primakoffV3_spectrum_n4_1e3.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-.', label="n=4")

#dat = loadtxt("data/primakoffV3_spectrum_n10_1e3.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='--', label="n=10")

#dat = loadtxt("data/primakoffV3_spectrum_n400_1e3.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', label="n=400")

#dat = loadtxt("data/primakoffV3_spectrum_cham_1e6.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='--', label=r'$\beta_m = 10^6$')

#dat = loadtxt("data/primakoffV3_spectrum_cham_1e7.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', label=r'$\beta_m = 10^7$')

#dat = loadtxt("data/primakoffV3_L-spectrum_fixed_1e0.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-', label='L-plasmon')


# B-field

#dat = loadtxt("data/scalarB_spectrum_cham_1e3.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-.', label='B-field 1e3 (test)')

#dat = loadtxt("data/scalarB_spectrum_fixed_1e-3-test.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='--', label='m = 1 meV')
#
#dat = loadtxt("data/scalarB_spectrum_fixed_1e0-test.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-.', label='m = 1 eV')
#
#dat = loadtxt("data/scalarB_spectrum_fixed_1e3.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', label='m = 1 keV')

#dat = loadtxt("data/scalarB_spectrum_cham_1e6.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', label='B-field 1e6')
#
#dat = loadtxt("data/scalarB_spectrum_cham_1e7.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', label='B-field 1e7')

#dat = loadtxt("data/primakoffV3_spectrum_fixed_1e-3--test.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-', label='Raffelt screened, full')

#dat = loadtxt("data/primakoffV3_spectrum_fixed_1e-3--test2.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-.', label='Raffelt screened, reduced')

#dat = loadtxt("data/primakoffV3_spectrum_fixed_1e-3--test5.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', label='Pi screened (B), full')

#dat = loadtxt("data/primakoffV3_spectrum_fixed_1e-3--test4.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='--', label='Pi screened, full')

#dat = loadtxt("data/primakoffV3_spectrum_fixed_1e-3--test3.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', label='Pi screened, reduced')

 
# axes
ax2.set_xlabel("Scalar energy [keV]")
ax2.set_ylabel(r'$\beta_\gamma^{-2}\; \frac{\mathrm{d}N}{\mathrm{d}\omega}$')	#[m-2 s-1 eV-1]")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/primakoffV3-coalescence--log.jpg')
plt.show()
