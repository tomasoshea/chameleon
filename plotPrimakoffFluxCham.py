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
ax2.set(xlim=(1e-2,2e1), ylim=(1e14, 5e19))


# Primakoff

#dat = loadtxt("data/primakoff_spectrum_1e-3_newer.dat")
##dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
##ax2.plot(dat[:,0]/1e3,dat[:,1], ls='--', label='m = 1 meV')
#
#dat = loadtxt("data/primakoff_spectrum_1e0_newer.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-', label='m = 1 eV (T)')
#
#dat = loadtxt("data/L_primakoff_spectrum_fixed_1e0_newer.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', label='m = 1 eV (L)')


#dat = loadtxt("data/primakoff_spectrum_cham_1e3_newer.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='--', label='Primakoff 1e3')
#
#dat = loadtxt("data/primakoff_spectrum_cham_1e6_newer.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='--', label='Primakoff 1e6')
#
#dat = loadtxt("data/primakoff_spectrum_cham_1e7_newer.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='--', label='Primakoff 1e7')


# B-field

#dat = loadtxt("data/scalarB_spectrum_cham_1e3.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-.', label='B-field 1e3 (test)')

dat = loadtxt("data/scalarB_spectrum_fixed_1e-3-test.dat")
ax2.plot(dat[:,0]/1e3,dat[:,1], ls='--', label='m = 1 meV')

dat = loadtxt("data/scalarB_spectrum_fixed_1e0-test.dat")
ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-.', label='m = 1 eV')

dat = loadtxt("data/scalarB_spectrum_fixed_1e3.dat")
ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', label='m = 1 keV')

#dat = loadtxt("data/scalarB_spectrum_cham_1e6.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', label='B-field 1e6')
#
#dat = loadtxt("data/scalarB_spectrum_cham_1e7.dat")
#ax2.plot(dat[:,0]/1e3,dat[:,1], ls=':', label='B-field 1e7')
 
# axes
ax2.set_xlabel("Scalar energy [keV]")
ax2.set_ylabel("dPhi/dw [eV2] * Lambda2")	#[m-2 s-1 eV-1]")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/scalarB_fixed-log.jpg')
plt.show()
