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
ax2.set(xlim=(1e-3,2e1), ylim=(1e-42, 1e-33))
Mpl = 2e27		# Planck mass [eV]


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

dat = loadtxt("data/primakoffV3_spectrum_cham_1e3.dat")
ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-', label="T")

dat = loadtxt("data/primakoffV3_L-spectrum_cham_1e3.dat")
ax2.plot(dat[:,0]/1e3,dat[:,1], ls='-', color='orange', label="L")
plt.vlines(dat[-1,0]/1e3,1e-99,dat[-1,1],color='orange')
print(dat[-1,0]/1e3)

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
ax2.set_ylabel(r'$\beta_\gamma^{-2}\; \frac{\mathrm{d}\Phi}{\mathrm{d}\omega}$ [eV$^2$]')	#[m-2 s-1 eV-1]")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/primakoffV3-LT--log.jpg')
plt.show()
