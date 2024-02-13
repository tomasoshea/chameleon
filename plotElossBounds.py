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
#ax2.set(xlim=(0,2e1), ylim=(0, 9e-35))
Mpl = 2e27		# Planck mass [eV]
Lsolar = 3.0128e28		# Solar luminosity [W]
s2eV = (6.582119569e-16)#	// Hz to eV
J2eV = (1. / 1.602176634e-19)#	// Joules to eV (1 / e)
Lsolar *= J2eV*s2eV


# Primakoff

dat = loadtxt("data/primakoff_Eloss_1e3.dat")
Bg = np.sqrt(Lsolar/10/dat[:,1])
ax2.plot(dat[:,0],Bg, ls='-', label=r'Primakoff, $\beta_m = 10^3$')

#dat = loadtxt("data/scalarB_Eloss_cham_1e3.dat")
#Bg = np.sqrt(Lsolar/10/dat[:,1])
#ax2.plot(dat[:,0],Bg, ls='--',label="n=1")

#dat = loadtxt("data/primakoff_Eloss_n10.dat")
#Bg = np.sqrt(Lsolar/10/dat[:,1])
#ax2.plot(dat[:,0],Bg, ls='--',label="n=10")

dat = loadtxt("data/scalarB_Eloss_cham_1e3.dat")
Bg = np.sqrt(Lsolar/10/dat[:,1])
ax2.plot(dat[:,0],Bg, ls='--',label=r'B-field, $\beta_m = 10^3$')

dat = loadtxt("data/scalarB_Eloss_cham_1e3--boost.dat")
Bg = np.sqrt(Lsolar/10/dat[:,1])
ax2.plot(dat[:,0],Bg, ls=':',label=r'B-field (boost), $\beta_m = 10^3$')



 
# axes
ax2.set_xlabel("Chameleon matter coupling")
ax2.set_ylabel("Chameleon photon coupling")	#[m-2 s-1 eV-1]")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/Eloss_B.jpg')
plt.show()
