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
ax2.set(xlim=(1e-6,1e-2), ylim=(1e10, 1e11))
Mpl = 2e27		# Planck mass [eV]
Lsolar = 3.0128e28		# Solar luminosity [W]
s2eV = (6.582119569e-16)#	// Hz to eV
J2eV = (1. / 1.602176634e-19)#	// Joules to eV (1 / e)
Lsolar *= J2eV*s2eV


# Primakoff

#dat = loadtxt("data/primakoff_Eloss_1e3.dat")
#Bg = np.sqrt(Lsolar/10/dat[:,1])
#ax2.plot(dat[:,0],Bg, ls='-', label=r'Primakoff, $\beta_m = 10^3$')

#dat = loadtxt("data/scalarB_Eloss_cham_1e3.dat")
#Bg = np.sqrt(Lsolar/10/dat[:,1])
#ax2.plot(dat[:,0],Bg, ls='--',label="n=1")

#dat = loadtxt("data/primakoff_Eloss_n10.dat")
#Bg = np.sqrt(Lsolar/10/dat[:,1])
#ax2.plot(dat[:,0],Bg, ls='--',label="n=10")

#dat = loadtxt("data/scalarB_Eloss_cham_1e3.dat")
#Bg = np.sqrt(Lsolar/10/dat[:,1])
#ax2.plot(dat[:,0],Bg, ls='--',label=r'B-field, $\beta_m = 10^3$')

dat = loadtxt("data/primakoff_total_Eloss_Lambda_10e2.dat")
Bg = np.sqrt(3*Lsolar/100/dat[:,1])
ax2.plot(dat[:,0],Bg, ls='-', label=r'$n=1$, $\beta_m=10^2$')

dat = loadtxt("data/primakoff_total_Eloss_Lambda_10e4.dat")
Bg = np.sqrt(3*Lsolar/100/dat[:,1])
ax2.plot(dat[:,0],Bg, ls='-', label=r'$n=1$, $\beta_m=10^4$')

dat = loadtxt("data/primakoff_total_Eloss_Lambda_10e6.dat")
Bg = np.sqrt(3*Lsolar/100/dat[:,1])
ax2.plot(dat[:,0],Bg, ls='-', label=r'$n=1$, $\beta_m=10^6$')

 
# axes
ax2.set_xlabel(r'$\Lambda$ (eV)')
ax2.set_ylabel("Chameleon photon coupling")	#[m-2 s-1 eV-1]")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/Eloss_total_Lambda.jpg')
plt.show()
