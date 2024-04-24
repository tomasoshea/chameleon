# Tom O'Shea 2023

# plot scalar primakoff spectrum

from numpy import loadtxt
from matplotlib import pyplot as plt
import numpy as np


plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set(xlim=(1e0,1e99), ylim=(1e9, 1e13))
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


dat = loadtxt("data/Eloss_Bg--1.dat")
Bg = np.sqrt(3*Lsolar/100/dat[:,1])
ax2.plot(dat[:,0],Bg, ls='-', label=r'$n=1$')

# lower limit
plt.vlines(1e10,1e0,1e20,color='r')
# upper limit from theory breakdown
plt.vlines(1e40,1e0,1e20,color='r')


 
# axes
ax2.set_xlabel(r'$\beta_\gamma$')
#ax2.set_ylabel(r'Energy loss rate [eV2 $\beta_\gamma^{-2}$]')
ax2.set_ylabel(r'$\beta_\gamma$ bound')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/Eloss_Bg.jpg')
plt.show()
