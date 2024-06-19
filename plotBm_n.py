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
ax2.set(xlim=(1,5), ylim=(1e0, 1e4))


Mpl = 2e27		# Planck mass [eV]
s2eV = (6.582119569e-16)#	// Hz to eV
J2eV = (1. / 1.602176634e-19)#	// Joules to eV (1 / e)
m2eV = (1.973269804e-7)#	// m-1 to eV
K2eV = (8.617333262e-5)#	// Kelvin to eV
kg2eV = 5.609588604e35#	// from hbar/c2
T2eV = 2e-16 * 1e18#		// Tesla to eV2 conversion [eV2/T]

dat = loadtxt("data/AI_n.csv",delimiter=',')
ax2.plot(dat[:,0],dat[:,1],color='blue', ls='--')
ax2.text(2,2e3,"AI")
plt.fill_between(dat[:,0],1e7,dat[:,1], facecolor='blue', alpha=0.5)

dat = loadtxt("data/TP_n.csv",delimiter=',')
ax2.plot(dat[:,0],dat[:,1],color='red', ls='--')
ax2.text(2,2e0,"TB")
plt.fill_between(dat[:,0],dat[:,1],0, facecolor='red', alpha=0.5)

dat = loadtxt("data/LFS_n.csv",delimiter=',')
ax2.plot(dat[:,0],dat[:,1],color='green', ls='--')
ax2.text(2,1e2,"LFS")
#dat2 = loadtxt("data/LFS_n-2.csv",delimiter=',')
#ax2.plot(dat2[:,0],dat2[:,1],color='red', ls='-')
#plt.fill_between(dat[:,0],dat2[:,1],dat[:,1], facecolor='red', alpha=0.5)


# axes
ax2.set_xlabel("n")
ax2.set_ylabel(r'$\beta_m$')
#ax2.set_xscale('log')
ax2.set_yscale('log')
#ax2.legend(loc='lower right')

plt.tight_layout()
name = "limits_n"
plt.savefig('plots/{}.jpg'.format(name))
plt.savefig('plots/pdfs/{}.pdf'.format(name))
plt.show()
