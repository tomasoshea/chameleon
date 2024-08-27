# Tom O'Shea 2024

# plot scalar primakoff profile

from numpy import loadtxt
from matplotlib import pyplot as plt
import numpy as np

# constants
R = 6.9598E+8/(1.973269804e-7);	# solar radius [eV-1]


plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.subplots()
ax2.set(xlim=(0,1), ylim=(5e-17, 1e-7))
#ax2.set(xlim=(0,1), ylim=(0, 1.02e49))

dat = loadtxt("data/T_profile_1e2.dat")
#dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
ax2.plot(dat[:,0]/R,dat[:,1], ls='--',label='T')

dat = loadtxt("data/B_profile_1e2.dat")
#dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
ax2.plot(dat[:,0]/R,dat[:,1], ls='-',label='B')

# axes
ax2.set_xlabel("Solar radius fraction")
ax2.set_ylabel(r'Emission rate $\beta_\gamma^{-2} \frac{\mathrm{d}N}{\mathrm{d}r}$ [eV]')
#ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.tight_layout()
name = "profile-log"
plt.savefig('plots/{}.jpg'.format(name))
plt.savefig('plots/pdfs/{}.pdf'.format(name))
plt.show()