# Tom O'Shea 2023

# plot scalar primakoff spectrum

from numpy import loadtxt
from matplotlib import pyplot as plt
import numpy as np

# constants
R = 6.9598E+8/(1.973269804e-7);	# solar radius [eV-1]


plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set(xlim=(0,1), ylim=(1e38, 1e48))
#ax2.set(xlim=(0,1), ylim=(0, 1.02e49))

dat = loadtxt("data/primakoff_profile_1e3.dat")
#dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
ax2.plot(dat[:,0]/R,dat[:,1], ls='-',label='Primakoff 1e3')

dat = loadtxt("data/scalarB_profile_cham-1e3.dat")
#dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
ax2.plot(dat[:,0]/R,dat[:,1], ls='-',label='B-field 1e3')

# axes
ax2.set_xlabel("Solar radius fraction")
ax2.set_ylabel("Emission rate dN/dr [eV Lambda2]")	#[m-2 s-1 eV-1]")
#ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/comparison_profile-log.jpg')
plt.show()
