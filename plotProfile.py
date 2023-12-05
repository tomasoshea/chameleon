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
#ax2.set(xlim=(0,1e1), ylim=(0, 1.01))
ax2.set(xlim=(0,1), ylim=(0, 1.02e49))

dat = loadtxt("data/primakoff_profile_1eV.dat")
#dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
ax2.plot(dat[:,0]/R,dat[:,1], ls='-',label='Fixed m = 1eV')

dat = loadtxt("data/primakoff_profile_cham_1e6.dat")
#dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
ax2.plot(dat[:,0]/R,dat[:,1], ls='-', label='Bm = 1e6')

dat = loadtxt("data/primakoff_profile_cham_1e7.dat")
#dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
ax2.plot(dat[:,0]/R,dat[:,1], ls='-', label='Bm = 1e7')

#dat = loadtxt("data/primakoff_profile_cham_1e8.dat")
##dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
#ax2.plot(dat[:,0]/R,dat[:,1], ls='-', label='Bm = 1e8')


# axes
ax2.set_xlabel("Solar radius fraction")
ax2.set_ylabel("Emission rate dN/dr [eV2 Lambda2]")	#[m-2 s-1 eV-1]")
#ax2.set_xscale('log')
#ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/primakoff_profile.jpg')
plt.show()
