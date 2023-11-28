# Tom O'Shea 2023

# plot scalar primakoff spectrum

from numpy import loadtxt
from matplotlib import pyplot as plt
import numpy as np


plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
#ax2.set(xlim=(0,1e1), ylim=(0, 1.01))
ax2.set(xlim=(1e-6,2e2), ylim=(0, 4e63))

dat = loadtxt("data/primakoff_Eloss.dat")
#dat[:,1] = dat[:,1]/np.nanmax(dat[:,1])
ax2.plot(dat[:,0]*1e-3,dat[:,1], ls='-')

# axes
ax2.set_xlabel("Scalar mass [eV]")
ax2.set_ylabel("Energy loss rate [eV2 Lambda2]")	#[m-2 s-1 eV-1]")
ax2.set_xscale('log')
#ax2.set_yscale('log')
#ax2.legend()

plt.savefig('plots/primakoff_Eloss_1.jpg')
plt.show()
