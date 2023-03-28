# Tom O'Shea 2023

# plot meff for slar chameleons

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set(xlim=(0,1))#, ylim=(1e0, 1.2e3))

# loop for various n
#for i in range(1,5):
#    c = i/4
#    dat = loadtxt("data/meff{}.dat".format(i))
#    ax2.plot(dat[:,0],dat[:,1], color=(1,c,0), label='n={}'.format(i))

dat = loadtxt("data/meff4.dat")
ax2.plot(dat[:,0],dat[:,1], color='magenta')#, label='Bm = 1e6, n=1')

# axes
ax2.set_xlabel("Solar radius fraction")
ax2.set_ylabel("Chameleon mass [eV]")
#ax2.set_xscale('log')
ax2.set_yscale('log')
#ax2.legend()

plt.savefig('plots/meff.jpg')
plt.show()
