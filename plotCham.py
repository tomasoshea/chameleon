# Tom O'Shea 2023
# plot IAXO chameleon-photon coupling limits

from numpy import *
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle

# setup plot

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set(xlim=(1e0, 1e6), ylim=(1e8, 1e11))

# add IAXO bits

dat = loadtxt("data/limits/chamstats-babyIAXO-tPlasmon.dat")
ax2.plot(dat[:,0],dat[:,1], color='magenta', label='n=1')

dat = loadtxt("data/limits/chamstats-babyIAXO4-tPlasmon.dat")
ax2.plot(dat[:,0],dat[:,1], color='cyan', ls=':', label='n=4')

dat = loadtxt("data/limits/chamstats-babyIAXO10-tPlasmon.dat")
ax2.plot(dat[:,0],dat[:,1], color='green', ls='--', label='n=10')

# add cast
ax2.hlines(5.74e10, 1e0, 1e6, color='black', ls=':', label='CAST')

# add other limits
ax2.add_patch( Rectangle( (2e1, 1e8), -1e8, 1e12, color='r', alpha=0.4, label='Torsion balance (n=1)') )
ax2.add_patch( Rectangle( (3e2, 1e8), 1e8, 1e12, color='b', alpha=0.4, label='Atom interferometry (n=1)') )
ax2.add_patch( Rectangle( (1e0, 1e9), 1e8, 1e12, color='g', alpha=0.4, label='Astronomical polarisation') )


# axes
ax2.set_xlabel("Matter coupling, Bm")
ax2.set_ylabel("Photon coupling, Bg")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/chamLimits-baby.jpg')
plt.show()
