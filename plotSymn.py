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
ax2.set(xlim=(1e3, 1e22))#, ylim=(1e8, 1e11))

# add IAXO bits
path = "data/limits/symstats"

dat = loadtxt("{}-babyIAXO4.dat".format(path))
ax2.plot(dat[:,0],dat[:,1], color='magenta', label='babyIAXO')
print("babyIAXO:	{}".format(average(dat[:,1])))

dat = loadtxt("{}-baselineIAXO4.dat".format(path))
ax2.plot(dat[:,0],dat[:,1], color='cyan', ls=':', label='baseline IAXO')
print("baseline IAXO:	{}".format(average(dat[:,1])))

dat = loadtxt("{}-upgradedIAXO4.dat".format(path))
ax2.plot(dat[:,0],dat[:,1], color='green', ls='--', label='upgraded IAXO')
print("upgraded IAXO:	{}".format(average(dat[:,1])))

# add cast
#ax2.hlines(5.74e10, 1e0, 1e6, color='black', ls=':', label='CAST')
#
## add other limits
#ax2.add_patch( Rectangle( (2e1, 1e8), -1e8, 1e12, color='r', alpha=0.4, label='Torsion balance (n=1)') )
#ax2.add_patch( Rectangle( (3e2, 1e8), 1e8, 1e12, color='b', alpha=0.4, label='Atom interferometry (n=1)') )
#ax2.add_patch( Rectangle( (1e0, 1e9), 1e8, 1e12, color='g', alpha=0.4, label='Astronomical polarisation') )


# axes
ax2.set_xlabel("Matter coupling, Bm")
ax2.set_ylabel("Photon coupling, Bg")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/symLimits4.jpg')
plt.show()