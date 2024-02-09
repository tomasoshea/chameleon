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
ax2.set(xlim=(1e0,1e10), ylim=(0, 1.4e10))

dat = loadtxt("data/primakoff_Eloss_n1.dat")
ax2.plot(dat[:,0],dat[:,1], ls='--',label="n=1")

dat = loadtxt("data/primakoff_Eloss_n10.dat")
ax2.plot(dat[:,0],dat[:,1], ls='--',label="n=10")

dat = loadtxt("data/primakoff_Eloss_n100.dat")
ax2.plot(dat[:,0],dat[:,1], ls='--',label="n=100")

# axes
ax2.set_xlabel("Chameleon matter coupling")
ax2.set_ylabel(r'$L_{\phi} \beta_\gamma^{-2}$ [eV2]')	#[m-2 s-1 eV-1]")
ax2.set_xscale('log')
#ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/primakoff_Eloss_n.jpg')
plt.show()
