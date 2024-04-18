# Tom O'Shea 2023

# plot scalar primakoff spectrum

from numpy import loadtxt
from matplotlib import pyplot as plt
import numpy as np


plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
#ax2.set(xlim=(0,4), ylim=(0, 8e-13))
n = 1
Lambda = 1
rho = 1

def V(x): return (Lambda**(n+4))/(x**n)
np.vectorize(V)

def U(x): return rho*x
np.vectorize(U)

phi = np.linspace(1e-1,10,1000)
ax2.plot(phi,V(phi), ls='--', label=r'$V(\phi)$')
ax2.plot(phi,U(phi), ls=':', label=r'$\rho \phi$')
ax2.plot(phi,V(phi)+U(phi), ls='-', label=r'$V_{\mathrm{eff}}(\phi)$')



 
# axes
ax2.set_ylabel(r'$V(\phi)$')
ax2.set_xlabel(r'$\phi$')	#[m-2 s-1 eV-1]")
ax2.set_xticks([])
ax2.set_yticks([])
#ax2.set_yscale('log')
#ax2.legend()

plt.savefig('plots/Vphi_Vrho_test--nolabel.jpg')
plt.show()
