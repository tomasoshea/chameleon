# comparing different contributions for babyIAXO

import numpy as np
from numpy import loadtxt
from matplotlib import pyplot as plt
from numpy import power as pow

plt.style.use("style.txt")	# import plot style

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8), projection='3d')
#ax2.set(xlim=(1e-2,1.01), ylim=(1e-3,1.01))
Mpl = 2e27		# planck mass [eV]
rho=6.486746444817414226e+20	# density in solar core [eV4]

# cham mass
def mCham2( Bm, E, n ):
	E4n = pow(E,4+n)
	return n*(n+1)*E4n* pow( Bm*rho/(n*Mpl*E4n), (n+2)/(n+1) )
np.vectorize(mCham2)

Beta = np.logspace(1e0,1e5,100)
Lambda = np.logspace(1e-6,1e0,100)
n = np.logspace(1e0,1e2,100)
x=[]
y=[]
z=[]

for B in Beta :
	for L in Lambda :
		for N in n :
			if mCham2(B,L,N)<1: continue
			x.append(B)
			y.append(L)
			z.append(N)
			break

			
ax2.plot_surface(x,y,z)

# axes
#ax2.set_xlabel("Fraction of solar radius")
#ax2.set_ylabel("Solar parameters (normalised)")
#ax2.set_xscale('log')
#ax2.set_yscale('log')
#ax2.legend()

plt.savefig('plots/3Dparams.jpg')
plt.show()
