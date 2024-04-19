# Tom O'Shea 2024

# plot scalar primakoff spectrum

from numpy import loadtxt
from matplotlib import pyplot as plt
import numpy as np
from numpy import power as pow


plt.style.use("style.txt")	# import plot style
m2eV = 1.973269804e-7		# m-1 in eV
rSolar = 6.957e8/m2eV		# solar radius in eV-1
nticks=100

# setup plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
ax2.set(xlim=(1e-1,1e18), ylim=(1e-8, 1e1))
v = np.logspace(1e-20,1e0,nticks)

# plot DE scale
ax2.hlines(2.4e-3, 1e-1, 1e18, color='black', ls='-',zorder=10)

# plot "cutoff"
Lam = np.logspace(-18,1,100)
def Beta(L):
	return 4.85532e9 * pow(L,1.53724)
np.vectorize(Beta)
#Bm = 1/Beta(Lam)
ax2.plot(Beta(Lam), Lam,ls='-',label='Solar',color='r')

# plot other bounds
dat = np.loadtxt("data/casimir.dat", delimiter=',')
print(dat[:,1])
#ax2.plot(dat[:,0],dat[:,1],ls=':',label='Casimir')

dat = np.loadtxt("data/interferometry.dat", delimiter=',')
ax2.plot(1/dat[:,0],dat[:,1],ls=':',label='Interferometry')

dat = np.loadtxt("data/lfs.dat", delimiter=',')
ax2.plot(1/dat[:,0],dat[:,1],ls=':',label='LFS')

dat = np.loadtxt("data/torsionbalance.dat", delimiter=',')
ax2.plot(1/dat[:,0],dat[:,1],ls=':',label='Torsion Balance')

# axes
ax2.set_xlabel(r'$\beta_m$')
ax2.set_ylabel(r'$\Lambda$ [eV]')	#[m-2 s-1 eV-1]")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()

plt.savefig('plots/casimirplot_Bm--full.jpg')
plt.show()












"""
# directly input integrand


Mpl = 2e27		# Planck mass [eV]
s2eV = (6.582119569e-16)#	// Hz to eV
J2eV = (1. / 1.602176634e-19)#	// Joules to eV (1 / e)
m2eV = (1.973269804e-7)#	// m-1 to eV
K2eV = (8.617333262e-5)#	// Kelvin to eV
kg2eV = 5.609588604e35#	// from hbar/c2
T2eV = 2e-16 * 1e18#		// Tesla to eV2 conversion [eV2/T]
pi = np.pi
alpha = 1/137
me = 510998.950
n=1
E = 2e-3


def mCham2( Bm, rho ):
	if n<1:	return 0
	E4n = pow(E,4+n)
	return n*(n+1)*E4n*pow( Bm*rho/(n*Mpl*E4n), (n+2)/(n+1) )
np.vectorize(mCham2)

def curlyI( u, v ):
	return (u*u - 1)/v*np.log((u-1)/(u+1)) - (pow(u+v,2) - 1)/v*np.log((u+v-1)/(u+v+1)) - 2

def curlyIapprox( u, v ):		#for u->1
	return u*u/v - (v+2)*np.log(v/(v+2)) - 2

# differential scalar production rate on earth d2N/dr/dw times Lambda2
# units eV Bm-2
def integrand(Bm, w, r, ne, T,rho):
	#if T==0:	return 0					# solves weird behaviour when ne = T = 0
	mg2 = 4*pi*alpha*ne/me		# assume mg2 = wp2
	ms2 = mCham2(Bm,rho)				# chameleon mass2 [eV2]
	#double ms2 = Bm*Bm						# fixed scalar mass2 [eV2]
	if w*w <= mg2:	return 0
	if w*w <= ms2:	return 0
	K2 = 8*pi*alpha*ne/T			# Debye screening scale ^2 [eV2]
	kgamma = np.sqrt(w*w - mg2)				# photon momentum [eV]
	kphi = np.sqrt(w*w - ms2)				# scalar momentum [eV]
	uArg = kgamma/(2*kphi) + kphi/(2*kgamma)# u for curlyI
	vArg = K2/(2*kphi*kgamma)		 		# v for curlyI
	Iuv = curlyI(uArg,vArg);
	if uArg < 1.01: Iuv = curlyIapprox(uArg,vArg)

	return alpha/(18*Mpl*Mpl*pi) * (r**2) * ne/(np.exp(w/T) - 1) * w*w*w * kphi/kgamma/kgamma * Iuv		# [eV Bg^-2]
np.vectorize(integrand)

T = loadtxt("data/T.dat")
ne = loadtxt("data/ne.dat")
rho = loadtxt("data/rho.dat")
"""
