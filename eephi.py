# Tom O'Shea 2024

# python script to plot resonance of w_p as fn of E, q
# from https://arxiv.org/pdf/2312.16306

import numpy as np
from scipy.special import dawsn
from matplotlib import pyplot as plt
from matplotlib import ticker
from scipy import integrate

plt.style.use("style.txt")

# constants
pi = 3.14159265359
alpha = 1/137.035999084
me = 510998.950		# e- mass [eV]

# solar densities etc
ne = np.loadtxt("data/ne.dat")
nbar = np.loadtxt("data/nbar.dat")
Ts = np.loadtxt("data/T.dat")
r = np.loadtxt("data/r.dat")
wp = np.loadtxt("data/wp.dat")

c = 500#195		# r = 0.1 rSolar
m = 0

# real part of self-E
def RePi(E, q, ni, mi, T):
	return -2*ni/q * np.sqrt(mi/2/T) * ( dawsn(np.sqrt(mi/2/T)*(E/q + q/2/mi)) - dawsn(np.sqrt(mi/2/T)*(E/q - q/2/mi)) )
RePi = np.vectorize(RePi)

# imaginary part of self-E
def ImPi(E, q, ni, mi, T):
	return -ni*np.sqrt(2*pi/mi/T)*mi/q * np.exp( -(np.power(mi*E/q,2) + np.power(q/2,2))/(2*mi*T) ) * np.sinh(E/2/T)
ImPi = np.vectorize(ImPi)

# electron F(E,q)
def F(E,q):
	Ve = 4*pi*alpha/q**2
	T = Ts[c]
	ni = ne[c]
	iPi = ImPi(E,q,ni,me,T)
	rPi = RePi(E,q,ni,me,T)
	return -2*Ve/(1-np.exp(-E/T)) * iPi / ( (1-Ve*rPi)**2 + (Ve*iPi)**2 )
#F = np.vectorize(F)

def xintegrand(x, E, q, w):
	k = np.sqrt(w*w - m*m)
	y = q/k
	u = 0.5*(y + 1/y)
	qi = np.sqrt(k*k + q*q - 2*x*q*k)
	return np.power(x-y,2)/(u-x) * F(w-E,qi)

def integrand(E, q, w):
	k = np.sqrt(w*w - m*m)
	I = integrate.quad(xintegrand,-1,1,args=(E, q, w))
	return 1/2/(2*pi)**3 *q*k*k*F(E,q)*F(w-E,k-q) * I[0]

def rate(w):
	I = integrate.dblquad(integrand, 0, np.inf, 0, np.inf, args=(w,))
	return I[0]

# plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.subplots()
ax2.set_xlabel(r'$\omega_\phi$ [eV]')
ax2.set_ylabel(r'$\Gamma/V\phi$ [eV3]')
#ax2.set_xscale('log')
ax2.set_yscale('log')

size=200
wphi = np.linspace(300,2000,size)
G = np.zeros(size)
for i in range(size):
	print(wphi[i])
	G[i] = rate(wphi[i])
#G = rate(300)
print(G)

plt.plot(wphi,G)
plt.savefig('plots/eephi_c500.jpg')
plt.tight_layout()
plt.show()