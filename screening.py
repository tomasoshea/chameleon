# Tom O'Shea 2024

# python script to plot resonance of w_p as fn of E, q
# from https://arxiv.org/pdf/2312.16306

import numpy as np
from scipy.special import dawsn
from matplotlib import pyplot as plt
from matplotlib import ticker

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

c = 195		# r = 0.1 rSolar


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

# plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.subplots()
#ax2.set(xlim=(2e-3,2e1), ylim=(1e14, 1e24))
size = 100
E = np.linspace(-250,-230,size)	# eV
q = np.linspace(1e3,2e3,size)	# eV
conts = np.zeros((size,size))
for i in range(size):
	for j in range(size):
		val = F(E[j],q[i])
		if val > 1e-5: conts[i,j] = F(E[i],q[j])

#print(conts)
lines = [1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4]
cs = ax2.contourf(q/1e3,E,conts,lines, locator=ticker.LogLocator())
ax2.set_xlabel(r'$q$ [keV]')
ax2.set_ylabel(r'$E$ [eV]')
fig2.colorbar(cs)
#plt.savefig('plots/screening-lowq.jpg')
plt.show()


