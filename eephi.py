# Tom O'Shea 2024

# python script to plot resonance of w_p as fn of E, q
# from https://arxiv.org/pdf/2312.16306

import numpy as np
from scipy.special import dawsn
from matplotlib import pyplot as plt
from matplotlib import ticker
from scipy import integrate

import warnings
warnings.filterwarnings("error")

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
	if E <=0: return 0
	if q <=0: return 0
	#if E > q: return 0
	#if (E/q < q/2/me): return 0
	Ve = 4*pi*alpha/q/q
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
	if qi <=0: return 0
	if E>=w: return 0
	if u<1.01: return (1-x) * F(w-E,qi)
	I = np.power(x-y,2)/(u-x) * F(w-E,qi)
	#print("x={}, k={}, E={}, I={}".format(x,k,w-E,I))
	if I<0: print("u={}, x={}".format(u,x))
	return I

def integrand(E, q, w):
	k2 = w*w - m*m
	I=[0,0]
	try: I = integrate.quad(xintegrand,-1,+1,args=(E, q, w))
	except: print(E, q, w)
	#print("value: {},	relative error: {}".format(I[0],I[1]/I[0]) )
	if I[0]==0: return 0
	if I[1]/I[0]>0.1: return 0
	if I[0]<0: print("E={}, q={}, w={}, I={}, err={}".format(E,q,w,I[0],I[1]/I[0]))
	return 1/2/(2*pi)**3 *q*k2*F(E,q) * I[0]

def integrand2(qi, q, E, w, k):
	if q <= 0: return 0
	if E <= 0: return 0
	if qi <= 0: return 0
	if w-E <= 0: return 0
	I = 1/4/np.power(2*pi,3) * np.power(k*k - q*q - qi*qi, 2)/q/qi * F(E,q) * F(w-E,qi)
	if I <= 0: return 0
	return I
	
def rate(w):
	I = integrate.dblquad(integrand, 0, np.inf, 0, np.inf, args=(w,))
	return I[0]

def rate2(w,cSolar):
	c = cSolar
	k = np.sqrt(w*w - m*m)
	qmax = lambda E,q: np.sqrt(q*q + k*k + 2*k*q)
	qmin = lambda E,q: np.sqrt(q*q + k*k - 2*k*q)
	I = integrate.tplquad(integrand2, 0, w, 0, np.inf, qmin, qmax, args=(w,k))
	return I[0]

# plot
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.subplots()
ax2.set_xlabel(r'$\omega_\phi$ [eV]')
ax2.set_ylabel(r'$\Gamma/V\omega$ [eV3]')
#ax2.set_xscale('log')
ax2.set_yscale('log')

size=100
G = np.zeros(size)
wphi = np.linspace(1,600,size)
for i in range(size):
	if wphi[i] > 220 and wphi[i] < 520: continue
	print(wphi[i])
	G[i] = rate2(wphi[i],195)
#G = rate2(220,195)
#print(wp[c])
print(G)
#G = F(245,170000)

"""
w = 300
E = 200
k = np.sqrt(w*w - m*m)
qarr = np.geomspace(10,1e6,size)
for i in range(size):
	q = qarr[i]
	#print(E)
	qmax = np.sqrt(q*q + k*k + 2*k*q)
	qmin = np.sqrt(q*q + k*k - 2*k*q)
	#G[i] = integrate.dblquad(integrand2,0,np.inf,qmin,qmax,args=(E, w, k))[0]
	G[i] = integrate.quad(integrand2,qmin,qmax,args=(q, E, w, k))[0]
print(G)
"""

ymin = np.nanmin(G)
ymax = np.nanmax(G)

plt.vlines(wp[c],ymin,ymax,color='orange')

#ax2.plot(qarr,G)
plt.plot(wphi,G)
#plt.savefig('plots/eephi_Eintg.jpg')
plt.tight_layout()
plt.show()

