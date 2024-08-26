# Tom O'Shea 2024

# plot Lambda against Bm chameleon bounds

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
#ax2 = fig2.add_axes((.15,.15,.8,.8))
ax2 = fig2.subplots()
#ax2.set(xlim=(1e-1,1e18), ylim=(1e-8, 1e1))
ax2.set(xlim=(1e-1,1e7), ylim=(1e-7, 1e-1))
v = np.logspace(1e-20,1e0,nticks)

# plot DE scale
ax2.hlines(2.4e-3, 1e-1, 1e18, color='black', ls='-',zorder=10)
ax2.text(1.1e-1,3e-3,r'$\Lambda=\Lambda_{\mathrm{DE}}$',zorder=10)
#ax2.text(1e11,3e-3,r'$\Lambda=\Lambda_{\mathrm{DM}}$')

# plot "cutoff"
Lam = np.logspace(-18,1,100)
def Beta2(L):
	return 4.85532e9 * pow(L,1.53724)
np.vectorize(Beta2)

def Beta(L):
	return pow(10,9.49285) * pow(L,1.66571)
np.vectorize(Beta)
#Bm = 1/Beta(Lam)
ax2.fill_between(Beta(Lam),Lam,1e15,ec='black',ls='-',label='This work',color='r',alpha=0.5,hatch='++',zorder=0.1)
ax2.text(3e-1,1e-5,"This work")
#ax2.text(1e0,1e-5,"This work")

# plot other bounds
#dat = np.loadtxt("data/casimir.dat", delimiter=',')
#print(dat[:,1])
#ax2.plot(dat[:,0],dat[:,1],ls=':',label='Casimir')

dat = np.loadtxt("data/interferometry.dat", delimiter=',')
ax2.fill_between(1/dat[:,0],dat[:,1],1e15,ec='black',ls='-',label='AI')
ax2.text(1e6,1e-4,"AI")

dat = np.loadtxt("data/lfs.dat", delimiter=',')
ax2.fill_between(1/dat[:,0],dat[:,1],1e15,ec='black',ls='-',label='LFS',zorder=0.2)
ax2.text(5e1,6e-4,"LFS")

dat = np.loadtxt("data/torsionbalance.dat", delimiter=',')
ax2.fill_between(1/dat[:,0],dat[:,1],1e15,ec='black',ls='-',label='TB')
ax2.text(3e1,2e-2,"TB")

# axes
ax2.set_xlabel(r'$\beta_m$')
ax2.set_ylabel(r'$\Lambda$ [eV]')	#[m-2 s-1 eV-1]")
ax2.set_xscale('log')
ax2.set_yscale('log')
#ax2.legend(loc='lower right')

plt.tight_layout()
name = "casimirplot_Bm--fill--reduced--2"
plt.savefig('plots/{}.jpg'.format(name))
plt.savefig('plots/pdfs/{}.pdf'.format(name))
plt.show()
