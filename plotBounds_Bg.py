# Tom O'Shea 2023
# plot IAXO chameleon-photon coupling limits

from numpy import *
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle


# constants
Mpl = 2e27		# Planck mass [eV]
Lsolar = 3.0128e28		# Solar luminosity [W]
s2eV = (6.582119569e-16)#	// Hz to eV
J2eV = (1. / 1.602176634e-19)#	// Joules to eV (1 / e)
m2eV = (1.973269804e-7)#		inverse meters to eV
Lsolar *= J2eV*s2eV		# solar luminosity [eV2]


# IAXO expected background	[eV3 keV-1]
babybkg = 1e-7*1e4*m2eV*m2eV*s2eV
basebkg = 1e-8*1e4*m2eV*m2eV*s2eV
plusbkg = 1e-9*1e4*m2eV*m2eV*s2eV
castbkg = 1e-5*1e4*m2eV*m2eV*s2eV


# setup plot
plt.style.use("style.txt")	# import plot style
#fig, (ax1,ax2) = plt.subplots(1, 2)	# display is 1920 x 1080 (16:9)
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
#ax1.set(xlim=(1e0, 1e4), ylim=(1e8, 1e11))
ax2.set(xlim=(1e0, 1e4), ylim=(1e7, 2e11))


########################################
############### n=1 ####################
########################################

# add CAST
ax2.hlines(5.74e10, 1e0, 1e11, ls=':', label='CAST (old)')

dat = loadtxt("data/primakoff_total_Eloss_n1.dat")
#dat = loadtxt("data/primakoff_total_Eloss_tach.dat")
Bg = sqrt(Lsolar*3/100/dat[:,1])
ax2.plot(dat[:,0],Bg, ls='--', color='m',label='Solar energy loss')

dat = loadtxt("data/CAST_totalflux.dat")
#dat = loadtxt("data/CAST_Brax_totalflux.dat")
#dat = loadtxt("data/CAST_old_totalflux.dat")
Bg = power(castbkg/dat[:,1], 1/4)
ax2.plot(dat[:,0],Bg, ls='-', label='CAST (new)')

dat = loadtxt("data/babyIAXO_totalflux.dat")
Bg = power(divide(babybkg,dat[:,1]), 1/4)
ax2.plot(dat[:,0],Bg, ls='-', label='BabyIAXO')

dat = loadtxt("data/baseIAXO_totalflux.dat")
Bg = power(divide(basebkg,dat[:,1]), 1/4)
ax2.plot(dat[:,0],Bg, ls='-', label='IAXO')

dat = loadtxt("data/plusIAXO_totalflux.dat")
Bg = power(divide(plusbkg,dat[:,1]), 1/4)
ax2.plot(dat[:,0],Bg, ls='-', label='IAXO+')


# add other limits
ax2.add_patch( Rectangle( (2e1, 1e-3), -1e11, 1e12, color='r', alpha=0.4, label='Torsion balance') )
ax2.add_patch( Rectangle( (3e2, 1e-3), 1e11, 1e12, color='b', alpha=0.4, label='Atom interferometry') )
ax2.add_patch( Rectangle( (5.88, 1e-3), 619, 1e12, color='g', alpha=0.4, label='Levitated force sensor') )


"""
########################################
############### n=4 ####################
########################################

dat = loadtxt("data/primakoff_total_Eloss_n4.dat")
Bg = sqrt(Lsolar*3/100/dat[:,1])
ax2.plot(dat[:,0],Bg, ls='-',color='r', label='Solar energy loss')

dat = loadtxt("data/babyIAXO_totalflux.dat")
Bg = power(divide(bkg,dat[:,1]), 1/4)
ax2.plot(dat[:,0],Bg, ls='-',color='g', label='BabyIAXO')

# add CAST
ax2.hlines(5.74e10 / 4, 1e0, 1e11, color='black', ls='-', label='CAST')

# add other limits
ax2.add_patch( Rectangle( (1e1, 1e-10), -1e11, 1e12, color='r', alpha=0.4, label='Torsion balance') )
ax2.add_patch( Rectangle( (3e3, 1e-10), 1e11, 1e12, color='b', alpha=0.4, label='Atom interferometry') )
#ax2.add_patch( Rectangle( (3e1, 1e7), 7e1, 1e12, color='g', alpha=0.4, label='Levitated force sensor') )
ax2.add_patch( Rectangle( (3e100, 1e-10), 7e1, 1e12, color='g', alpha=0.4, label='Levitated force sensor') )
"""
# axes
ax2.set_xlabel(r'Matter coupling $\beta_m$')
ax2.set_ylabel(r'Photon coupling $\beta_\gamma$')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()#loc='lower right')

plt.savefig('plots/chamLimits_IAXO.jpg')
plt.show()
