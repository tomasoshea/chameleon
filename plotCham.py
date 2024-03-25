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
Lsolar *= J2eV*s2eV

# babyIAXO expected background
babybkg = 1e-7*1e4*m2eV*m2eV*s2eV
basebkg = 1e-8*1e4*m2eV*m2eV*s2eV
plusbkg = 1e-9*1e4*m2eV*m2eV*s2eV
castbkg = 1e-5*1e4*m2eV*m2eV*s2eV


# average function
def avg(array1, array2, cutoff):
    total = 0
    c = 0
    for i in range(len(array1)):
        if array1[i] < cutoff:	total += array2[i]
        else:
            c = i
            break
    return total / c 

# setup plot
plt.style.use("style.txt")	# import plot style
#fig, (ax1,ax2) = plt.subplots(1, 2)	# display is 1920 x 1080 (16:9)
fig2 = plt.figure(1)	# display is 1920 x 1080 (16:9)
ax2 = fig2.add_axes((.1,.1,.8,.8))
#ax1.set(xlim=(1e0, 1e4), ylim=(1e8, 1e11))
#ax2.set(xlim=(1e0, 1e4), ylim=(1e-1, 1e11))

# add IAXO bits
path = "data/limits/chamstats"

# loop over n values
#for n in range(1,11):
#	dat = loadtxt("{}-babyIAXO{}-tPlasmon.dat".format(path, n))
#	ax1.plot(dat[:,0],dat[:,1], label='n = {}'.format(n))
#	print("babyIAXO:	{}".format(avg(dat[:,0], dat[:,1], 1e3)))

#dat = loadtxt("{}-babyIAXO1-tPlasmon.dat".format(path))
#ax1.plot(dat[:,0],dat[:,1], color='magenta', ls='-.', label='baby IAXO')
#print("baby IAXO:	{}".format(avg(dat[:,0], dat[:,1], 1e3)))

#dat = loadtxt("{}-baselineIAXO1-tPlasmon.dat".format(path))
#ax1.plot(dat[:,0],dat[:,1], color='cyan', ls=':', label='baseline IAXO')
#print("baseline IAXO:	{}".format(avg(dat[:,0], dat[:,1], 1e3)))

#dat = loadtxt("{}-upgradedIAXO1-tPlasmon.dat".format(path))
#ax1.plot(dat[:,0],dat[:,1], color='green', ls='--', label='upgraded IAXO')
#print("upgraded IAXO:	{}".format(avg(dat[:,0], dat[:,1], 1e3)))

########################################
############### n=1 ####################
########################################

dat = loadtxt("data/primakoff_total_Eloss_n1.dat")
Bg = sqrt(Lsolar*3/100/dat[:,1])
ax2.plot(dat[:,0],Bg, ls='-', color='m',label='Solar energy loss')


#dat = loadtxt("data/CAST_totalflux.dat")
dat = loadtxt("data/CAST_old_totalflux.dat")
Bg = power(divide(castbkg,dat[:,1]), 1/4)
ax2.plot(dat[:,0],Bg, ls='-', label='CAST (old, calcd)')

dat = loadtxt("data/babyIAXO_totalflux.dat")
Bg = power(divide(babybkg,dat[:,1]), 1/4)
ax2.plot(dat[:,0],Bg, ls='-', label='BabyIAXO')

dat = loadtxt("data/baseIAXO_totalflux.dat")
Bg = power(divide(basebkg,dat[:,1]), 1/4)
ax2.plot(dat[:,0],Bg, ls='-', label='IAXO')

dat = loadtxt("data/plusIAXO_totalflux.dat")
Bg = power(divide(plusbkg,dat[:,1]), 1/4)
ax2.plot(dat[:,0],Bg, ls='-', label='IAXO+')

# add CAST
ax2.hlines(5.74e10 / 4, 1e0, 1e11, ls='-', label='CAST (old)')

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

plt.savefig('plots/chamLimits_CASTold.jpg')
plt.show()
