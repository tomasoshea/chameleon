# Tom O'Shea 2022

# filter - create CSV files from solar model to feed to C++ program

from numpy import power as pw
import numpy as np

############ constants ################

# physical constants
a = (1./137.)
pi = np.pi
me = 0.51099895000e6	# electron mass in eV
Rsun= 6.9598e+8	# solar radius in meters

# conversion factors to eV
K2eV = 8.617333262e-5	# from Boltzmann
kg2eV = 5.609588604e35	# from hbar/c2
amu = 931.49410242e6	# atomic mass unit in eV
m2eV = 1.973269804e-7	# m-1 to eV by hc = 197 MeV fm


########### open file ###############

# creates arrays of [r, T, Rho, 1H, 4He, 3He, 12C, 14N, 16O]
cols1 = np.array((1,2,3))
cols2 = np.arange(6, 35, 1)
cols = np.append(cols1, cols2)

rFrac, T, rho, H1, He4, He3, C12, C13, N14, N15, O16, O17, O18, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni = np.loadtxt('data/AGSS09_solar_model_stripped.dat', usecols=cols).swapaxes(0,1)

# convert values to eV
T = T * K2eV	# temp from K to eV
rho = rho * 1e3 * pw(m2eV, 3) * kg2eV	# density from g cm-3 to eV4
rReal = rFrac * Rsun / m2eV	# radial distance from m to eV-1


############# compute values ###############

# array of mass fractions
fracs = [H1, He4, He3, C12, C13, N14, N15, O16, O17, O18, Ne, Na, Mg, Al, Si, P, S, Cl,  Ar,  K,   Ca,  Sc,   Ti,  V,   Cr,  Mn, Fe, Co, Ni]
A = [1.00782503223, 4.00260325413, 3.01602932007, 12.0, 13.00335483507, 14.00307400443, 15.00010889888, 15.99491461957, 16.99913175650, 17.99915961286, 20.18, 22.99, 24.305, 26.982, 28.086, 30.974, 32.066, 35.453, 39.948, 39.098, 40.078, 44.956, 47.876, 50.942, 51.996, 54,938, 55.845, 58.933, 58.693]	# masses of ions https://wwwndc.jaea.go.jp/NuC/ 
Z = [1, 2, 2, 6, 6, 7, 7, 8, 8, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28]


# compute ion densities & electron density
ne = np.zeros(len(rReal))
ionDensities = np.empty((len(fracs), len(fracs[0])))

for i in range(len(rReal)):
	const = rho[i] / amu
	total = 0
	for j in range(len(fracs)):
		ions = const * fracs[j][i] / A[j]
		ionDensities[j][i] = ions	# ion number densities in eV3
		
		electrons = ions * Z[j]
		total += electrons
		
	ne[i] = total

nH, nHe4, nHe3 = ionDensities[0:3, :]


# compute plasma frequency
wp2 = np.multiply(( 4 * pi * a / (me) ) , ne)
wp = np.sqrt(wp2)


# save computed data
np.savetxt('data/rFrac.dat', rFrac, delimiter=',')
np.savetxt('data/r.dat', rReal, delimiter=',')
np.savetxt('data/T.dat', T, delimiter=',')
np.savetxt('data/ne.dat', ne, delimiter=',')
np.savetxt('data/wp.dat', wp, delimiter=',')
np.savetxt('data/nH.dat', nH, delimiter=',')
np.savetxt('data/nHe4.dat', nHe4, delimiter=',')
np.savetxt('data/nHe3.dat', nHe3, delimiter=',')
np.savetxt('data/rho.dat', rho, delimiter=',')


