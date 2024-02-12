import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
from scipy import interpolate
from scipy.optimize import fsolve
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')

## Numerical factors
pi    = 1.*np.pi
pi2   = 2.*np.pi
pi4   = 4.*np.pi
pi8   = 8.*np.pi
sq8   = np.sqrt(pi8)
sqp   = np.sqrt(np.pi)

## Physical constants
kB    = 8.61733034e-8   # Boltzmann constant  [keV/K]
hbar  = 6.582119569e-19 # Planck's constant   [keV s]
cs    = 2.99792458e10   # speed of light      [cm/s]
hbarc = hbar*cs         #                     [keV*cm]
me    = 510.9989461     # Electron mass       [keV]
RSun  = 6.9634e10       # Solar radius        [cm]
ASun  = 1.5e13          # Earth-Sun distance  [cm]
ng    = 1.e21           # Photon flux at tachocline [s^-1 cm^-2; Eq.32 in 1110.2583]
MPl   = 1.221e25/sq8    # Reduced Planck mass [keV]
alpha = 1./137.036      # Fine-structure constant
X     = 0.7381
Y     = 0.2485
Z     = 1.-X-Y
ergtokeV = 6.24150934e8 # Conversion from erg to keV
solarlum = 6.e10        # Solar luminosity in erg s^-1 cm^-2
mu    = 1.67377e-24 
mH    = 1.00784
mHe4  = 4.00260325415
mHe3  = 3.01602931914
mC12  = 12
mC13  = 13.00335483778
mN14  = 14.00307400478
mN15  = 15.00010889823
mO16  = 15.99491461956
mO17  = 16.999131703
mO18  = 17.999161001

## Conversion units
pc      = 3.086e18      # cm
gkeV    = 1.782662e-30  # conversion from g to keV
year    = 3.15e7        # conversion from year to s
TtokeV2 = 7.e-4         # conversion from T to keV^2


# Solar magnetic field from astro-ph/0203107
# Magnetic fields in T
B0core  = 1.e4
B0tacho = 30.
B0outer = 2.
x0      = 0.712
cc      = 1.+10.*x0
Kappa   = (1.+cc)*(1.+1./cc)**cc

# x is r/RSun
def Bcore(x):
    res = 0.
    if x < x0:
      a = (x/x0)**2
      res = Kappa*B0core*a*(1.-a)**cc
    return res
Bcore = np.vectorize(Bcore)

def Btacho(x):
    btacho = 0
    if abs(x - x0) < 0.02:
      a      = (x-x0)/0.02
      btacho = B0tacho*(1.-a**2)
    return btacho
Btacho = np.vectorize(Btacho)

def Bouter(x):
    bouter = 0
    if abs(x - 0.96) < 0.035:
      a = (x - 0.96)/0.035
      bouter = B0outer*(1.-a**2)
    return bouter
Bouter = np.vectorize(Bouter)

# Total magnetic field in keV^2
def BSun(x, ccos):
    Btot = Bcore(x) + Btacho(x) + Bouter(x)
#    Btot = Btacho(x)
    return 3.*Btot*np.sqrt(1-ccos**2)*ccos*TtokeV2

## Define the location of the files containing the solar model data
## From Serenelli et al. 2009 - arxiv:0909.2668
## Model from Vinyoles et al. 2016 - arxiv:1611.09867
agss09 = '../data/struct_b16_agss09.dat'

## Read the solar model data
df = pd.read_csv(agss09, delim_whitespace=True, skiprows=9, header=None)
nf = len(df)
np.zeros(nf)

## Load the profile.
## Density dSun is in keV^4
MSun = np.zeros(nf)
xSun = np.zeros(nf)
TSun = np.zeros(nf)
dSun = np.zeros(nf)
PSun = np.zeros(nf)
LSun = np.zeros(nf)
nSun = np.zeros(nf)
delnSun = np.zeros(nf)
deldSun = np.zeros(nf)
for i in range(nf):
    MSun[i] = df[0][i]
    xSun[i] = df[1][i]
    TSun[i] = df[2][i]
    dSun[i] = df[3][i]
    PSun[i] = df[4][i]
    LSun[i] = df[5][i]
    Zp      = df[6][i]+2*(df[7][i]+df[8][i]+df[9][i]+df[10][i]+df[11][i]+df[12][i]+df[13][i]+df[14][i]+df[15][i]) 
    nSun[i] = dSun[i] / Zp / mu
    if i < nf-1:
      deldSun[i] = (dSun[i+1] - dSun[i]) / (xSun[i+1] - xSun[i])
      delnSun[i] = (nSun[i+1] - nSun[i]) / (xSun[i+1] - xSun[i])

xcore = np.log10(xSun[0])

## Density rhoSun is in keV^4
rhoSun = interpolate.interp1d(xSun, dSun/gkeV*hbarc**3, fill_value=(0,0), bounds_error=False)

#xT  = np.linspace(0.01, 1, 1000)
#B1  = rhoSun(xT)
#plt.plot(xT, B1, 'k-')
#plt.xscale('log')
#plt.yscale('log')
#plt.show()
#exit()


## derivative of rhoSun in keV^5
drhoSun = interpolate.interp1d(xSun, deldSun/gkeV*hbarc**4/RSun, fill_value=(0,0), bounds_error=False)

## Number density of electrons in keV^3
ne = interpolate.interp1d(xSun, nSun*hbarc**3, fill_value=(0,0), bounds_error=False)

## derivative of ne in keV^4
dne = interpolate.interp1d(xSun, delnSun*hbarc**4/RSun, fill_value=(0,0), bounds_error=False)

## Temperature in keV
T = interpolate.interp1d(xSun, kB*TSun, fill_value=(0,0), bounds_error=False)

#Kramers formula
def kappae(x):
    conv = 4.e25*(kB)**3.5/hbarc**3*gkeV
    return conv*(1.+X)*(1.e-3+Z)*rhoSun(x)/T(x)**3.5

#Klein-Nishina formula
def kappaTh(x):
    cnv1 = 2.7e11*kB**2/hbarc**3*gkeV
    cnv2 = 4.5e8*kB
    return 0.2*(1.+X)/(1+cnv1*rhoSun(x)/T(x)**2)/(1+(T(x)/cnv2)**0.86)

#photon mean free path in cm
def lam(x):
    return 1./gkeV*hbarc**3/(kappaTh(x) + kappae(x))/rhoSun(x)

# Chameleon field
# Parameters: betam, betag, n, Md
Lambda = 2.e-6 # keV

# Chameleon mass squared in keV^2
# Eq.3-4 in 1110.2583
def mphi2(x, ccos, betam, betag, n):
    rhoeff = betam*rhoSun(x) + betag*BSun(x, ccos)**2/2.
    phimin = Lambda*(n*MPl*Lambda**3/rhoeff)**(1./(1.+n))
    return abs((1.+n)*rhoeff/MPl/phimin)

# Plasma frequency squared in keV^2
def omPl2(x):
    return pi4*alpha*ne(x)/me

plotprofiles=0
if plotprofiles==1:
  xT  = np.linspace(0.01, 1, 1000)
  xT1  = np.linspace(0.925, 1, 1000)
  B1  = Bcore(xT)*TtokeV2
  B2  = Btacho(xT)*TtokeV2
  B3  = np.sqrt(Bouter(xT1)*TtokeV2)
  oB  = np.sqrt(omPl2(xT))
  mu1 = np.sqrt(mphi2(xT, 0.7, 1.e5, 1.e10, 1))
  mu2 = np.sqrt(mphi2(xT, 0.7, 1.e6, 1.e10, 1))
  mu3 = np.sqrt(mphi2(xT, 0.7, 1.e7, 1.e10, 1))
  Tx  = T(xT)
  plt.plot(xT, np.sqrt(B1+B2), 'k-', label=r"$B$ field")
  plt.plot(xT1, B3, 'k-')
  plt.plot(xT, oB, 'b-', label=r"$\omega_p$")
  plt.plot(xT, Tx, 'c-', label=r"$T$")
  plt.plot(xT, mu1, 'r-',  label=r"$\mu_{\rm eff},\,\beta_m=10^5$")
  plt.plot(xT, mu2, 'r--', label=r"$\mu_{\rm eff},\,\beta_m=10^6$")
  plt.plot(xT, mu3, 'r-.', label=r"$\mu_{\rm eff},\,\beta_m=10^7$")
  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel(r"Solar radius $r/R_\odot$")
  plt.ylabel(r"Profiles [$keV$]")
  plt.ylim((1.e-4, 1.e1))
  plt.legend(loc="lower left")
  plt.tight_layout()
  plt.savefig("../plots/SolarProfile.pdf", format="pdf")
  plt.show()
  exit()

def Deltagp2(x, ccos, betag):
    #b2     = (BSun(x, ccos)/M**2)**2
    return (betag*BSun(x, ccos)/MPl)**2 #/(1.+b2)
Deltagp2 = np.vectorize(Deltagp2)

def Deltaphi(om, x, ccos, betam, betag, n):
    # b2 = (BSun(x, ccos)/M**2)**2
    return mphi2(x, ccos, betam, betag, n)/2./om  # (a1-b2*om)/(1.+b2)

def Deltap(om, x):
    return omPl2(x)/2./om

def Deltaosc2(om, x, ccos, betam, betag, n):
    del2 = (Deltaphi(om, x, ccos, betam, betag, n)-Deltap(om, x))**2
    return Deltagp2(x, ccos, betag) + del2
Deltaosc2 = np.vectorize(Deltaosc2)

#xT = np.linspace(0.01,1,100)
#B1 = Deltagp2(xT, 0.7, 1.e7, 1.e10, 1.e7, 1)
#B2 = Deltaosc2(xT, 0.7, 1.e7, 1.e10, 1.e7, 1)
#B3 = B1/B2
#plt.plot(xT, B3, 'k-')
#plt.plot(xT, B2, 'b-')
#plt.plot(xT, B3, 'r-')
#plt.xscale('log')
#plt.yscale('log')
#plt.show()
#exit()


# Eq.31 in 1110.2583
# pg is in keV^-1
def pg(om, x):
    Tx = T(x)
    return om**2/Tx**3/(np.exp(om/Tx)-1.)/2.404

# Photon flux at x in cm^-2 s^-1
def flux(x):
    return ng*(x0/x)**2

def II(y0):
    return np.sqrt(pi/2.)*np.sqrt(1.+np.sqrt(1.+4.*y0**2))-np.sqrt(2.)

# Integrand function
def dprofnonres(om, x, ccos, betam, betag, n):
    mu2  = mphi2(x, ccos, betam, betag, n)
    Dos2 = Deltaosc2(om, x, ccos, betam, betag, n)
    l0   = lam(x)
    y0   = l0*np.sqrt(Dos2)/2./hbarc 
    dd   = Deltagp2(x, ccos, betag)/Dos2
    intg = 0.
    if om > np.sqrt(mu2):
      intg = flux(x)*pg(om, x)*np.sqrt(cs/l0)*RSun/l0*II(y0)*dd
    return intg
dprofnonres = np.vectorize(dprofnonres)

dplotres=0
if dplotres==1:
  xT  = np.geomspace(1.e-4, 1., 200)
  mu1 = dprofnonres(1., xT, 0.7, 1.e5, 1.e10, 1)
  mu2 = dprofnonres(1., xT, 0.7, 1.e6, 1.e10, 1)
  mu3 = dprofnonres(1., xT, 0.7, 1.e7, 1.e10, 1)
  plt.plot(xT, mu1, 'r-')
  plt.plot(xT, mu2, 'b-')
  plt.plot(xT, mu3, 'g-')
  plt.xscale('log')
  plt.yscale('log')
  plt.xlim((1.e-4, 1.))
  plt.show()
  exit()


# Return the spectrum in cm^-2 s^-1 keV^-1
def profnonres(om, betam, betag, n):
    nn = 50
    nx = 100
    cT = np.linspace(0.,  1., nn)
    xT = np.linspace(0.1, 1., nx)
    ll = np.zeros(nn)
    for i, c in enumerate(cT):
      ll[i] = np.trapz(dprofnonres(om, xT, c, betam, betag, n), xT)
    return np.trapz(ll, cT)
profnonres = np.vectorize(profnonres)

# Luminosity in erg s^-1 cm^-2
def luminosity(betam, betag, n):
    nn   = 100
    omT  = np.geomspace(0.005, 20., nn)
    intg = omT*profnonres(omT, betam, betag, n)/ergtokeV
    return np.trapz(intg, omT)
luminosity = np.vectorize(luminosity)

gridresults=1
if gridresults==1:
    Nm = 25
    Ng = 25
    Ntot   = Nm*Ng
    betamT = np.geomspace(1.e3, 1.e8,  Nm)
    betagT = np.geomspace(1.e6, 1.e10, Ng)
    lumT   = np.zeros([Ntot, 3])
    for i in range(Nm):
      for j in range(Ng):
        lumT[j+Ng*i][0] = betamT[i]
        lumT[j+Ng*i][1] = betagT[j]
        lumT[j+Ng*i][2] = luminosity(betamT[i], betagT[j], 1)
        print(j+Ng*i)
    abs_path  = './'
    data = lumT
    int_file = abs_path + 'solar_chameleon_luminosity.txt'
    np.savetxt(int_file, data)
    exit()

#print(luminosity(1.e5, 1.e10, 1))
#print(luminosity(1.e6, 1.e10, 1))
print(luminosity(1.e7, 1.e10, 1))
exit()

#print(profnonres(1, 0., 1.e10, 1))
#print(profnonres(1, 1.e1, 1.e10, 1))
#print(profnonres(1, 1.e2, 1.e10, 1))
#print(profnonres(1, 1.e3, 1.e10, 1))
#print(profnonres(1, 1.e4, 1.e10, 1))
#print(profnonres(1, 1.e5, 1.e10, 1))
#print(profnonres(1, 1.e6, 1.e10, 1))
#print(profnonres(1, 1.e7, 1.e10, 1))
#exit()

plotres=1
if plotres==1:
  omp = np.geomspace(0.01, 20, 20)
  mu1 = omp*profnonres(omp, 1.e6, 1.e10, 1)*(RSun/ASun)**2
#  mu2 = profnonres(omp, 1.e8, 1.e10, 1)*(RSun/ASun)**2
  plt.plot(omp, mu1, 'r-')
#  plt.plot(omp, mu2, 'b-')
#  plt.xscale('log')
#  plt.yscale('log')
#  plt.xlim((0.01, 30))
#  plt.ylim((1.e2, 1.e8))
  plt.show()
  exit()

##until here!!!!!

# Conversion length in units of RSun
def Lc(om, x, ccos, betag):
    return abs(2.*om*betag*BSun(x, ccos)*ne(x)*hbarc/(MPl*omPl2(x)*dne(x))/RSun)

def xconv(ccos, betam, betag, n):
    f  = lambda x: mphi2(x, ccos, betam, betag, n) - omPl2(x)
    xc = fsolve(f, [1.e-1, 1.])[0]
    return xc

# Return the spectrum in cm^-2 s^-1 keV^-1
def PhiSpectrum(om, ccos, betam, betag, n):
    xc   = xconv(ccos, betam, betag, n)
    PhiSpectrum = 0
    if xc < 1:
      lc = Lc(om, xc, ccos, betag)
      l0 = lam(xc)
      Fenh = np.sqrt(cs/l0)*RSun*lc/l0
      Pres = 0.5*Fenh*(lc*betag*BSun(xc, ccos)/MPl)**2
      flx  = flux(xc)*pg(om, xc)
      PhiSpectrum = flx*Pres
    return PhiSpectrum
PhiSpectrum = np.vectorize(PhiSpectrum)

def dprofnonres(om, x, ccos, betam, betag, n):
    mu2 = mphi2(x, ccos, betam, betag, n)
    lom = 4.*om*hbarc / mu2
    intg = 0.
    if om > np.sqrt(mu2):
        intg = (2.*om*betag*BSun(x, ccos)/MPl/mu2)**2*(RSun/lam(x))*np.sqrt(cs*pi/lom)*flux(x)*pg(om, x)
    return intg
dprofnonres = np.vectorize(dprofnonres)

def xmin(om, ccos, betam, betag, n):
    eps = 1.e-4
    f  = lambda x: abs(1. - mphi2(10.**x, ccos, betam, betag, n)/om**2) - eps
    xc = fsolve(f, [-1., 0.])[0]
    return xc

def profnonres(om, betam, betag, n):
    cT = np.linspace(0., 1., 50)
    xT = np.geomspace(0.01, 1., 500)
    ll = np.zeros(50)
    for i, c in enumerate(cT):
      ll[i] = np.trapz(dprofnonres(om, xT, c, betam, betag, n), xT)
    return np.trapz(ll, cT)
profnonres = np.vectorize(profnonres)

plotres=1
if plotres==1:
  omp = np.linspace(-1, 1.5, 100)
  omp = 10.**omp
  mu1 = profnonres(omp, 1.e6, 1.e10, 1)
  mu2 = profnonres(omp, 1.e7, 1.e10, 1)
  mu3 = profnonres(omp, 1.e8, 1.e10, 1)
  mu4 = profnonres(omp, 1.e9, 1.e10, 1)
  plt.plot(omp, mu1, 'r-')
  plt.plot(omp, mu2, 'b-')
  plt.plot(omp, mu3, 'm-')
  plt.plot(omp, mu4, 'g-')
  plt.xscale('log')
  plt.yscale('log')
  plt.show()
  exit()


def Deltaosc(x, th, betam, betag, M, n):
    del2 = (Deltaphi(x, th, betam, betag, n, M)-Deltap(x, th, betam, betag, n))**2
    return np.sqrt(4.*Deltagp(x, th, betam, betag, n, M)**2 + del2)

def theta(x, th, betam, betag, n, M):
    dgp = Deltagp(x, th, betam, betag, n, M)
    del2 = (Deltaphi(x, th, betam, betag, n, M)-Deltap(x, th, betam, betag, n))**2
    return dgp/np.sqrt(4.*dgp**2 + del2)

#print(Deltagp(0.7, np.pi, 1.e2, 1.e1, 1, 1.e-3))
#print(Deltaphi(0.7, np.pi, 1.e2, 1.e1, 1, 1.e-3))
#print(Deltap(0.7, np.pi, 1.e2, 1.e1, 1))
print(RSun*Deltaosc(0.7, np.pi, 1.e2, 1.e1, 1, 1.e-3)/hbarc)
#print(theta(0.7, np.pi, 1.e2, 1.e1, 1, 1.e-3))
#exit()

# Function accounting for the conversion length of a chameleon at tachocline
def Int(a):
    return np.sqrt(np.pi/2.)*(np.sqrt(a + np.sqrt(a**2+4.)) - np.sqrt(2.*a))

## logMg = Log10[M_gamma/keV]
def bg2(logMg):
    return Bt**4/10**(8.*logMg)

# Eq.24 in 1505.01020
# returns keV
def DeltaB(gg, logMg):
    return 2.*10**gg*Bt/mPl/np.sqrt(1.+bg2(logMg) )

# Eq.22 in 1505.01020
# returns keV
def DeltaPl(om):
    return 0.5*omPl2/om

# Eq.23 in 1505.01020
# returns keV
def Delta_a(om, ge, gc, n, logMg):
    temp = mphi2(ge, gc, n)/2./om + bg2(logMg)*om*(2./3.)
    return temp/(1.+ bg2(logMg) )

# Denominator in Eq.25 of 1505.01020
# returns keV^2
def DEN(om, ge, gg, gc, n, logMg):
    return 4.*DeltaB(gg, logMg)**2 + (DeltaPl(om) - Delta_a(om, ge, gc, n, logMg))**2

# Eq.43 in 1110.2583
# om is in keV
# returns cm
def lom(om, ge, gg, gc, n, logMg):
    return 2.*hbarc/np.sqrt(DEN(om, ge, gg, gc, n, logMg))

# Eq.31 in 1110.2583 
# pg is in keV^-1
def pg(om):
    return om**2/TSun**3/(np.exp(om/TSun)-1.)/2.404

# Eq.47 in 1110.2583
# Pcham is dimensionless
def Pcham(om, ge, gg, gc, n, logMg):
    lo = lom(om, ge, gg, gc, n, logMg)
    dd = DEN(om, ge, gg, gc, n, logMg)
    return Dx*np.sqrt(cs/lo)*RSun/lam*4.*DeltaB(gg, logMg)**2/dd*Int(lo/lam)

# Eq.48 in 1110.2583
# Returns the differential flux at Earth in cm^-2 s^-1 keV^-1
def dPhidomega_cham(om, ge, gg, gc, n, logMg):
    return ng*pg(om)*Pcham(om, ge, gg, gc, n, logMg)*(RSun/dSun)**2/4.

## Differential luminosity in cm^-2 s^-1
def dLdom(om, ge, gg, gc, n, logMg):
    return ng*pg(om)*Pcham(om, ge, gg, gc, n, logMg)*om
dLdom = np.vectorize(dLdom)

## Total luminosity in units of Lsun
def Ltot(ge, gg, gc, n, logMg):
    omT = np.geomspace(0.01, 10.0, 100)
    cov = 1.6021e-16*pi4*RSun**2
    Lsun = 4e26
    return cov*np.trapz(dLdom(omT, ge, gg, gc, n, logMg), omT)/Lsun
Ltot = np.vectorize(Ltot)


##
## DETECTION
##

## Disformal cross section at detection in cm^2/g
def sigma_cham_dis(om, logMe): # om in keV
    return 0.5*NXe*(hbarc*me*om**2/pi2/10**(4.*logMe))**2

# photoelectric cross section in cm^2/g vs energy in keV
EMeV, sigme = np.loadtxt(fs, dtype='f8', delimiter = ',', usecols=(0,1), unpack=True)
sigmae = interpolate.interp1d(EMeV*1.e3, sigme, fill_value=(sigme[0],0), bounds_error=False)

## Absorption cross section at detection in cm^2/g
def sigma_cham_abs(om, ge): # om in keV
    return (10**ge*om/mPl)**2*sigmae(om)*4./alpha

## Total cross section at detection in cm^2/g
def sigma_cham(om, ge, logMe): # om in keV
    return sigma_cham_dis(om, logMe) + sigma_cham_abs(om, ge)

####
## We multiply dPhidomega_cham * sigma_cham * conv
## dPhidomega_cham [cm^-2 s^-1 keV^-1]
## sigma_cham [cm^2/g]
## conv       [s*g*(year*ton)^-1]
## So Rate = year^-1*ton^-1*keV^-1 as in 2006.09721 Fig.3
####

# Rate before convolution
# Input energy is om in keV
# Output rate is in year^-1*ton^-1*keV^-1
def Rate_bare_cham(om, ge, gg, gc, n, logMe, logMg):
    ret = 0
    if om > 0:
      ret = np.nan_to_num(conv*sigma_cham(om, ge, logMe)*dPhidomega_cham(om, ge, gg, gc, n, logMg))
    return ret
Rate_bare_cham = np.vectorize(Rate_bare_cham)

def Rate_cham(om, ge, gg, gc, n, logMe, logMg):
    omp = np.geomspace(0.01, 10.0, 100)
    omX = Rate_bare_cham(omp, ge, gg, gc, n, logMe, logMg)
    cnv = epsilon(om)*np.trapz(res_sgm(omp, om)*omX, omp)
    return cnv
Rate_cham = np.vectorize(Rate_cham)

#AXION
#
#
xs, TF, rho = np.loadtxt('../data/bp2004stdmodel.txt', dtype='f8', usecols=(1,2,3), unpack=True)
xM=max(xs)
rhokeV4 = interpolate.interp1d(xs, rho*hbarc**3/gkeV,fill_value=(rho[0]*hbarc**3/gkeV,0.0), bounds_error=False)
tsun = interpolate.interp1d(xs, kB*TF,fill_value=(kB*TF[0],0.0), bounds_error=False)
nHSun = interpolate.interp1d(xs, X*rho/mH,fill_value=(X*rho[0]/mH,0.0), bounds_error=False)
nHeSun = interpolate.interp1d(xs,Y*rho/mHe,fill_value=(Y*rho[0]/mHe,0.0), bounds_error=False)
Tr , neF = np.loadtxt('../data/nevsT.csv', dtype='f8', delimiter = ',', usecols=(0,1), unpack=True)
neSunT = interpolate.interp1d(kB*Tr, neF, fill_value=(0.0,0.0), bounds_error=False)
def neSun(x):
    return neSunT(tsun(x))

def k_Debye(x): # T is in keV
    return np.sqrt(pi4*alpha*(neSun(x)+nHSun(x)+4.*nHeSun(x))*hbarc/tsun(x))

def k_th(x): # in keV^-1
    return np.sqrt(2.*me*tsun(x))/hbarc

def y(x):
    return k_Debye(x)/k_th(x)

def f1(z, za):
    return z*np.exp(-z**2)*np.log((za+z)/(za-z))

def Fga(w0):
    NT = 100
    zr = np.linspace(0., 50., NT)
    ll = np.zeros(NT)
    for i, z in enumerate(zr):
        za = np.sqrt(z**2+w0)
        ll[i] = f1(z, za)
    return np.trapz(ll, zr)

def Gamma_ax(x, om): #s^-1
    TT=tsun(x)
    w0=om/TT
    y0=y(x)
    den = hbarc**5/np.sqrt(TT*me**7)/om # cm^5
    ex  = np.exp(-w0)
    return cs*alpha**2*pi8/(3.*np.sqrt(pi2))*(nHSun(x)+4.*nHeSun(x))*neSun(x)*den*ex*Fga(w0)
Gamma_ax = np.vectorize(Gamma_ax)

def f_ax(x, om): # cm^-2 s^-1 keV^-1
    return pi4*(om/dSun)**2*Gamma_ax(x, om)/(pi2**3)*x**2*(RSun/hbarc)**3
f_ax = np.vectorize(f_ax)

def dPhidomega_ax(om): # cm^-2 s^-1 keV^-1
    xt=np.linspace(0.0, xM, 30)
    return np.trapz(f_ax(xt, om), xt)
dPhidomega_ax = np.vectorize(dPhidomega_ax)

def sigma_ax(om): # om in keV
    pref = (om/me)**2/pi8/alpha
    return pref*sigmae(om)

def Rate_bare_ax(om, g):
    ret = 0
    if om > 0:
      ret = np.nan_to_num(g**4*conv*sigma_ax(om)*dPhidomega_ax(om))
    return ret
Rate_bare_ax = np.vectorize(Rate_bare_ax)

def Rate_ax(om, g):
    omp = np.geomspace(0.01, 10.0, 100)
    omX = Rate_bare_ax(omp, g)
    cnv = epsilon(om)*np.trapz(res_sgm(omp, om)*omX, omp)
    return cnv
Rate_ax = np.vectorize(Rate_ax)

### MCMC

### Values chosen for the MCMC
### n_in is the exponent of the potential
### gc_in = Log10(Lambda/eV) where Lambda^4 is the DE energy density
n_in  =  1.
gc_in = -3. # fix Lambda = 1meV
gax_in = 3.e-12

### Starting points of the MCMC
logMe0 = 3.67
logMg0 = 0.
ge0    = 0.
gg0    = 9.8

omp = np.geomspace(0.1, 4.0, 20)
bk  = back(omp)
cnv_ax = Rate_ax(omp, gax_in)
cnv_ch = Rate_cham(omp, ge0, gg0, gc_in, n_in, logMe0, logMg0)
plt.plot(enb, backb, 'r-')
plt.plot(omp, bk+cnv_ax+cnv_ch, 'k-')
plt.plot(omp, bk+cnv_ax, 'b-')
plt.plot(omp, bk+cnv_ch, 'y-')
plt.errorbar(xd, yC, yerr=sigD, xerr=None, fmt='.k')
plt.xlim((0., 4.))
#plt.xscale('log')
plt.show()
