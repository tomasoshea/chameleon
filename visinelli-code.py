import numpy as np

## Physical constants
hbar  = 6.582119569e-19 # Planck's constant  [keV s]
cs    = 2.99792458e10   # speed of light     [cm/s]
hbarc = hbar*cs         #                    [keV*cm]
me    = 511.            # Electron mass      [keV]
RSun  = 6.9634e10       # Solar radius       [cm]
Rt   = 0.7*RSun         # Tachocline radius  [cm]
dSun  = 1.5e13          # Earth-Sun distance [cm]
mPl   = 1.221e25        # Planck mass        [keV]
alpha = 1./137.036      # Fine-structure constant

## Values at tachocline
xt    = 0.7             # Relative position 
Bt    = 0.021           # Solar magnetic field    [keV^2]
rhot  = 891921.86598    # Plasma density          [keV^4]
TSun  = 0.1995759863    # Temperature             [keV]
ne    = 1.06051586e+23  # Electron number density [cm^-3]
lam   = 0.3             # Photon mean free path   [cm]
ng    = 1.e21           # Photon flux [s^-1 cm^-2; Eq.32 in 1110.2583]
Dx    = 0.01            # Thickness of the tachocline

## Values at Sun core
rhoc   = 6.5e8  # Plasma density  [keV^4]
Tcore  = 1.3    # Temperature     [keV]

## Conversion units
pc    = 3.086e18     # cm
gkeV  = 1.782662e-30 # conversion from g to keV
year  = 3.15e7       # conversion from year to s
T2eV = 2e-16 * 1e18  # Tesla in eV2

## Numerical factors
pi2   = 2.*np.pi
pi4   = 4.*np.pi
pi8   = 8.*np.pi

# CAST params
L = 9.26 # CAST length [m]
B = 9 * T2eV    # CAST B-field [eV2]
phi = 1e-5 * 1e4 * 14.5 # CAST bg flux [m-2 s-1]

# Plasma frequency squared at the tachocline in keV^2
omPl2 = pi4*alpha*ne*hbarc**3/me

# chameleon mass squared in keV^2
# Eq.3-4 in 1110.2583
def mphi2(ge, gc, n):
    phim = (n*mPl*10**((4.+n)*gc-ge)/rhot)**(1./(1.+n))
    return (1.+n)*10**ge*rhot/mPl/phim

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
    return ng*pg(om)*Pcham(om, ge, gg, gc, n, logMg)*(Rt/dSun)**2/4.


###################################################################################################
####################################### CAST BIT ##################################################
###################################################################################################

def wIntegral( n, Bm, ge, gg, gc, logMg ):

    # integrate wrt w over CAST energies (0.5 - 15 keV) by trapezia
    dw = 1e0
    item = 0.
    w = 5e2
    while w < 15e3:     # w in eV
        item += ( dw * pow(Rt/R,2) * pow(B*(L/hbarc)/(2*mPl),2) * ( dPhidomega_cham(w, ge, gg, gc, n, logMg) + dPhidomega_cham(om, ge, gg, gc, n, logMg) ) / 2 )
        w += dw
    
    return item

# set model parameter n
n = 1
# initialise vectors
BgVec = []
BmVec = []
# scan over various Bm
Bm = 1
while Bm < 1e6:
    BgVec.append( pow( phi / wIntegral(n,Bm), 0.25) )
    BmVec.append(Bm)
    Bm*=2

# set path for writeout
path = "data/limits/CAST-n"
ext = "-python.dat"
np.savetxt( path + int(n) + ext, [BmVec, BgVec] )