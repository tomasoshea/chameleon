// Tom O'Shea 2024

// Primakoff production of scalars in the sun
// V3 - the one using full ee scatter with full disp. rel. but no Raffelt screening=

#include "utils.h"

using namespace std;


// solar parameters
double B0 = 3e3*T2eV;	//200*T2eV;						// radiative zone max B [eV2]
double B1 = 50*T2eV;							// tachocline max B [eV2]  4*T2eV;//
double B2 = 3*T2eV;								// outer region max B [eV2]  3*T2eV;//
double r0 = 0.712*rSolar;						// [eV-1]
double r1 = 0.732*rSolar;						// [eV-1]
double d1 = 0.02*rSolar;						// [eV-1]
double r2 = 0.96*rSolar;						// [eV-1]
double d2 = 0.035*rSolar;						// [eV-1]
double ngamma0 = 1e25*m2eV*m2eV*s2eV;			// photon flux at r0 [eV3]

// solar model
vector<double> ne = read("data/ne.dat");		// electron number density [eV3]
vector<double> nbar = read("data/nbar.dat");	// Z2-summed number density [eV3]
vector<double> nbar2 = read("data/nbar2.dat");	// Z2-summed number density minus electrons [eV3]
vector<double> wp = read("data/wp.dat");		// plasma frequency [eV]
vector<double> T = read("data/T.dat");			// solar temperature [eV]
vector<double> r = read("data/r.dat");			// radial distance [eV-1]
vector<double> rho = read("data/rho.dat");		// solar density [eV4]
vector<double> nH = read("data/nH.dat");		// H number density [eV3]
vector<double> nHe3 = read("data/nHe3.dat");	// He3 number density [eV3]
vector<double> nHe4 = read("data/nHe4.dat");	// He4 number density [eV3]
// get gaunt factors
vector<vector<double>> z1 = readGaunt("data/Z1.dat");	// gaunt factors for Z=1
vector<vector<double>> z2 = readGaunt("data/Z2.dat");	// gaunt factors for Z=2

// parameters
double E = 2.4e-3;    							// cham potential energy scale [eV] (2.4e-3 for cosmological const)
double n = 1;									// chameleon potential 1/phi^n
double Bm = 1e2;								// cham matter coupling
bool tachoclining = false;						// for CAST comparison



// chameleon mass as a function of solar radius and model parameters
// cham model params n (phi-n potential), Bm (matter coupling)
// assume rho dominated by matter density
double mCham2( int c, double Bm ) {
	//if(n==0) { cout<<"ERROR! n = 0"<<endl; return 0; }
	double E4n = pow(E,4+n);
	if(n<0) { return n*(n+1)*E4n*pow( pow( Bm*rho[c]/(n*Mpl*E4n), (n+2) ) , 1/(n+1) ); }
	else { return n*(n+1)*E4n* pow( Bm*rho[c]/(n*Mpl*E4n), (n+2)/(n+1) ); }
}


///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// ELECTRON-ION PRIMAKOFF ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


// integral I(u,v) solved analytically in cross section calc
// dimensionless, with dimensionless arguments
double curlyI( double u, double v ) {
	return (u*u - 1)/v*log((u-1)/(u+1)) - (pow(u+v,2) - 1)/v*log((u+v-1)/(u+v+1)) - 2;
}

// I(u,v) in u=>1 limit
double curlyIapprox( double u, double v ) {		// for u->1
	return u*u/v - (v+2)*log(v/(v+2)) - 2;
}


// similarly solved integral to I(u,v) for for L-case
// D(y,u,v) in the text but only one of u or y is needed
double curlyD( double y, double v) {
	double u = 1/(2*y) + y/2;
	return pow(u-y,2)/v *log((u+1)/(u-1))
			- pow(u+v-y,2)/v*log((u+v+1)/(u+v-1))
			+ 2;
}


// differential scalar production rate on earth d2N/dr/dw divided by beta_gamma^2
// transverse case
// units eV Bg-2
double T_integrand( int c, double Bm, double w ) {
	if( T[c]==0 ) { return 0; }				// solves weird behaviour when ne = T = 0
	double mg2 = 4*pi*alpha*ne[c]/me;		// assume mg2 = wp2
	double ms2 = mCham2(c,Bm);				// chameleon mass2 [eV2]
	//cout<<ms2<<endl;
	//double ms2 = Bm*Bm;						// fixed scalar mass2 [eV2]
	if( w*w <= mg2 ) { return 0; }
	if( w*w <= ms2 ) { return 0; }
	double K2 = 8*pi*alpha*nbar[c]/T[c];		// Debye screening scale ^2 [eV2]
	double kgamma = sqrt(w*w - mg2);		// photon momentum [eV]
	double kphi = sqrt(w*w - ms2);			// scalar momentum [eV]
	double uArg = kgamma/(2*kphi) + kphi/(2*kgamma);	// u for curlyI
	double vArg = K2/(2*kphi*kgamma);		// v for curlyI
	double Iuv = curlyI(uArg,vArg);
	// explicitally put in the u=>1 limit to avoid badnesses in the code
	if(uArg < 1.01) { Iuv = curlyIapprox(uArg,vArg); }

	return alpha/(8*Mpl*Mpl*pi) * pow(r[c], 2) * nbar[c]/(exp(w/T[c]) - 1) 
			* w*w * kphi/kgamma * Iuv;		// [eV Bg-2]
}


// differential scalar production rate on earth d2N/dr/dk divided by beta_gamma^2
// pure longitudinal component only
// units eV Bg-2
double L_integrand( int c, double Bm, double kgamma ) {
	if( T[c]==0 ) { return 0; }				// solves weird behaviour when ne = T = 0
	//double w = wp[c];						// omega is plasma freq [eV]
	double w = sqrt( pow(wp[c],2) + 3*T[c]*kgamma*kgamma/me );	// corrected omega [eV]
	double ms2 = mCham2(c,Bm);				// chameleon mass2 [eV2]
	//double ms2 = Bm*Bm;					// fixed scalar mass2 [eV2]
	if( w*w <= ms2 ) { return 0; }
	double K2 = 8*pi*alpha*nbar[c]/T[c];		// Debye screening scale ^2 [eV2]
	//double kgamma = sqrt(w*w - mg2);		// photon momentum [eV]
	double kphi = sqrt(w*w - ms2);			// scalar momentum [eV]
	double yArg = kgamma/kphi;				// y for curlyD
	double vArg = K2/(2*kphi*kgamma);		// v for curlyD
	double Dyuv = curlyD(yArg,vArg);

	return alpha/(4*Mpl*Mpl*pi) * pow(r[c], 2) * nbar2[c] * T[c]
			* kphi * Dyuv;		// [eV Bg-2]
	// MODIFIED
	//return alpha/(4*Mpl*Mpl*pi) * pow(r[c], 2) * nbar2[c] * kphi
	//		* pow(w,3/2) * Dyuv / (exp(w/T[c]) - 1);
}

// integral over solar volume, for a given scalar mass and energy
// returns dN/dw Bg-2
// units Bg-2
double T_solarIntg( double w, double Bm ) {
	double total = 0;
	for( int c = 0; c < r.size() - 1; c++ ) {
		total += 0.5 * (r[c+1] - r[c]) * (T_integrand(c+1, Bm, w) + T_integrand(c, Bm, w));
	}
	return total;		// [Bg-2]
}

// integral over k_gamma for l-plasmon
// integrate from 0 to kappa
/// units eV2 Bg-2
double kIntg( double Bm, int c ) {
	double total = 0;
	double kD = sqrt(8*pi*alpha*nbar[c]/T[c]);
	//double kD = sqrt(8*pi*alpha*nbar[c]/me);
	double kmax = kD/10;
	double dk = kmax/1000;
	for( double k = dk; k < kmax; k+= dk ) {
		//cout << k << endl;
		total += 0.5 * dk * (L_integrand(c, Bm, k+dk) + L_integrand(c, Bm, k));
	}
	
	return total;
}

// MODIFIED
// differential scalar production rate on earth d2N/dr/dw divided by beta_gamma^2
// AS ABOVE BUT IN w NOT K
// units eV Bg-2
double L_integrand_omega( int c, double Bm, double w ) {
	double Tc = T[c];
	if( Tc==0 ) { return 0; }				// solves weird behaviour when ne = T = 0
	double wpc = wp[c];
	if(w<wpc) { return 0; }
	//double w = wp[c];						// omega is plasma freq [eV]
	double ms2 = mCham2(c,Bm);				// chameleon mass2 [eV2]
	//double ms2 = Bm*Bm;					// fixed scalar mass2 [eV2]
	if( w*w <= ms2 ) { return 0; }
	double K2 = 4*pi*alpha*nbar[c]/Tc;		// Debye screening scale ^2 [eV2]
	double kD = sqrt(4*pi*alpha*ne[c]/Tc);
	//double kgamma = sqrt(w*w - mg2);		// photon momentum [eV]
	double kphi = sqrt(w*w - ms2);			// scalar momentum [eV]
	double kgamma = sqrt( me/3/Tc * (w*w - wpc*wpc) );	// kL from corrected omega [eV]
	if( kgamma > sqrt(K2)/10 ) { return 0; }
	double yArg = kgamma/kphi;				// y for curlyD
	double vArg = K2/(2*kphi*kgamma);		// v for curlyD
	double Dyuv = curlyD(yArg,vArg);
	//if( kgamma == 0 ) { Dyuv = 0;}//32*kgamma*kgamma/3/(kphi*kphi + K2); }
	//cout << "kgamma = "<<kgamma<<endl; 
	//cout<<kgamma<<endl;
	return alpha/(4*Mpl*Mpl*pi) * pow(r[c], 2) * nbar2[c] * me/3/Tc
			* kphi *w*w/(exp(w/Tc) - 1) * Dyuv/kgamma;		// [eV Bg-2]
}

// integral over solar volume, for a given scalar mass and energy
// returns dN/dw Bg-2
// units Bg-2
double L_solarIntg( double w, double Bm ) {
	double total = 0;
	for( int c = 0; c < r.size() - 1; c++ ) {
		total += 0.5 * (r[c+1] - r[c]) * (L_integrand_omega(c+1, Bm, w) + L_integrand_omega(c, Bm, w));
	}
	return total;		// [Bg-2]
}



///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// MAGNETIC FIELD PRODUCTION /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


// Gamma_photon
// simplified to only contain plasma and free-free effects
// units eV
double GammaPhoton( double w, int c, double g1, double g2 ) {

	double p1 = 64 * pow(pi,2) * pow(alpha,3);
	double p2 = 3 * pow(me,2) * pow(w,3);
	double p3 = me * pow(ne[c],2) / (2*pi*T[c]);
	double p4 = 1 - exp(- w / T[c]);
	double p5 = 8 * pi * pow(alpha,2) * ne[c] / (3 * pow(me,2) );

	// sum of ion densities
	double ions = (nH[c] * g1) + g2 * ( (4 * nHe4[c]) + (4 * nHe3[c]) );

	return p1 * pow(p2, -1) * pow(p3, 0.5) * p4 * ions + p5;
}


// solar B field
// units eV2
double Bfield( int c ) {
	// B field in solar radiative zone
	if ( r[c] <= r0 ) {
		double l = (10 * (r0/rSolar)) + 1;
		double K = (1+l) * pow( 1 + pow(l,-1) , l ) * B0;
		return K * pow(r[c]/r0 , 2) * pow(1 - pow(r[c]/r0 , 2), l);
		}
	// B-field in tachocline
	else if ( r[c] > (r1 - d1) and r[c] < (r1 + d1)) {
		return B1 * ( 1 - pow( (r[c] - r1)/d1 , 2 ) );
		}
	// B-field in outer region
	else if ( r[c] > (r2 - d2) and r[c] < (r2 + d2) ) {
		return B2 * ( 1 - pow( (r[c] - r2)/d2 , 2 ) );
		}
	else { return 0; }
}


// selects Gaunt factor from matrix for Gamma
// returns Gamma [eV]
double selectG( int c, double w ) {
		// select g(w, T) value from matrix
		int indexT1;
		int indexT2;
		int indexX1;
		int indexX2;
		for( int i = 1; i < 200; i++ ) {
			if( z1[0][i] < T[c] and z1[0][i+1] > T[c] ) { indexT1 = i; }
			if( z2[0][i] < T[c] and z2[0][i+1] > T[c] ) { indexT2 = i; }
		}
		for( int i = 1; i < 500; i++ ) {
			if( (z1[i][0] * T[c]) < w and (z1[i+1][0] * T[c]) > w ) { indexX1 = i; }
			if( (z2[i][0] * T[c]) < w and (z2[i+1][0] * T[c]) > w ) { indexX2 = i; }
		}
		double g1 = z1[ indexT1 ][ indexX1 ];
		double g2 = z2[ indexT2 ][ indexX2 ];
		return GammaPhoton(w, c, g1, g2);
}


// differential scalar production rate d2N/dr/dw times Lambda2
// B-field conntribution
// units eV Bg-2
double B_integrand( int c, double Bm, double w ) {
	if( T[c]==0 ) { return 0; }					// solves weird behaviour when ne = T = 0
	double mg2 = 4*pi*alpha*ne[c]/me;			// assume mg2 = wp2 [eV2]
	double ms2 = mCham2(c,Bm);					// chameleon mass2 [eV2]
	//double ms2 = Bm*Bm;							// fixed scalar mass2 [eV2]
	if( w*w <= mg2 ) { return 0; }
	if( w*w <= ms2 ) { return 0; }
	double kgamma = sqrt(w*w - mg2);			// photon momentum [eV]
	double kphi = sqrt(w*w - ms2);				// scalar momentum [eV]
	double B = Bfield(c);						// solar B field [eV2]
	double G = selectG(c,w);

	return 2/(pi*Mpl*Mpl) * pow(r[c], 2) * B*B 
			* w*pow(w*w - ms2, 3/2)/( pow(ms2 - mg2, 2) + (w*w*G*G) )
			* G/(exp(w/T[c]) - 1);	// [eV Bg-2]
}


// integral over solar volume, for a given scalar mass and energy
// B-field contribution, with option to fix to tachocline for comparison
// returns dN/dw Bg-2
// units Bg-2
double B_solarIntg( double w, double Bm ) {
	double total = 0;
	for( int c = 0; c < r.size() - 1; c++ ) {
		// tachoclining
		if(tachoclining) {
			if(r[c]/rSolar < r1/rSolar-0.05) { continue; }
			else if(r[c]/rSolar > r1/rSolar+0.05) { continue; }
		}
		total += 0.5 * (r[c+1] - r[c]) * (B_integrand(c+1, Bm, w) + B_integrand(c, Bm, w));
	}
	return total;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// PHOTON COALESCENCE //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////



// differential scalar production rate on earth dN/dr times Lambda2
// LL coalescence
// units eV2 Bm-2
double integrand_ll( int c, double Bm, double kgamma ) {
	double Tc = T[c];
	double rc = r[c];
	double nec = ne[c];
	if( Tc==0 ) { return 0; }				// solves weird behaviour when ne = T = 0
	double ms2 = mCham2(c,Bm);				// chameleon mass2 [eV2]
	double mg2 = 4*pi*alpha*nec/me;		// assume mg2 = wp2
	//double ms2 = Bm*Bm;					// fixed scalar mass2 [eV2]
	double K2 = 8*pi*alpha*nec/Tc;		// Debye screening scale ^2 [eV2]
	if( 2*mg2 <= ms2 ) { cout<<ms2-(2*mg2)<<endl; return 0; }
	//return 1/(36*Mpl*Mpl*pi*pi) * pow(K2,3/2) * Tc*Tc * rc*rc * sqrt(4*mg2 - ms2);
	return 1/(8*Mpl*Mpl*pi*pi) * kgamma*kgamma * Tc*Tc * rc*rc * sqrt(4*mg2 - ms2);
			//* pow(exp(kgamma/wp[c])-1, -1);
}


// differential scalar production rate on earth d2N/dr/dw times Lambda2
// LT coalescence
// units eV Bm-2
double integrand_lt( int c, double Bm, double w ) {
	double Tc = T[c];
	double rc = r[c];
	if( Tc==0 ) { return 0; }				// solves weird behaviour when ne = T = 0
	double ms2 = mCham2(c,Bm);				// chameleon mass2 [eV2]
	double mg2 = 4*pi*alpha*ne[c]/me;		// assume mg2 = wp2
	if( w*w <= mg2 ) { return 0; }			// w > m_s requirement
	if( w*w <= ms2 ) { return 0; }			// w > m_g requirement
	//if( ms2 <= mg2 ) { return 0; }			// requirement for coalescence
	//double ms2 = Bm*Bm;					// fixed scalar mass2 [eV2]
	double wt = w - sqrt(mg2);				// t-photon omega [eV]
	double kphi = sqrt(w*w - ms2);			// scalar wavenumber [eV]
	double kt = sqrt(w*(w - 2*sqrt(mg2)));	// t-photon wavenumber [eV]
	if(w - 2*sqrt(mg2) <= 0) { return 0; }
	//cout<<kt<<endl;
	return 2*rc*rc/(Mpl*Mpl*12*pi*pi) * wt*wt * kphi * kt * Tc/(exp(wt/Tc) - 1);
}


// differential scalar production rate on earth d2N/dr/dw times Lambda2
// units eV Bm-2
double integrand_decay( int c, double Bm, double w ) {
	double Tc = T[c];
	double rc = r[c];
	if( Tc==0 ) { return 0; }				// solves weird behaviour when ne = T = 0
	double ms2 = mCham2(c,Bm);				// chameleon mass2 [eV2]
	double mg2 = 4*pi*alpha*ne[c]/me;		// assume mg2 = wp2
	if( w*w <= mg2 ) { return 0; }			// w > m_s requirement
	if( w*w <= ms2 ) { return 0; }			// w > m_g requirement
	//double ms2 = Bm*Bm;					// fixed scalar mass2 [eV2]
	double wt = w + sqrt(mg2);				// t-photon omega [eV]
	double kphi = sqrt(w*w - ms2);			// scalar wavenumber [eV]
	double kt = sqrt(w*(w - 2*sqrt(mg2)));	// t-photon wavenumber [eV]
	if(w - 2*sqrt(mg2) <= 0) { return 0; }
	//cout<<kt<<endl;
	return rc*rc/(Mpl*Mpl*12*pi*pi) * wt*wt * kphi * kt * Tc/(exp(wt/Tc) - 1);
}


// integral over solar volume, for a given scalar mass and energy
// summing LT and decay contributions
// returns dPhi/dw Bg-2 [eV2]
double solarIntg_lt( double w, double Bm ) {
	double total = 0;
		for( int c = 0; c < r.size() - 1; c++ ) {
			total += 0.5 * (r[c+1] - r[c]) * (integrand_decay(c+1, Bm, w) + integrand_decay(c, Bm, w)
												+ integrand_lt(c+1, Bm, w) + integrand_lt(c, Bm, w) );
		}
	return total;
}


// integral over k_gamma for l-plasmon
double kIntg_ll( double Bm, int c ) {
	double total = 0;
	double kMax = sqrt(2*me*wp[c])/1;
	//double kMax = 4*pi
	double dk = kMax/1000;
	for( double k = dk; k < kMax; k+= dk ) {
		//cout << k << endl;
		total += 0.5 * dk * (integrand_ll(c, Bm, k+dk) + integrand_ll(c, Bm, k));
	}
	return total;
}

// A(y,u) for coalescence
// dimensionless
double A( double y, double u ){
	return pow(u-y,2)*log((u+1)/(u-1)) - 2*u + 4*y;
}

// differential scalar production rate on earth d2N/dr/dk times Lambda2
// LL COALESCENCE WITH w1 INTG.
// units eV2 Bm-2
double integrand_ll_omega( int c, double Bm, double w, double w1 ) {
	double Tc = T[c];
	double rc = r[c];
	double nec = ne[c];
	double wpc = wp[c];
	if( Tc==0 ) { return 0; }				// solves weird behaviour when ne = T = 0
	if( w1<wpc ) { return 0; }					// min w is wp
	if( w<2*wpc ) { return 0; }					// min w is wp
	double ms2 = mCham2(c,Bm);				// chameleon mass2 [eV2]
	double kphi = sqrt(w*w - ms2);			// phi momentum [eV]
	double kgamma = sqrt( me/3/Tc * (w*w - wpc*wpc) );
	//double ms2 = Bm*Bm;					// fixed scalar mass2 [eV2]
	double K2 = 8*pi*alpha*nec/Tc;		// Debye screening scale ^2 [eV2]
	//if( kgamma >= sqrt(K2)/10 ) { cout<<"wut"<<endl;return 0; }
	double yArg = kgamma/kphi;
	double uArg = yArg/2 + 1/yArg/2;
	double Ayu = A(yArg,uArg);
	double w2 = w - w1;						// [eV]
	return 1/(16*Mpl*Mpl*pi*pi) * me/3/Tc * rc*rc * kphi*kphi
			* w1*w1/(exp(w1/Tc)-1) * w2/(exp(w2/Tc)-1) * Ayu;
}


// integral over all w1
// returns dN/dr [eV2]
double w1Intg( int c, double w, double Bm ) {
	double total = 0;
	double dw1 = 1;
	double w1 = wp[c];
	while( true ) {
		double I1 = integrand_ll_omega(c, Bm, w, w1+dw1);
		double I0 = integrand_ll_omega(c, Bm, w, w1);
		if( I1 + I0 == 0 ) { break; }
		total += 0.5 * dw1 * (I0 + I1);
		w1+=dw1;
		//cout<<I1+I0<<endl;
	}
	return total;	// [eV2]
}


// integral over solar volume, for a given scalar mass and energy
// returns N Bg-2 [eV]
double solarIntg_ll_omega( double w, double Bm ) {
	double total = 0;
	for( int c = 0; c < r.size() - 1; c++ ) {
		total += 0.5 * (r[c+1] - r[c]) * (w1Intg(c+1, Bm, w) + w1Intg(c, Bm, w));
	}
	return total;		// [eV Bg-2]
}




///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// COMPARISON WITH CAST LIMITS ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////



// get old CAST back-converted flux for old limits
// input detector params
// units eV3/keV
void CAST_old() {
	vector<double> beta, flux;
	string model = "CAST_old";
	n = 1;
	E = 2.4e-3;
	double Emin = 1e3;
	double Emax = 2e4;
	double B_cast = 9*T2eV;		// CAST B-field [eV2]
	double L_cast = 9.26/m2eV;	// CAST length [eV-1]
	double P = pow(B_cast*L_cast/2/Mpl, 2);	// back conversion prob for low m
	double dw = 1e1;
	double w1, w2, r1, r2 = 0;
	tachoclining = true;		// only look in tachocline
	for( double Bm = 1e-1; Bm <= 1e8; Bm*=10 ) {
		double total = 0;
		for( double w = Emin; w < Emax; w+=dw ){
			total += 0.5*dw*P*(B_solarIntg(w+dw,Bm)+B_solarIntg(w,Bm));
		}
		beta.push_back(Bm);
		flux.push_back(total/4/pi/dSolar/dSolar / ((Emax-Emin)/1e3));	// eV3 keV-1
		cout << "Bm = " << Bm << " of 1e18" << endl;
	}
	// write to file
	string name = "data/"+model+"_totalflux.dat";
	write2D( name, beta, flux );
}


// Brax integral
// dimensionless
double Ibrax( double a ) {
	return sqrt(pi/2)*( sqrt(a + sqrt(a*a + 4)) - sqrt(2*a) );
}


// Brax tachocline calculation
// differential emission rate (d2N/dwdr)
// units eV Bg-2
double Brax( int c, double Bm, double w ){
		if(r[c]/rSolar < r1/rSolar-0.05) { return 0; }
		else if(r[c]/rSolar > r1/rSolar+0.05) { return 0; }
	double Tc = T[c];

	// tachocline values
	double mfp_t = 0.3/100/m2eV;				// tachocline mean free path [eV-1]
	double flux_t = 1e21*1e4*m2eV*m2eV*s2eV;	// tachocline photon flux [eV3]
	double n_t = 2*zeta3*Tc*Tc*Tc/pi/pi;		// tachocline photon density [eV3]

	double mg2 = 4*pi*alpha*ne[c]/me;			// assume mg2 = wp2 [eV2]
	double ms2 = mCham2(c,Bm);					// chameleon mass2 [eV2]
	//double ms2 = Bm*Bm;						// fixed scalar mass2 [eV2]
	if( w*w <= mg2 ) { return 0; }
	if( w*w <= ms2 ) { return 0; }
	double kgamma = sqrt(w*w - mg2);			// photon momentum [eV]
	double kphi = sqrt(w*w - ms2);				// scalar momentum [eV]
	double q = abs(kgamma - kphi);				// mtm transfer [eV]
	double B = Bfield(c);						// solar B field [eV2]
	//B = 30*T2eV;

	return 2/(pi*Mpl*Mpl) * pow(r[c], 2) * B*B /(exp(w/T[c]) - 1)
			* pow(w*w/mg2, 2) * flux_t/(2*pi*n_t*mfp_t)
			* sqrt(q/2/s2eV) * Ibrax(2/mfp_t/q);	// [eV Bg-2]
}


// integral over solar volume, for a given scalar mass and energy
double solarIntgBrax( double w, double Bm ) {
	double total = 0;
	for( int c = 0; c < r.size() - 1; c++ ) {
		// tachoclining
		if(r[c]/rSolar < r1/rSolar-0.05) { continue; }
		else if(r[c]/rSolar > r1/rSolar+0.05) { continue; }

		total += 0.5 * (r[c+1] - r[c]) * (Brax(c+1, Bm, w) + Brax(c, Bm, w));
		//cout<<"r = "<<r[c]/rSolar<<"	dN/dw = "<<Brax(c, Bm, w)<<endl;
	}
	return total;
}


// get old CAST back-converted flux for old limits
// input detector params
// units eV3/keV
void CAST_Brax() {
	vector<double> beta, flux;
	string model = "CAST_Brax";
	n = 1;
	E = 2.4e-3;
	double Emin = 1e3;
	double Emax = 2e4;
	double B_cast = 9*T2eV;		// babyIAXO B-field [eV2]
	double L_cast = 9.26/m2eV;	// babyIAXO length [eV-1]
	double P = pow(B_cast*L_cast/2/Mpl, 2);
	double dw = 1e1;
	double w1, w2, r1, r2 = 0;
	for( double Bm = 1e-1; Bm <= 1e12; Bm*=10 ) {
		double total = 0;
		for( double w = Emin; w < Emax; w+=dw ){
			total += 0.5*dw*P*(solarIntgBrax(w+dw,Bm)+solarIntgBrax(w,Bm));
		}
		beta.push_back(Bm);
		flux.push_back(total/4/pi/dSolar/dSolar / ((Emax-Emin)/1e3));
		cout << "Bm = " << Bm << " of 1e18" << endl;
	}
	// write to file
	string name = "data/"+model+"_totalflux.dat";
	write2D( name, beta, flux );
}



///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// CHECKING Bg ASSUMPTIONS ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


// chameleon mass as a function of solar radius and model parameters
// cham model params n (phi-n potential), Bg (photon coupling)
// assume rho dominated by photon coupling
double mCham2_B( int c, double Bg ) {
	//if(n==0) { cout<<"ERROR! n = 0"<<endl; return 0; }
	double E4n = pow(E,4+n);
	//double B = Bfield(c);
	double B = B0;
	return n*(n+1)*E4n* pow( (2*Bm*rho[c] + Bg*B*B)/(2*n*Mpl*E4n), (n+2)/(n+1) );
}


// differential scalar production rate on earth d2N/dr/dw divided by beta_gamma^2
// transverse case
// units eV Bg-2
double T_integrand_B( int c, double Bm, double w ) {
	if( T[c]==0 ) { return 0; }				// solves weird behaviour when ne = T = 0
	double mg2 = 4*pi*alpha*ne[c]/me;		// assume mg2 = wp2
	double ms2 = mCham2_B(c,Bm);				// chameleon mass2 [eV2]
	//cout<<ms2<<endl;
	//double ms2 = Bm*Bm;						// fixed scalar mass2 [eV2]
	if( w*w <= mg2 ) { return 0; }
	if( w*w <= ms2 ) { return 0; }
	double K2 = 8*pi*alpha*nbar[c]/T[c];		// Debye screening scale ^2 [eV2]
	double kgamma = sqrt(w*w - mg2);		// photon momentum [eV]
	double kphi = sqrt(w*w - ms2);			// scalar momentum [eV]
	double uArg = kgamma/(2*kphi) + kphi/(2*kgamma);	// u for curlyI
	double vArg = K2/(2*kphi*kgamma);		// v for curlyI
	double Iuv = curlyI(uArg,vArg);
	// explicitally put in the u=>1 limit to avoid badnesses in the code
	if(uArg < 1.01) { Iuv = curlyIapprox(uArg,vArg); }

	return alpha/(8*Mpl*Mpl*pi) * pow(r[c], 2) * nbar[c]/(exp(w/T[c]) - 1) 
			* w*w * kphi/kgamma * Iuv;		// [eV Bg-2]
}


// differential scalar production rate on earth d2N/dr/dk divided by beta_gamma^2
// pure longitudinal component only
// units eV Bg-2
double L_integrand_B( int c, double Bm, double kgamma ) {
	if( T[c]==0 ) { return 0; }				// solves weird behaviour when ne = T = 0
	double w = wp[c];						// omega is plasma freq
	double ms2 = mCham2_B(c,Bm);				// chameleon mass2 [eV2]
	//double ms2 = Bm*Bm;					// fixed scalar mass2 [eV2]
	if( w*w <= ms2 ) { return 0; }
	double K2 = 8*pi*alpha*nbar[c]/T[c];		// Debye screening scale ^2 [eV2]
	//double kgamma = sqrt(w*w - mg2);		// photon momentum [eV]
	double kphi = sqrt(w*w - ms2);			// scalar momentum [eV]
	double yArg = kgamma/kphi;				// y for curlyD
	double vArg = K2/(2*kphi*kgamma);		// v for curlyD
	double Dyuv = curlyD(yArg,vArg);

	return alpha/(4*Mpl*Mpl*pi) * pow(r[c], 2) * nbar[c] * T[c]
			* kphi * Dyuv;		// [eV Bg-2]
}


// integral over solar volume, for a given scalar mass and energy
// returns dN/dw Bg-2
// units Bg-2
double T_solarIntg_B( double w, double Bm ) {
	double total = 0;
	for( int c = 0; c < r.size() - 1; c++ ) {
		total += 0.5 * (r[c+1] - r[c]) * (T_integrand_B(c+1, Bm, w) + T_integrand_B(c, Bm, w));
	}
	return total;		// [Bg-2]
}


// integral over k_gamma for l-plasmon
// integrate from 0 to kappa
/// units eV2 Bg-2
double kIntg_B( double Bm, int c ) {
	double total = 0;
	double kD = sqrt(8*pi*alpha*nbar[c]/T[c]);
	double dk = kD/1000;
	for( double k = dk; k < kD; k+= dk ) {
		//cout << k << endl;
		total += 0.5 * dk * (L_integrand_B(c, Bm, k+dk) + L_integrand_B(c, Bm, k));
	}
	return total;
}




///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// PLOTTING AND BITS /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////



// calculate emission rate profile over solar radius (dN/dr)
// integrated over relevant energies up to 20 keV
// units eV2 Bg-2
void profile( char option ) {
	vector<double> radius;
	vector<double> rate;
	string name;
	Bm = 1e2;		// cham matter coupling
	n = 1;					// cham model n
	double dw = 1e0;
	for( int c = 0; c < r.size(); c++ ) {
		if (option=='L') {
			radius.push_back(r[c]);
			rate.push_back( kIntg(Bm, c) );
		}
		else{
			double total = 0;
			for( double w = dw; w < 2e4; w+=dw ){
				if(option=='T') { total += 0.5*dw*( T_integrand(c,Bm,w+dw) + T_integrand(c,Bm,w) ); }
				if(option=='B') { total += 0.5*dw*( B_integrand(c,Bm,w+dw) + B_integrand(c,Bm,w) ); }
			}
			radius.push_back(r[c]);
			rate.push_back(total);
		}

	}
	// write to file
	if (option=='T') { name = "data/T_profile_1e2.dat"; }
	else if (option=='L') { name = "data/L_profile_1e2.dat"; }
	else if (option=='B') { name = "data/B_profile_1e2.dat"; }
	
	write2D( name , radius, rate );
}


// calculate differential particle flux spectrum dN/dw by intg over solar volume
// units Bg-2
// (dN is really dN/dt, sorry)
void spectrum( char option ) {
	vector<double> count, energy;
	string name;
	Bm = 1e2;		// cham matter coupling
	n = 1;					// cham model n
	if (option=='L') {
		double w1, w2 = 0;
		double r1, r2 = rSolar;
		for( int j = wp.size()-1; j >=0; j-- ){
			w1 = wp[j];
			if(w2 > w1) { continue; }
			else{
			r1 = r[j];
			energy.push_back(w1);
			count.push_back( kIntg(Bm, j) * abs((r2-r1)/(w2-w1)) );
			r2 = r[j];
			w2 = wp[j];
			}
		}
	}
	else{
		double dw = 1e0;
		for( double w = dw; w < 2e4; w+=dw ){
			energy.push_back(w);					// eV
			if (option=='T') { count.push_back( T_solarIntg(w,Bm) );	}	// Bg-2
			else if (option=='B') { count.push_back( B_solarIntg(w,Bm) );	}	// Bg-2
			if((int)(w) % (int)(1e3) == 0) { cout<<"w = "<<w/1e3<<"keV of 20keV"<<endl; }
		}
	}
	// write to file
	if (option=='T') { name = "data/T_spectrum_1e2.dat"; }
	else if (option=='L') { name = "data/L_spectrum_1e2--test3.dat"; }
	else if (option=='B') { name = "data/B_spectrum_1e2.dat"; }
	
	write2D( name , energy, count );
}

// as above for L-w case
void spectrumL() {
	vector<double> count, energy;
	string name;
	Bm = 1e2;		// cham matter coupling
	n = 1;					// cham model n
	double dw = 1e0;
	for( double w = dw; w < 2e4; w+=dw ){
		energy.push_back(w);					// eV
		count.push_back( L_solarIntg(w,Bm) );	// Bg-2
		if((int)(w) % (int)(1e3) == 0) { cout<<"w = "<<w/1e3<<"keV of 20keV"<<endl; }
	}
	name = "data/L_spectrum_1e2--omega_kD.dat";
	write2D( name , energy, count );	
}


// calculate differential particle flux spectrum dN/dw by intg over solar volume
// units Bg-2
// total L+T spectrum
void total_spectrum() {
	vector<double> count, energy;
	string name;
	Bm = 1e2;				// cham matter coupling
	n = 1;					// cham model n
	double w1, w2 = 0;
	double r1, r2 = rSolar;
	for( int j = wp.size()-1; j >=0; j-- ){
		w1 = wp[j];
		if(w2 > w1) { continue; }
		else{
		r1 = r[j];
		energy.push_back(w1);
		count.push_back( kIntg(Bm, j) * abs((r2-r1)/(w2-w1))
						+ T_solarIntg(w1,Bm) );
		r2 = r[j];
		w2 = wp[j];
		}
	}
	double dw = 1e0;
	for( double w = w2; w < 2e4; w+=dw ){
		energy.push_back(w);					// eV
		count.push_back( T_solarIntg(w,Bm) );	// Bg-2
		if((int)(w) % (int)(1e3) == 0) { cout<<"w = "<<w/1e3<<"keV of 20keV"<<endl; }
	}
	// write to file
	name = "data/total_spectrum_1e2.dat";
	write2D( name , energy, count );
}

// calculate total energy loss rate as a function of Bm
// both L and T
// units eV2 Bm-2
void Eloss() {
	vector<double> mass;
	vector<double> Q;
	double dw = 1e2;
	n = 1;
	for( double Bm = 1e-1; Bm <= 1e4; Bm*=10 ) {
		double total = 0;
		for( int j = wp.size()-1; j > 0; j-- ){
			total += 0.5*abs(r[j]-r[j-1])*( wp[j]*kIntg(Bm,j) + wp[j-1]*kIntg(Bm,j-1) )
					+ 0.5*abs(wp[j]-wp[j-1])*(wp[j]*T_solarIntg(wp[j],Bm) + wp[j-1]*T_solarIntg(wp[j-1],Bm));
			//total += 0.5*(w1-w2)*(w1*T_solarIntg(w1,Bm) + w2*T_solarIntg(w2,Bm));	// only T
			//total += 0.5*abs(r[j]-r[j-1])*( wp[j]*kIntg(Bm,j) + wp[j-1]*kIntg(Bm,j-1) );	// only L
		}
		for( double w = wp[0]; w < 2e4; w+=dw ){
			total += 0.5*dw*( (w+dw)*T_solarIntg(w+dw,Bm) + w*T_solarIntg(w,Bm) );
		}
		mass.push_back(Bm);
		Q.push_back(total);
		cout<<"Bm = "<<Bm<<endl;
	}
	// write to file
	string name = "data/Eloss_Bm_TL.dat";
	write2D( name , mass, Q );
}


// calculate total energy loss rate as a function of Bg
// both L and T
// units eV2 Bm-2
void Eloss_Bg() {
	vector<double> mass;
	vector<double> Q;
	double dw = 1e1;
	n = 1;
	double Bm = 1e2;
	double w1, w2, r1, r2 = 0;
	for( double Bg = 1e0; Bg <= 1e100; Bg*=10 ) {
		double total = 0;
		for( int j = wp.size()-1; j >= 0; j-- ){
			w1 = wp[j];
			if(w2 >= w1) { continue; }
			else{
			r1 = r[j];
			total += 0.5*(w1-w2)*( (w1+w2)*( kIntg_B(Bg, j) * abs((r2-r1)/(w2-w1)) )
								 + w1*T_solarIntg_B(w1,Bg) + w2*T_solarIntg_B(w2,Bg) );
			r2 = r[j];
			w2 = wp[j];
			}
		}
		for( double w = w1+dw; w < 2e4; w+=dw ){
			total += 0.5*dw*( (w+dw)*T_solarIntg_B(w+dw,Bg) + w*T_solarIntg_B(w,Bg) );
		}
		mass.push_back(Bg);
		Q.push_back(total);
		cout<<"Bg = "<<Bg<<endl;
	}
	// write to file
	string name = "data/Eloss_Bg--constB.dat";
	write2D( name , mass, Q );
}


// calculate energy loss as a function of n
// units eV2 Bg-2
void Eloss_n() {
	vector<double> nvec;
	vector<double> Q;
	double dw = 1e1;
	double Bm = 1e2;
	n = -10;
	double w1, w2, r1, r2 = 0;
	while ( n <= 100 ) {
		if( (n<0) && ((int)n%2 != 0) ) { continue; }
		double total = 0;
		for( int j = wp.size()-1; j >= 0; j-- ){
			w1 = wp[j];
			if(w2 >= w1) { continue; }
			else{
			r1 = r[j];
			total += 0.5*(w1-w2)*( (w1+w2)*( kIntg(Bm, j) * abs((r2-r1)/(w2-w1)) )
								 + w1*T_solarIntg(w1,Bm) + w2*T_solarIntg(w2,Bm) );
			if(isnan(total)) {cout<<"Bm = "<<Bm<<"	w = "<<w1<<"	j = "<<j<<endl;}
			r2 = r[j];
			w2 = wp[j];
			}
		}
		for( double w = w1+dw; w < 2e4; w+=dw ){
			total += 0.5*dw*( (w+dw)*T_solarIntg(w+dw,Bm) + w*T_solarIntg(w,Bm) );
			//if((int)(w) % (int)(1e3) == 0) { cout<<"w = "<<w/1e3<<"keV of 20keV"<<endl; }
		}
		nvec.push_back(n);
		Q.push_back(total);
		//if((int)(log10(Bm)) % 1 == 0) { cout<<"Bm = 1e"<<(int)(log10(Bm))<<" of 1e8"<<endl; }
		cout<<"n = "<<n<<endl;
		n+=4;
	}
	// write to file
	string name = "data/Eloss_n--1e2.dat";
	write2D( name, nvec, Q );
}


// calculate energy loss as a function of Lambda
// units eV2 Bg-2
void Eloss_Lambda() {
	vector<double> Evec;
	vector<double> Q;
	double dw = 1e1;
	double Bm = 1e6;
	n = 1;
	E = 1e-8;
	//double w1, w2, r1, r2 = 0;
	while ( E <= 1e1 ) {
		//if( (n<0) && ((int)n%2 != 0) ) { continue; }
		double total = 0;
		//for( int j = wp.size()-1; j >= 0; j-- ){
		//	w1 = wp[j];
		//	if(w2 >= w1) { continue; }
		//	else{
		//	r1 = r[j];
		//	total += 0.5*(w1-w2)*( (w1+w2)*( kIntg(Bm, j) * abs((r2-r1)/(w2-w1)) )
		//						 + w1*T_solarIntg(w1,Bm) + w2*T_solarIntg(w2,Bm) );
		//	if(isnan(total)) {cout<<"Bm = "<<Bm<<"	w = "<<w1<<"	j = "<<j<<endl;}
		//	r2 = r[j];
		//	w2 = wp[j];
		//	}
		//}
		for( double w = dw; w < 2e4; w+=dw ){
			total += 0.5*dw*( (w+dw)*T_solarIntg(w+dw,Bm) + w*T_solarIntg(w,Bm) );
			//if((int)(w) % (int)(1e3) == 0) { cout<<"w = "<<w/1e3<<"keV of 20keV"<<endl; }
		}
		Evec.push_back(E);
		Q.push_back(total);
		//if((int)(log10(Bm)) % 1 == 0) { cout<<"Bm = 1e"<<(int)(log10(Bm))<<" of 1e8"<<endl; }
		cout<<"E = "<<E<<endl;
		E*=1.1;
	}
	// write to file
	string name = "data/Eloss_Lambda_T--1e6.dat";
	write2D( name, Evec, Q );
}


// chameleon mass squared as a function of solar radius for given model parameters
// units eV2
void mass_profile() {
	vector<double> radius, mass;
	double Bm = 1e2;		// cham matter coupling
	n = 1;					// cham model n
	E = 2.4e-3;				// Lambda [eV]
	for( int c = 0; c < r.size(); c++ ) {
		radius.push_back(r[c]);
		mass.push_back(mCham2(c,Bm));
	}
	// write to file
	string name = "data/mass_profile_1e2.dat";
	write2D( name , radius, mass );
}


// get m as function of Lambda and B_m for given n
// m in eV
// down columns: Bm from 1e0 to 1e10
// right along rows: Lambda from 1e-10 to 1e0
// ie. dat[r,w]
void mass_contour() {
	vector<double> mass, BOut, LOut;
	double mB = 1.1;
	double mL = mB;
	n = 1;
	int c = 0;
	string name = "data/massregion_n1--1e3.dat";
	for( double Bm = 1e-1; Bm <= 1e18; Bm*=mB ) {
		BOut.push_back(Bm);
		for( double L = 1e-8; L <= 1e1; L*=mL ) {
			E = L;
			double total = sqrt( mCham2(0,Bm) );
			if(total >= 1e3) {total = 1e3;}
			mass.push_back(total);
			//cout<<total<<endl;
			if(c==0) { LOut.push_back(E); }
		}
		writeREST( name, mass, c );
		//if(int(10*log10(Bm))%10 == 0) { cout<<"Bm = "<<Bm<<" (of "<<1e10<<") written!"<<endl; }
		mass.clear();
		c++;
	}
	write2D("data/massregion_Bm.dat",BOut,BOut);
	write2D("data/massregion_Lambda.dat",LOut,LOut);
	cout<<"\acompleted it mate"<<endl;
}


// calculate differential particle flux spectrum by intg over solar volume
// LT coalescence and plasmon decay
// dN/dw, units Bg-2
void spectrum_lt() {
	vector<double> count, energy;
	Bm = 1e2;		// cham matter coupling (or fixed scalar mass)
	n = 1;					// cham model n
	double dw = 1e0;
	for( double w = dw; w < 2e4; w+=dw ){
		energy.push_back(w);					// eV
		count.push_back( solarIntg_lt(w,Bm) );		// Bg-2
		if((int)(w) % (int)(1e3) == 0) { cout<<"w = "<<w/1e3<<"keV of 20keV"<<endl; }
	}
	// write to file
	string name = "data/coalescence_lt_spectrum_1e2--uncapped.dat";
	write2D( name , energy, count );
}


void profile_ll() {
	vector<double> radius, rate;
	//double ms = 1e-3;		// scalar mass 1 eV
	//double ms2 = ms*ms;
	double Bm = 1e2;		// cham matter coupling
	n = 1;					// cham model n
	for( int c = 0; c < r.size(); c++ ) {
		radius.push_back(r[c]);
		rate.push_back( kIntg_ll(Bm, c) );
	}
	// write to file
	string name = "data/coalescence_ll_profile_1e2.dat";
	write2D( name , radius, rate );
}


// calculate differential particle flux spectrum by intg over solar volume
// dN/dw, units Bg-2
// (where dN is really dN/dt, sorry)
void spectrum_ll() {
	vector<double> count, energy;
	//double ms = 1e-3;		// scalar mass 1 eV
	//double ms2 = ms*ms;
	double Bm = 1e2;		// cham matter coupling
	n = 1;					// cham model n
	double w1, w2 = 0;
	double r1, r2 = rSolar;
	for( int j = wp.size()-1; j >=0; j-- ){
		w1 = wp[j];
		if(w2 > w1) { continue; }
		else{
		r1 = r[j];
		energy.push_back(2*w1);
		count.push_back( kIntg_ll(Bm, j) * abs((r2-r1)/(w2-w1))/2 );
		r2 = r[j];
		w2 = wp[j];
		}
	}
	// write to file
	string name = "data/coalescence_ll_spectrum_1e2--test.dat";
	write2D( name , energy, count );
}


// ll spectrum for omega intg
void spectrum_ll_omega() {
	vector<double> count, energy;
	Bm = 1e2;		// cham matter coupling (or fixed scalar mass)
	n = 1;					// cham model n
	double dw = 1e0;
	for( double w = dw; w < 2e4; w+=dw ){
		energy.push_back(w);					// eV
		count.push_back( (solarIntg_ll_omega(w+dw,Bm) - solarIntg_ll_omega(w,Bm))/dw );		// Bg-2
		//if((int)(w) % (int)(1e3) == 0) { cout<<"w = "<<w/1e3<<"keV of 20keV"<<endl; }
		//cout<<"w = "<<w<<"eV of 20keV"<<endl;
	}
	// write to file
	string name = "data/coalescence_ll_spectrum_omega.dat";
	write2D( name , energy, count );
}


///////////////////////////////////////////////////////////////////////////////////////////////////


int main() { 
	// convert Gaunt factor Theta to T in eV
	//for( int i = 1; i < 201; i++ ) { z1[0][i] = z1[0][i] * me; }
	//for( int i = 1; i < 201; i++ ) { z2[0][i] = z2[0][i] * me; }

	//spectrum('B');
	//spectrumL();
	//profile('B');
	Eloss_Lambda();
	//spectrum_ll();
	return 0;
}

