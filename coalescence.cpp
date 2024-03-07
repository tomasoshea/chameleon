// Tom O'Shea 2024

// Primakoff production of scalars in the sun
// V3 - the one using full ee scatter with full disp. rel. but no Raffelt screening=

#include "utils.h"

using namespace std;

// constants
double pi = 3.14159265359;
double alpha = 1/137.035999084;
double me = 510998.950;		// e- mass [eV]
double Mpl = 2e27;   		// planck mass [eV]
double E = 2e-3;    		// cham potential energy scale [eV] (2.4e-3 for cosmological const)
double n = 0;				// chameleon potential 1/phi^n
double rSolar = 6.957e8/m2eV;					// solar radius [eV-1]
double dSolar = 149.5978707e9/(1.973269804e-7);		// mean earth-sun distance [eV-1]
int nancount = 0;
int ucount = 0;

// solar params
vector<double> ne = read("data/ne.dat");		// electron number density [eV3]
vector<double> wp = read("data/wp.dat");		// electron number density [eV3]
vector<double> T = read("data/T.dat");			// solar temperature [eV]
vector<double> r = read("data/r.dat");			// radial distance [eV-1]
vector<double> rho = read("data/rho.dat");		// solar density [eV4]
vector<double> nH = read("data/nH.dat");		// H number density [eV3]
vector<double> nHe3 = read("data/nHe3.dat");	// He3 number density [eV3]
vector<double> nHe4 = read("data/nHe4.dat");	// He4 number density [eV3]
// get gaunt factors
vector<vector<double>> z1 = readGaunt("data/Z1.dat");	// gaunt factors for Z=1
vector<vector<double>> z2 = readGaunt("data/Z2.dat");	// gaunt factors for Z=2



// chameleon mass as a function of solar radius and model parameters
// cham model params n (phi-n potential), Bm (matter coupling)
// assume rho dominated by matter density
double mCham2( int c, double Bm ) {
	if(n<1) { cout<<"ERROR! n < 1"<<endl; return 0; }
	double E4n = pow(E,4+n);
	return n*(n+1)*E4n*pow( Bm*rho[c]/(n*Mpl*E4n), (n+2)/(n+1) );
}


// differential scalar production rate on earth dN/dr times Lambda2
// units eV2 Bm-2
double integrand_ll( int c, double Bm ) {
	double Tc = T[c];
	double rc = r[c];
	double nec = ne[c];
	if( Tc==0 ) { return 0; }				// solves weird behaviour when ne = T = 0
	double ms2 = mCham2(c,Bm);				// chameleon mass2 [eV2]
	double mg2 = 4*pi*alpha*nec/me;		// assume mg2 = wp2
	//double ms2 = Bm*Bm;					// fixed scalar mass2 [eV2]
	double K2 = 8*pi*alpha*nec/Tc;		// Debye screening scale ^2 [eV2]
	if( 2*mg2 <= ms2 ) { cout<<ms2-(2*mg2)<<endl; return 0;}
	return 1/(144*Mpl*Mpl*pi*pi) * pow(K2,3/2) * Tc*Tc * rc*rc * sqrt(4*mg2 - ms2);
}


// differential scalar production rate on earth d2N/dr/dw times Lambda2
// units eV Bm-2
double integrand_lt( int c, double Bm, double w ) {
	double Tc = T[c];
	double rc = r[c];
	if( Tc==0 ) { return 0; }				// solves weird behaviour when ne = T = 0
	double ms2 = mCham2(c,Bm);				// chameleon mass2 [eV2]
	double mg2 = 4*pi*alpha*ne[c]/me;		// assume mg2 = wp2
	if( w*w <= mg2 ) { return 0; }			// w > m_s requirement
	if( w*w <= ms2 ) { return 0; }			// w > m_g requirement
	if( ms2 <= mg2 ) { return 0; }
	//double ms2 = Bm*Bm;					// fixed scalar mass2 [eV2]
	double wt = w - sqrt(mg2);				// t-photon omega [eV]
	double kphi = sqrt(w*w - ms2);			// scalar wavenumber [eV]
	double kt = sqrt(w*(w - 2*sqrt(mg2)));	// t-photon wavenumber [eV]
	if(w - 2*sqrt(mg2) <= 0) { return 0; }
	//cout<<kt<<endl;
	return rc*rc/Mpl/Mpl * wt*wt * kphi * kt * Tc/(exp(wt/Tc) - 1);
}

// integral over solar volume, for a given scalar mass and energy
// returns dPhi/dw Bg-2 [eV2]
double solarIntg_lt( double w, double Bm ) {
	double total = 0;
	for( int c = 0; c < r.size() - 1; c++ ) {
		total += 0.5 * (r[c+1] - r[c]) * (integrand_lt(c+1, Bm, w) + integrand_lt(c, Bm, w));
	}
	return total;
}


// calculate differential particle flux spectrum by intg over solar volume
// dN/dw, units Bg-2
// (where dN is really dN/dt, sorry)
void spectrum_lt() {
	vector<double> count, energy;
	double Bm = 1e4;		// cham matter coupling (or fixed scalar mass)
	n = 1;					// cham model n
	double dw = 1e0;
	for( double w = dw; w < 2e4; w+=dw ){
		energy.push_back(w);					// eV
		count.push_back( solarIntg_lt(w,Bm) );		// Bg-2
		//cout<<solarIntg_lt(w,Bm)<<endl;
		if((int)(w) % (int)(1e3) == 0) { cout<<"w = "<<w/1e3<<"keV of 20keV"<<endl; }
	}
	// write to file
	string name = "data/coalescence_lt_spectrum_1e4.dat";
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
		rate.push_back( integrand_ll(Bm, c) );
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
		energy.push_back(w1);
		count.push_back( integrand_ll(Bm, j) * abs((r2-r1)/(w2-w1)) );
		//cout << integrand_ll(Bm, j) << endl;
		r2 = r[j];
		w2 = wp[j];
		}
	}
	// write to file
	string name = "data/coalescence_ll_spectrum_1e2.dat";
	write2D( name , energy, count );
}




int main() { 
	//L_profile();
	spectrum_ll();
	//Eloss();
	//contour();
	return 0;
	}
