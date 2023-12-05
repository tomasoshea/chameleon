// Tom O'Shea 2023

// Primakoff production of scalars in the sun

#include "utils.h"

using namespace std;

// constants
double pi = 3.14159265359;
double alpha = 1/137.035999084;
double me = 510998.950;		// e- mass [eV]
double Mpl = 2e27;   		// planck mass [eV]
double E = 2e-3;    		// cham potential energy scale [eV] (2.4e-3 for cosmological const)
double n = 1;				// chameleon potential 1/phi^n

// solar params
vector<double> ne = read("data/ne.dat");		// electron number density [eV3]
vector<double> T = read("data/T.dat");			// solar temperature [eV]
vector<double> r = read("data/r.dat");			// radial distance [eV-1]
vector<double> rho = read("data/rho.dat");		// solar density [eV4]


// chameleon mass as a function of solar radius and model parameters
// cham model params n (phi-n potential), Bm (matter coupling)
// assume rho dominated by matter density
double mCham2( int c, double Bm ) {
	double E4n = pow(E,4+n);
	return n*(n+1)*E4n*pow( Bm*rho[c]/(n*Mpl*E4n), (n+2)/(n+1) );
}


// differential particle emssion rate d2N/dr/dw times Lambda2
// units Lambda2
double integrand( int c, double Bm, double w ) {
	double mg2 = 4*pi*alpha*ne[c]/me;	// assume mg2 = wp2
	//double ms2 = mCham2(c,Bm);			// chameleon mass2 [eV2]
	double ms2 = 1e1;					// fixed scalar mass2 [eV2]
	if( w*w <= mg2 ) { return 0; }
	if( w*w <= ms2 ) { return 0; }
	double nt = 2*ne[c];
	double K2 = 4*pi*alpha*nt/T[c];
	double l = (K2*K2 - 2*K2*(mg2 + ms2) + pow(mg2 - ms2,2))*pow(T[c],-4);
	double ln = log(pow(2*w/T[c]/T[c],2)*K2 + l);

	if( isnan(ln) ) { return 0; }	// get nan for T = 0
	if( ln < 0 ) { return 0; }	//
	return 2 * alpha / pi * nt * r[c]*r[c] * (w*w/(exp(w/T[c]) - 1)) * ln;
}


// integral over solar volume, for a given scalar mass and energy
double solarIntg( double w, double Bm ) {
	double total = 0;
	for( int c = 0; c < r.size(); c++ ) {
		total += 0.5 * (r[c+1] - r[c]) * (integrand(c+1, Bm, w) + integrand(c, Bm, w));
	}
	return total;
}

// interal over solar volume of interand times energy, for energy loss calc
double solarIntgw( double w, double Bm ) {
	double total = 0;
	for( int c = 0; c < r.size(); c++ ) {
		total += 0.5 * (r[c+1] - r[c]) * w * (integrand(c+1, Bm, w) + integrand(c, Bm, w));
	}
	return total;
}

// integral to get flux spectrum rather than number (divide by 4 pi r2)
double solarIntgr( double w, double Bm ) {
	double total = 0;
	for( int c = 0; c < r.size(); c++ ) {
		total += 0.5 * (r[c+1] - r[c]) * w * (integrand(c+1, Bm, w) + integrand(c, Bm, w)) / (4*pi*r[c]*r[c]);
	}
	return total;
}

// calculate emission rate profile over solar radius (dN/dr)
// integrated over relevant energies up to 20 keV
void profile() {
	vector<double> radius;
	vector<double> rate;
	double Bm = 1e6;
	double mw = 1.1;
	double dw = 1;		// [eV]
	for( int c = 0; c < r.size(); c++ ) {
		double total = 0;
		for( double w = dw; w < 2e4; w+=dw ) {
			//total += 0.5*w*(mw-1)*( integrand(c,Bm,w*mw) + integrand(c,Bm,w) );
			total += 0.5*dw*( integrand(c,Bm,w*mw) + integrand(c,Bm,w) );
		}
		radius.push_back(r[c]);
		rate.push_back(total);
		//double w = 1e3;
		//rate.push_back(integrand(c,Bm,w));
	}
	// write to file
	string name = "data/primakoff_profile_1eV.dat";
	write2D( name , radius, rate );
}

// calculate differential particle flux spectrum by intg over solar volume
void spectrum() {
	vector<double> count, energy;
	//double ms = 1e-3;		// scalar mass 1 eV
	//double ms2 = ms*ms;
	double Bm = 1e7;		// cham matter coupling
	for( double w = 1; w < 2e4; w+=1 ){
		energy.push_back(w);
		count.push_back( solarIntg(w,Bm) );
	}
	// write to file
	string name = "data/primakoff_spectrum_cham_1e7.dat";
	write2D( name , energy, count );
}


// calculate energy loss as a function of m
void Eloss() {
	vector<double> mass;
	vector<double> Q;
	double mw = 1.1;
	for( double ms = 1e-6; ms < 1e6; ms*=1.1 ) {
		double ms2 = ms*ms;
		double total = 0;
		for( double w = ms; w < 1e6; w*=mw ) {
			total += 0.5*w*(mw-1)*( solarIntgw(w*mw,ms2) + solarIntgw(w,ms2) );
		}
		mass.push_back(ms);
		Q.push_back(total);
	}
	// write to file
	string name = "data/primakoff_Eloss.dat";
	write2D( name , mass, Q );
}


int main() { 
	spectrum();
	return 0;
	}