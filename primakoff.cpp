// Tom O'Shea 2023

// Primakoff production of scalars in the sun

#include "utils.h"

using namespace std;

// constants
double pi = 3.14159265359;
double alpha = 1/137.035999084;
double me = 510998.950;		// e- mass [eV]

// solar params
vector<double> ne = read("data/ne.dat");	// electron number density [eV3]
vector<double> T = read("data/T.dat");	// solar temperature [eV]
vector<double> r = read("data/r.dat");	// radial distance [eV-1]


// differential particle emssion rate dN/dr/dw times Lambda2
// units Lambda2
double integrand( int c, double ms2, double w ) {
	double mg2 = 4*pi*alpha*ne[c]/me;
	if( w*w <= mg2 ) { return 0; }
	double nt = 2*ne[c];
	double K2 = 4*pi*alpha*nt/T[c];
	double l = (K2*K2 - 2*K2*(mg2 + ms2) + pow(mg2 - ms2,2))*pow(T[c],-4);
	double ln = log(pow(2*w/T[c]/T[c],2)*K2 + l);

	if( isnan(ln) ) { return 0; }	// get nan for T = 0
	if( ln < 0 ) { return 0; }	//
	return alpha / pi * nt * r[c]*r[c] * (w/(exp(w/T[c]) - 1)) * ln;
}


// integral over solar volume, for given w
double solarIntg( double w, double ms2 ) {
	double total = 0;
	for( int c = 0; c < r.size(); c++ ) {
		total += 0.5 * (r[c+1] - r[c]) * (integrand(c+1, ms2, w) + integrand(c, ms2, w));
	}
	return total;
}

double solarIntgw( double w, double ms2 ) {
	double total = 0;
	for( int c = 0; c < r.size(); c++ ) {
		total += 0.5 * (r[c+1] - r[c]) * w * (integrand(c+1, ms2, w) + integrand(c, ms2, w));
	}
	return total;
}


// calculate differential particle flux spectrum by intg over solar volume
void spectrum() {
	vector<double> count, energy;
	double ms = 1e-3;		// scalar mass 1 eV
	double ms2 = ms*ms;
	for( double w = 1e-6; w < 1e5; w*=1.01 ){
		energy.push_back(w);
		if(w*w <= ms2) { count.push_back(0); continue; }
		else { count.push_back( solarIntg(w,ms2) ); }
	}
	// write to file
	string name = "data/primakoff_spectrum_1e-3.dat";
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
	Eloss();
	return 0;
	}