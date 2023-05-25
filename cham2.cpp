// Tom O'Shea 2023
// simple script to find chameleon limits for CAST just from tachocline
// some values taken from https://github.com/lucavisinelli/XENONCHAM/blob/main/code/ChameleonFlux.py

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <math.h> 
#include <filesystem>
#include <bits/stdc++.h>

using namespace std;

// physical constants
double hbar = 6.582119569e-16;  // [eV s]
double ls =  2.99792458e8;  // 1 lightsecond [m]
double hbarc = hbar * ls;   // m2eV [eV m]
double pi = 3.14159265358979323846;	// pi
double m_e = 510998.950;		// m_e [eV]
double m_p = 938272088;	// m_p [eV]
double a = (1. / 137.);		// alpha
double mpl = 2e27;   // planck mass [eV]

// conversion factors
double s2eV = (6.582119569e-16);	// Hz to eV (hbar)
double J2eV = (1. / 1.602176634e-19);	// Joules to eV (1 / e)
double m2eV = (1.973269804e-7);	// m-1 to eV
double K2eV = (8.617333262e-5);	// Kelvin to eV
double kg2eV = 5.609588604e35;	// from hbar/c2
double T2eV = 2e-16 * 1e18; // Tesla to eV2

// other constants
double z3 = 1.202056903159594;  // Riemann zeta(3) 
double Lam = 1e-3;    // cosmological constant [eV] 2.4e-3

// solar params
double R = 149.5978707e9;	// mean earth-sun distance [m]
double rSolar = 6.9598E+10;	// solar radius [m]
double xt = 0.7;    //radial location of tachocline [frac]
double Dx = 0.01;   // tachocline thickness [frac]
double rt = xt * rSolar;   // radial location of tachocline [m]
double Dr = Dx * rSolar;  // tachocline thickness [m]
double Bt = 30 * T2eV;    // tachocline B-field [eV2]
double rho = 891921.86598 * 1e12; // tachocline density [eV4]
double T = 0.1995759863 * 1e3;  // tachocline temperature [eV]
double ne = 1.06051586e+23 * 1e6 * pow(m2eV,3); // e- number density [eV3]
double mfp = 0.003;  // photon mean free path [m]
double ng = 1e21 * 1e4;  // tachocline photon flux [m-2 s-1]
double wp2 = 4*pi*a*ne/(m_e); // plasma frequency squared [eV2]

// CAST params
//double L = 9.26; // CAST length [m]
//double B = 9 * T2eV;    // CAST B-field [eV2]
//double phi = 0.12;   // CAST bg flux [m-2 s-1] (95% CL)


// write out 2 column datafile
void write2D( string name, vector<double> data1, vector<double> data2) {

	// delete old file if extant
	if ( remove(name.c_str()) == 0 ) {
		cout << "file " << name << " deleted." << endl;
		}
	
	// create new file
	fstream fout;
	fout.open(name, ios::out | ios::trunc);
	
	// get size of vectors
	int len1 = data1.size();
	int len2 = data2.size();
	
	// check vectors of equal length
	if ( len1 == len2 ) {
		
		// Read the input from vector
		for ( int c = 0; c < len1; c++ ) {
			// Insert the data to file
			fout << data1[c] << "	" << data2[c] << endl;
		}
		
		cout << "file " << name << " created succesfully." << endl;
	}
	
	// if vectors not equal - error
	else { cout << "ERROR - ensure vectors are of equal length" << endl; }
	
	fout.close();
	
}

// chameleon mass squared
double mCham( double w, double n, double Bm ){
    // calculate chameleon mass squared [eV2]
    double wd2 = ( (n+1) * rho / mpl ) 
    * pow( rho / (n * mpl * pow( Lam, (n+4) ) ) , ( 1 / (n+1) ) );

    double m2 = pow( Bm, (n+2)/(n+1) ) * wd2;

    return m2;
}

// symmetron mass squared (for high density approx) [eV2]
double mSym( double Bm ){
    return Bm * rho / pow(mpl,2);
}

// density dependant dilaton (Brax) for negligible Vc [eV2]
double mDil( double Bm ){
    return Bm * rho / pow(mpl,2);
}

// coherence length [eV-1]
double lom( double w, double n, double Bm ) {
    
    // get m2 from corresponding theory
    double m2 = mCham( w, n, Bm );

    // calculate coherence length (pos def)
    double item = 4 * w / sqrt( pow( m2 - wp2, 2 ));
    return item;
}

// function I(a)
double I( double a ) {
    return sqrt(pi/2) * ( sqrt(a + sqrt( pow(a,2) + 4 ) ) - sqrt( 2 * a ) );
}

// tachocline photon distribution [eV-1]
double pg( double w ) {
    return pow(w,2) / ( ( 2 * z3 * pow(T,3) ) * ( exp(w/T) - 1 ) );
}

/*
// phi from dw integral
double wIntegral( double n, double Bm ) {

    // integrate wrt w over CAST energies (0.5 - 15 keV) by trapezia
    double dw = 1e1;
    double item = 0.;
    for( double w = 5e2; w < 15e3; w+=dw) {

        double l = lw(w,n,Bm);
        double p1 = ng * pg(w) * sqrt(ls/l) * pow(rSolar,3) * pow( B * Bt * l * xt / ( 4 * mpl*mpl * R ), 2 ) * Dx * I(l/mfp) / mfp;
        l = lw(w+dw,n,Bm);
        double p2 = ng * pg(w+dw) * sqrt(ls/l) * pow(rSolar,3) * pow( B * Bt * l * xt / ( 4 * mpl*mpl * R ), 2 ) * Dx * I(l/mfp) / mfp;

        item += ( dw * (p1+p2) / 2 );
    }

    return item;
}
*/


// differential flux emitted from sun [m-2 s-1 eV-1]
double solarFlux( double w, double n, double Bm ) {

    double lw = lom(w,n,Bm);  // [eV-1]
    return ng * pg(w) * (Dr / mfp) * pow( Bt * lw / mpl*2, 2 ) * sqrt(ls/lw*hbarc) * I(lw*hbarc/mfp);
}


// detector phi from dw integral [m-2 s-1]
double wIntegral( double n, double Bm, double L, double B ) {

    // integrate wrt w over CAST energies (0.5 - 15 keV) by trapezia
    double dw = 1e0;
    double item = 0.;
    for( double w = 1e2; w < 1e4; w+=dw) {     // omega in eV
        item += ( dw * pow(rt/R,2) * pow(B*(L/hbarc)/(2*mpl),2) * (solarFlux( w+dw, n, Bm ) + solarFlux( w, n, Bm )) / 2 );
    }

    return item;
}

// calculate limits
void calc ( double n, double L, double B, string detector ) {
    
    // initialise vectors
    vector<double> BgVec;
    vector<double> BmVec;
    // scan over various Bm
    for ( double Bm = 1e-10; Bm < 1e6; Bm*=2 ) {

        //BgVec.push_back( pow( phi / wIntegral(n,Bm), 0.25) );
        BgVec.push_back(wIntegral(n, Bm, L, B));   // output coupling=1 flux
        BmVec.push_back(Bm);
    }

	// set path for writeout
	string path = "data/limits/" + detector;
	string ext = "-cham-flux.dat";
	write2D( path + to_string((int)n) + ext, BmVec, BgVec );
}


int main(){

    // set model parameter n
    double n = 1;

    // run over each

    // babyIAXO
    double L = 10;  // babyIAXO bore length [m]
    double B = 2 * T2eV;   // babyIAXO B-field [eV2]
    thread t1(calc, n, L, B, "babyIAXO");

    // baseline IAXO
    L = 2.5; B = 20;
    thread t2(calc, n, L, B, "baselineIAXO");

    // upgraded IAXO
    L = 3.5; B = 22;
    thread t3(calc, n, L, B, "upgradedIAXO");

    t1.join();
    t2.join();
    t3.join();

    return 0;
}

// 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