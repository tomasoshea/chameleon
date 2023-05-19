// Tom O'Shea 2023
// script to find chameleon limits for IAXO over whole sun

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

// solar params
double R = 149.5978707e9;	// mean earth-sun distance [m]
double rSolar = 6.957e8;	// solar radius [m]
double B0 = 3e3;	// radiative zone max B [T]
double B1 = 50;	// tachocline max B [T]
double B2 = 3;	// outer region max B [T]
double r0 = 0.712;	// [R0]
double r1 = 0.732;	// [R0]
double d1 = 0.02;	// [R0]
double r2 = 0.96;	// [R0]
double d2 = 0.035;	// [R0]


// other constants
double z3 = 1.202056903159594;  // Riemann zeta(3) 
double Lam = 1e-3;    // cosmological constant [eV] 2.4e-3

// CAST params
//double L = 9.26; // CAST length [m]
//double B = 9 * T2eV;    // CAST B-field [eV2]
//double phi = 0.12;   // CAST bg flux [m-2 s-1] (95% CL)

// babyIAXO
double L = 10;  // babyIAXO bore length [m]
double B = 2 * T2eV;   // babyIAXO B-field [eV2]
double phi = 2.89557e-06;   // babyIAXO bg flux [m-2 s-1] (95% CL)



// read in datafiles in csv form
vector<double> read( string name ) {

	//cout << "Reading file " << name << "..." << endl;

	// open file defined in argument
	fstream file;
	file.open(name,ios::in);
	
	char delim('\n');	// define delimiter for file parsing
	
	if (file.is_open()){   // checking whether the file is open
		string temp;	// define temporary storage string
		vector<double> row;	// define vector to store input values and return
		
		while(getline(file, temp, delim)){  // read data from file object and put it into string
			double item = stod(temp);	// convert string to double
			row.push_back(item);	// append double to vector
		}
		
	file.close();   // close the file object.
	
	return row;	// returns vector of values
	
	}
	else {
		cout << "file " << name << " doesnt exist..." << endl;
		return {69.};
	}
}

// read in gaunt factors from matlab matrix files
vector<vector<double>> readGaunt( string name ) {

	//cout << "Reading file " << name << "..." << endl;

	// open file defined in argument
	fstream file;
	file.open(name,ios::in);
	
	char delim(' ');	// define delimiter for file parsing
	
	if (file.is_open()){   // checking whether the file is open
		string temp;	// define temporary storage string
		vector<vector<double>> g2D; // 2D matrix to store rows
		vector<double> row;	// define vector to store input values and return
		vector<double> x;	// define vector to store input values and return
		
		int c = 0;	// define counter for elements
		int line = 0;	// counter for lines
		
		while( !file.eof() ) {  // loop until end of file
		
			getline(file, temp, delim);	// get one row
			
			// check data is not empty
			if( temp == "\n" ) { continue; }
			
			double item = stod(temp);	// convert string to double

			// add a zero on first line
			if( c == 0 ) { row.push_back(0.); }
			
			row.push_back(item);	// append double to vector
			c++;
			
			// when row is full append to 2D vector and wipe row				
			if( row.size() == 201 ) {
			
				g2D.push_back(row);
				row.clear();
				
				line++;
			}

			if( line == 501 ) { break; }

		}
		
	file.close();   // close the file object.
	
	return g2D;	// returns vector of values
	}
	
	else{ cout << "file " << name << " doesnt exist..." << endl;
			vector<vector<double>> err = {{69.}} ; return err ; }
}


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

// B field in solar radiative zone [T]
double radiativeZone( double r ) {

	double l = (10 * r0) + 1;
	double K = (1+l) * pow( 1 + pow(l,-1) , l ) * B0;
	
	double part1 = pow( r/r0 , 2 );
	double part21 = 1 - pow( r/r0 , 2 );
	double part2 = pow( part21, l );
	
	double a = K * part1 * part2;
	return a;	
}


// B field in tachocline and outer regions [T]
double tachocline( double r, double rmax, double d, double B ) {

	double a = B * ( 1 - pow( (r - rmax)/d , 2 ) );
	return a;
}

// T plasmon production rate
double Gamma( double w, double number, double T, double nH, double nHe4, double nHe3, double g1, double g2 ) {

	double p1 = 16 * pow(pi,2) * pow(a,3);
	double p2 = 3 * pow(m_e,2) * pow(w,3);
	double p3 = 2 * pi * m_e * pow(number,2) / (3 * T);
	double p4 = 1 - exp(- w / T);
	double p5 = 8 * pi * pow(a,2) * number / (3 * pow(m_e,2) );
	
	// sum of ion densities
	double ions = (nH * g1) + g2 * ( (4 * nHe4) + (4 * nHe3) );

	double item = p1 * pow(p2, -1) * pow(p3, 0.5) * p4 * ions;
	
	item += p5;

	return item;
}

// calculate photon mfp from gamma
double mfpath( double w, double ne, double T, double nH,
			double nHe4, double nHe3,
			vector<vector<double>> z1, vector<vector<double>> z2 ){
		
	// select g(w, T) value from matrix
	int indexT1;
	int indexT2;
	int indexX1;
	int indexX2;
	
	for( int i = 1; i < 200; i++ ) {
		if( z1[0][i] < T and z1[0][i+1] > T ) { indexT1 = i; }
		if( z2[0][i] < T and z2[0][i+1] > T ) { indexT2 = i; }
	}
	
	for( int i = 1; i < 500; i++ ) {
		if( (z1[i][0] * T) < w and (z1[i+1][0] * T) > w ) { indexX1 = i; }
		if( (z2[i][0] * T) < w and (z2[i+1][0] * T) > w ) { indexX2 = i; }
	}
	
	double g1 = z1[ indexT1 ][ indexX1 ];
	double g2 = z2[ indexT2 ][ indexX2 ];

	double G = Gamma( w, ne, T, nH, nHe4, nHe3, g1, g2 );	// eV
	return ( m2eV / G );	// [m]
}


// chameleon mass squared
double mCham( double w, double n, double Bm, double rho ){
    // calculate chameleon mass squared [eV2]
    double wd2 = ( (n+1) * rho / mpl ) 
    * pow( rho / (n * mpl * pow( Lam, (n+4) ) ) , ( 1 / (n+1) ) );

    double m2 = pow( Bm, (n+2)/(n+1) ) * wd2;

    return m2;
}

// symmetron mass squared (for high density approx) [eV2]
double mSym( double Bm, double rho ){
    return Bm * rho / pow(mpl,2);
}

// density dependant dilaton (Brax) for negligible Vc [eV2]
double mDil( double Bm, double rho ){
    return Bm * rho / pow(mpl,2);
}

// coherence length [eV-1]
double lom( double w, double n, double Bm, double rho, double wp ) {
    
    // get m2 from corresponding theory
    double m2 = mCham( w, n, Bm, rho );

    // calculate coherence length (pos def)
    double item = 4 * w / sqrt( pow( m2 - pow(wp,2) , 2 ));
    return item;
}

// function I(a)
double I( double a ) {
    return sqrt(pi/2) * ( sqrt(a + sqrt( pow(a,2) + 4 ) ) - sqrt( 2 * a ) );
}

// tachocline photon distribution times flux [eV2]
double npg( double w, double T ) {
    return pow(w,2) / ( pow(pi,2) * ( exp(w/T) - 1 ) );
}


// dr integral of solar flux []
double solarFlux( double w, double n, double Bm, vector<double> ne,
		 vector<double> rho,
        vector<double> wp, vector<double> r, vector<double> T,
		vector<double> nH, vector<double> nHe4, vector<double> nHe3,
		vector<vector<double>> z1, vector<vector<double>> z2 ) {
    
	// loop over whole sun with B-field
	int len = r.size();
	double total = 0;
    for ( int c = 0; c < len - 1; c+=10 ) {

		// get mfp
		double mfp1 = mfpath( w, ne[c], T[c], nH[c], nHe4[c],
							nHe3[c], z1, z2 );
		double mfp2 = mfpath( w, ne[c+1], T[c+1], nH[c+1], nHe4[c+1],
							nHe3[c+1], z1, z2 );

		// set values of solar B-field
		double B;
		double rFrac = r[c] / rSolar;
		// impoved model
		if ( rFrac <= r0 ) { B = radiativeZone( rFrac ); 
		if ( B > 1e10 ) { cout << B << endl; } }	// [T]
		else if ( rFrac > (r1 - d1) and rFrac < (r1 + d1)) {
			B = tachocline( rFrac, r1, d1, B1 );
			if ( B > 1e10 ) { cout << B << endl; } }	// [T]
		else if ( rFrac > (r2 - d2) and rFrac < (r2 + d2) ) {
			B = tachocline( rFrac, r2, d2, B2 ); 
			if ( B > 1e10 ) { cout << B << endl; } }	// [T]
		cout<< B << endl;
		
		
		B *= T2eV;	// [eV2]

		
		if( B == 0 ) { continue; }
    	double lw1 = lom( w, n, Bm, rho[c], wp[c] );  // [eV-1]
		double lw2 = lom( w, n, Bm, rho[c+1], wp[c+1] );  // [eV-1]
		double dr = r[c+1] - r[c];	// [m]
		total += ( ( ( npg(w,T[c]) * (1 / mfp1) * pow( B * lw1 / (2*mpl), 2 )
             * sqrt(ls/lw1*hbarc) * I(lw1*hbarc/mfp1) * pow(r[c]/R,2)  )
			 + ( npg(w,T[c+1]) * (1 / mfp2) * pow( B * lw2 / (2*mpl), 2 )
             * sqrt(ls/lw2*hbarc) * I(lw2*hbarc/mfp2) * pow(r[c+1]/R,2) ) )
			 * pow(B*(L/hbarc)/(2*mpl),2) * dr / 2 );
	}
	return total;
}


// detector phi from dw integral [m-2 s-1]
double wIntegral( double n, double Bm, vector<double> ne, vector<double> rho,
		vector<double> wp, vector<double> r,
		vector<double> T,
		vector<double> nH, vector<double> nHe4, vector<double> nHe3,
		vector<vector<double>> z1, vector<vector<double>> z2 ) {

    // integrate wrt w over CAST energies (0.5 - 15 keV) by trapezia
    double dw = 1e2;
    double item = 0.;
    for( double w = 1e2; w < 1e4; w+=dw) {     // omega in eV
        item += ( dw * (solarFlux( w+dw, n, Bm, ne, rho, wp, r, T,
						nH, nHe4, nHe3, z1, z2 )
			+ solarFlux( w, n, Bm, ne, rho, wp, r, T, nH,
						nHe4, nHe3, z1, z2 )) / 2 );
    }
    return item;
}

// calculate limits
int main(){

    vector<double> raw = read("data/r.dat");	// sun radial distance [eV-1]
	vector<double> T = read("data/T.dat");	// solar temperature [eV]
	vector<double> rho = read("data/rho.dat");	// solar density [eV4]
	vector<double> wp = read("data/wp.dat");	// plasma freq [eV]
	vector<double> ne = read("data/ne.dat");	// electron number density [eV3]
	vector<double> nH = read("data/nH.dat");	// H ion density [eV3]
	vector<double> nHe4 = read("data/nHe4.dat");	// He4 ion density [eV3]
	vector<double> nHe3 = read("data/nHe3.dat");	// He3 ion density [eV3]
	
	// get gaunt factors
	vector<vector<double>> z1 = readGaunt("data/Z1.dat");	// gaunt factors for Z=1
	vector<vector<double>> z2 = readGaunt("data/Z2.dat");	// gaunt factors for Z=2

	// convert Gaunt factor Theta to T in eV
	for( int i = 1; i < 201; i++ ) { z1[0][i] = z1[0][i] * m_e; }
	for( int i = 1; i < 201; i++ ) { z2[0][i] = z2[0][i] * m_e; }

	vector<double> r;
	for ( double i : raw ) { r.push_back( i * m2eV ); }	// [m]

    // set model parameter n
    double n = 1;
    // initialise vectors
    vector<double> BgVec;
    vector<double> BmVec;
    // scan over various Bm
    for ( double Bm = 1e0; Bm < 1e6; Bm*=2 ) {

        BgVec.push_back( pow( phi / wIntegral(n, Bm, ne, rho, wp, r,
									T, nH, nHe4, nHe3, z1, z2), 0.25) );
        BmVec.push_back(Bm);
		cout << log10((Bm+1)/1e6) << "\% complete..." << endl;
	// << "\33[A\33[2K\r"
    }

	// set path for writeout
	string path = "data/limits/babyIAXO";
	string ext = "-full-cham.dat";
	write2D( path + to_string((int)n) + ext, BmVec, BgVec );
    //write2D( path + ext, BmVec, BgVec );
    return 0;
}
