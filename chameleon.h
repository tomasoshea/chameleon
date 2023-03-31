// Tom O'Shea 2023

// chameleon header

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <math.h> 
#include <filesystem>
#include <bits/stdc++.h>

using namespace std;

// constants
double pi = 3.14159265358979323846;	// pi
double m_e = 510998.950;		// m_e [eV]
double m_p = 938272088;	// m_p [eV]
double a = (1. / 137.);		// alpha
double mpl = 1.2e29;   // planck mass [eV]

// solar params
double R_raw = 149.5978707e9;	// mean earth-sun distance [m]
double rSolar_raw = 6.9598E+10;	// solar radius [m]
double B0 = 3e3;	// radiative zone max B [T]
double B1 = 50;	// tachocline max B [T]
double B2 = 3;	// outer region max B [T]
double r0 = 0.712;	// [R0]
double r1 = 0.732;	// [R0]
double d1 = 0.02;	// [R0]
double r2 = 0.96;	// [R0]
double d2 = 0.035;	// [R0]
double mfp_raw = 0.1;	// approximate photon tachocline mean free path [m]


// conversion factors
double s2eV = (6.582119569e-16);	// Hz to eV
double J2eV = (1. / 1.602176634e-19);	// Joules to eV (1 / e)
double m2eV = (1.973269804e-7);	// m-1 to eV
double K2eV = (8.617333262e-5);	// Kelvin to eV
double kg2eV = 5.609588604e35;	// from hbar/c2

// converted values
double R = R_raw / m2eV;	// eV-1
double rSolar = rSolar_raw / m2eV;	// eV-1
double mfp = mfp_raw / m2eV;	// eV-1

// utility constants
double wRange = 1e3;	// range of w integral

// model params
//double n = 1;  // n=4 chameleon
double L = 0.1e-3;    // Lambda [eV]
double Bm = 1e6;	// matter mixing



///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// READ & WRITE ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


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
	else{ cout << "some error i guess..." << endl ; vector<double> err = {69.} ; return err ; }
}


// write out simple datafile
void write( string name, vector<double> data ) {

	// delete old file if extant
	if ( remove(name.c_str()) == 0 ) {
		cout << "file " << name << " deleted." << endl;
		}
	
	// file pointer
	fstream fout;

	// creates new csv file
	fout.open(name, ios::out | ios::trunc);

	// Read the input from vector
	for ( double item : data ) {
		// Insert the data to file
		fout << item << endl;
	}
	cout << "file " << name << " created succesfully." << endl;
	
	fout.close();
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
	
	// if vectors not equal return error
	else { cout << "ERROR - ensure vectors are of equal length" << endl; }
	
	fout.close();
	
}


// merge and write out integration output
void merge( string name ) {
	
	// read in data
	vector<double> m1 = read("data/" + name + "-mass-suppressed.dat");
	vector<double> m2 = read("data/" + name + "-mass-resonance.dat");
	vector<double> m3 = read("data/" + name + "-mass-unsuppressed.dat");
	
	vector<double> chi1 = read("data/" + name + "-chi-suppressed.dat");
	vector<double> chi2 = read("data/" + name + "-chi-resonance.dat");
	vector<double> chi3 = read("data/" + name + "-chi-unsuppressed.dat");
	
	// merge
	vector<double> m;
	m.insert( m.end(), m1.begin(), m1.end() );
	m.insert( m.end(), m2.begin(), m2.end() );
	m.insert( m.end(), m3.begin(), m3.end() );
	
	vector<double> chi;
	chi.insert( chi.end(), chi1.begin(), chi1.end() );
	chi.insert( chi.end(), chi2.begin(), chi2.end() );
	chi.insert( chi.end(), chi3.begin(), chi3.end() );
	
	// write out
	string path = "data/limits/";
	string ext = ".dat";
	write2D( path + name + ext, m, chi );	

}


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// GENERAL FUNCS ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


// complex sqrt of (a + ib), starting points a0, b0 and results stored at pointa, pointb
void csqrt( double a, double b, double* pointa, double* pointb ) {

	double x,y;
	double e = 1;
	if( a == 0 ) { x = 1; }
	else if ( a == 1 ) { x = 0.5; }
	else{ x = sqrt( abs(a) ); }
	if( b == 0 ) { y = 1; }
	else if ( b == 1 ) { y = 0.5; }
	else{ y = sqrt( abs(b) ); }
	int n = 0;
	double TOL = 1e-20;

	while( e > TOL ) {

		e = abs(  0.5 * ( x + (a*x + b*y) / ( pow(x,2) + pow(y,2) ) ) - x );
		//cout << e << endl;
		x = 0.5 * ( x + (a*x + b*y) / ( pow(x,2) + pow(y,2) ) );
		y = 0.5 * ( y + (b*x - a*y) / ( pow(x,2) + pow(y,2) ) );
		n++;
	}

	if( abs(x) < TOL ) { *pointa = 0; }
	else{ *pointa = x; }
	if( abs(y) < TOL ) { *pointb = 0; }
	else{ *pointb = y; }
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


// calculate m_eff^2
double meff( double Bm, double n, double rho, double wp ){

	// set constants
	//double Bg = pow(10,10.29);
	double B = 30 * 2e2;	// 30 T in eV2
	double Bg = pow(10,8);
	Bm = pow(10,8);
	n = 4;

	// add EM cpt to rho
	//cout << "rho1: " << rho << endl;
	rho += ( (Bg / Bm) * (B * B / 2) );
	//cout << "rho2: " << rho << endl;
	double wp2 = 4*pi*a*rho/(m_e*m_p);
	//cout << "plasma freq part: " << pow(wp,2) / rho << endl;

	//cout << "density part: " << ((n+1)/mpl)*pow( rho/(n*mpl*pow(L,n+4)) , 1/(n+1)) << endl;


    // compute omega_rho^2
    double wd = ( (n+1) * rho / mpl ) * pow( rho / (n * mpl * pow( L, (n+4) ) ) , ( 1 / (n+1) ) );
    //cout<<wd<<endl;
	//cout << pow(wp,2) << endl;

    // compute m_eff^2
    double item = pow( Bm, (n+2)/(n+1) ) * wd - pow(wp,2);

	return item;
    //return sqrt(item);
}

// CAST prob p
double pCAST( double w, double n, double rho, double wp ){

	// define m^2_eff and l_w
	double m2 = meff( Bm, n, rho, wp );
	double lw = 4 * w / m2;

	// set constants
	double Bg = pow(10,10.29);
	double B = 30 * 2e2;	// 30 T in eV2
	double l = 3 / m2eV;	// 3 m in eV

	double item = pow( Bg * B * lw / (2*mpl) , 2 ) * pow( sin(l/lw) , 2 );
	return item;
}


// plot conversion prob
double Pcham( double w, double n, double rho, double wp ){

	// B field in eV2
	double B = 30 * 2e2;	// 30 T in eV2
	double L = 0.1 * rSolar;	// tachocline length in eV

	// set B_gamma
	double Bg = pow(10,10.29);

	double m2 = meff( Bm, n, rho, wp );
	double theta = 0.5 * atan( 2 * w * B * Bg / ( m2 * mpl ) );
	cout << theta << endl;
	//double theta = w * B * Bg / ( m2 * mpl )
	double Delta = m2 * L / ( 4 * w );

	double P = pow( sin(2*theta) * sin( Delta / cos(2*theta) ) , 2 );
	return P;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// IAXO FLUX ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


// integrate over solar radius and E range
void integral( double Bm, double n, vector<double> rho, vector<double> wp, 
				vector<double> r, double* address ) {
	
	double len = r.size();
	double item = 0;

	// E integral
	for( double w = 1e3, w < 1e5, w = w + 1e3 ) {
	// solar radius integral
	for ( int c = 0, c < len, c++ ){

		double x = r[c] / rSolar;

		if ( x < (r1 - d1) or x > (r1 + d1) ) { continue; }

		double lw = 4 * w / meff( Bm, n, rho[c], wp[c] );
		double B = tachocline( x, r1, d1, B1 ) * 2e-16;	// eV2
		double temp = pow(x*B,2) * pow(lw,1.5) / mfp;
		item += temp;
	}
	}
	
	*address = item;
}

// calculate and output limits
void chameleon( double n, vector<double> rho, vector<double> wp, vector<double> r,
				double Biaxo, double L, double phi, string name ) {

	// define vectors
	vector<double> vecBg vecBm;

	// set path for writeout
	string path = "data/limits/";
	string ext = ".dat";

	// calculate for various Bm values
	for ( double Bm = 1e0, Bm < 1e6, Bm *= 2 ) {

		double* I;
		integral( Bm, n, hro, wp, r, I );
		double J = (*I) * C * pow(rSolar,3) * pow(L*Biaxo/R,2) * pow(2*mpl,-4);
		vecBg.push_back( pow(phi/J, 0.25) );
		vecBm.push_back( Bm );
	}

	write2D( path + name + ext, vegBm, vecBg );
}


///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// SPECTRA ETC ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


// function to output E spectrum for T-plasmon
void plotMeff( double n ) {
	
    // read in solar params
	vector<double> r = read("data/rFrac.dat");	// sun radial distance [eV-1]
	vector<double> rho = read("data/rho.dat");	// density [eV3]
	vector<double> wp = read("data/wp.dat");	// plasma frequency [eV]

    // set Bm
    //double Bm = 1e6;

    // set up vectors
    vector<double> rvec;
    vector<double> mvec;

    int len = wp.size();
	for ( int c = 0; c < len; c++ ) {

		// restrict to tachocline
		//if( r[c] < 0.6 ) { continue; }
		//else if( r[c] > 0.8 ) { continue; }

		// compute meff for each value of r
		double m = meff(Bm, n, rho[c], wp[c]);
		rvec.push_back(r[c]);
		//mvec.push_back( sqrt(m) );
		mvec.push_back(m);

	}

    //int Bint = (int)log10(Bm);
	//string name = "data/meff" + to_string( Bint ) + ".dat";
	string name = "data/meff" + to_string( (int)n ) + ".dat";

	write2D( name , rvec, mvec );
	cout << "length of r vector: " << rvec.size() << endl;
	cout << "length of m vector: " << mvec.size() << "\n" << endl;
}

void plotSpectrum(){

	// read in solar params
	vector<double> rho = read("data/rho.dat");	// electron number density [eV3]
	vector<double> wp = read("data/wp.dat");	// plasma frequency [eV]

	// set n
	double n = 1;

	// set up vectors
	vector<double> E;
	vector<double> P;

	for( double w = 1; w < 2e4; w += 1 ){

		//double prob = pCAST( w, n, rho[1400], wp[1400] );
		double prob = Pcham( w, n, rho[0], wp[0] );
		P.push_back(prob);
		E.push_back(w);
	}
	
	string name = "data/Espectrum" + to_string( (int)n ) + ".dat";

	write2D( name , E, P );
	cout << "length of E vector: " << E.size() << endl;
	cout << "length of P vector: " << P.size() << "\n" << endl;
}
