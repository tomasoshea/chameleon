// Tom O'Shea 2023

// script to find the upper limit on phi for a given number of observed events
// using a poisson distribution, for use in IAXO chameleon analysis

#include <bits/stdc++.h>

using namespace std;

double CL = 0.95;	// confidence level
double dE = 10;	// E range [keV]
int samplesize = 1e3;		// size of random sample

//string n = "10"; // model parameter

// conversion factors
double s2eV = (6.582119569e-16);	// Hz to eV
double J2eV = (1. / 1.602176634e-19);	// Joules to eV (1 / e)
double m2eV = (1.973269804e-7);	// m-1 to eV
double K2eV = (8.617333262e-5);	// Kelvin to eV
double kg2eV = 5.609588604e35;	// from hbar/c2

// read 2nd column from 2 column datafile
vector<double> loadtxt( string name, int col ) {

	//cout << "Reading file " << name << "..." << endl;

	// open file defined in argument
	fstream file;
	file.open(name,ios::in);
	
	char delim('	');	// define delimiter for file parsing (tab)
	
	if ( file.is_open() ){   // checking whether the file is open
		string temp;	// define temporary storage string
		vector<double> row1, row2;	// define vector to store input values and return
		vector<string> v;
		
		while( getline(file, temp) ){ v.push_back( temp ); }

		for ( string i : v ) {

			stringstream X(i);
			string temp;
			vector<string> vec;
			while ( getline(X, temp, delim ) ) { vec.push_back(temp); }
			row1.push_back( stod(vec[0]) );
			row2.push_back( stod(vec[1]) );
		}
		
	file.close();   // close the file object.
	
	// choose column
	if ( col == 0 ) { return row1; }
	else if ( col == 1 ) { return row2; }
	else { cout << "put 0 or 1 for columns" << endl; return {69.}; }	
	}
	else{ cout << "file " << name << " doesn't exist" << endl; return {69.}; }
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


double L( double x, double b, double s, double n ) {

	double item;

	if ( n == 0 ) { item = ( b + pow(x,4)*s ); }
	else if ( n == 1 ) { item = ( ( b + pow(x,4)*s ) - log( b + pow(x,4)*s ) - 1 ); }
	else { item = ( ( b + pow(x,4)*s ) - n * ( log( b + pow(x,4)*s ) + 1 - log(n) ) ); }

	return exp(-item);
}

// integral for getting CL
double integral( double b, double s, double n ) {

	double x = 1e7;
	double x2 = x;
	//double dx = x;
	//cout << x2 << endl;
	double mx = 1.001;
	//double total = 0.5 * x * ( exp( - f( x, b, s, n ) ) + exp( - f( 0, b, s, n ) ) );
	double total= 0.5 * x * ( L( x, b, s, n ) + L( 0, b, s, n ) );
	double norm = total;

	// normalise with intg to inf
	while(true) {
		double dx = x2 * (mx - 1);
		double L1 = L( x2, b, s, n );
		//double L2 = L( x2*mx, b, s, n );
		double L2 = L( x2+dx, b, s, n );
		if ( isnan(L1 + L2) ) { continue; }
		else { norm += 0.5 * dx * ( L1 + L2 ); }
		if ( L2 + L1 == 0 ) { break; }
		//x2 += dx;
		x2 *= mx;
	}

	// integrate up to CL
	while ( ( total / norm ) < CL ) {
		double dx = x * (mx - 1);
		double L1 = L( x, b, s, n );
		double L2 = L( x+dx, b, s, n );
		//double L2 = L( x*mx, b, s, n );
		if ( isnan(L1 + L2) ) { continue; }
		else { total += 0.5 * dx * ( L1 + L2 ); }
		//cout << total * dx / norm << endl;
		x *= mx;
		//x2 += dx;
	}

	return x;
}


void chis( int detector, int nint ) {

	// initialise parameters
	double A, phiBg, a, t, effD, effO, effT, len;
	string name;
	vector<double> flux, m;
	string ext = "-cham-flux.dat";
	string path = "data/limits/";
	string nModel = to_string(nint);

	// choose detector

	if ( detector==0 ) {
		// babyIAXO parameters
		name="babyIAXO" + nModel;
		A = 0.77;	// detector area [m2]
		phiBg = 1e-7 * 1e4 * dE;	// background flux [m-2 s-1]
		a = 0.6 * 1e-4;	// XRay detection area [m2]
		t = 1.5 * 365.25 * 24 * 3600;	// 1.5 years
		effD = 0.7;	// detectior efficiency
		effO = 0.35;	// optical efficiency
		effT = 0.5;	// time efficiency (proportion pointed at sun)
		string load = path + name + ext;
		m = loadtxt(load,0);
		flux = loadtxt(load,1);
		len = flux.size();
		}

	else if ( detector==1 ) {
		// baseline IAXO parameters
		name="baselineIAXO" + nModel;
		A = 2.3;	// detector area [m2]
		phiBg = 1e-8 * 1e4 * dE;	// background flux [m-2 s-1]
		a = 1.2 * 1e-4;	// XRay detection area [m2]
		t = 3 * 365.25 * 24 * 3600;	// 3 years
		effD = 0.8;	// detectior efficiency
		effO = 0.7;	// optical efficiency
		effT = 0.5;	// time efficiency (proportion pointed at sun)
		string load = path + name + ext;
		m = loadtxt(load,0);
		flux = loadtxt(load,1);
		len = flux.size();
		}

	else if ( detector==2 ) {
		// upgraded IAXO parameters
		name="upgradedIAXO" + nModel;
		A = 3.9;	// detector area [m2]
		phiBg = 1e-9 * 1e4 * dE;	// background flux [m-2 s-1]
		a = 1.2 * 1e-4;	// XRay detection area [m2]
		t = 5 * 365.25 * 24 * 3600;	// 5 years
		effD = 0.8;	// detectior efficiency
		effO = 0.7;	// optical efficiency
		effT = 0.5;	// time efficiency (proportion pointed at sun)
		string load = path + name + ext;
		m = loadtxt(load,0);
		flux = loadtxt(load,1);
		len = flux.size();
		}

	// calculate background count
	double Nbg = phiBg * a * t * effT;
	cout << "Detector: "<< name << endl;
	cout << "Expected background count: " << Nbg << endl << endl;

	// get poisson random number
	default_random_engine generator;
	poisson_distribution<int> distro( Nbg );

	vector<double> chi;	// initialise chi vector

	// get 95% chi for each m value my minimisation
	for ( int c = 0; c < len; c++ ) {
		
		// signal flux for chi = 1 for small dt
		double Nsig = ( flux[c] ) * A * effO * effD * effT * t;
		double total = 0;	// keep total of all runs
		
		// get sample of random N from poisson
		for ( int i = 0; i < samplesize; i++ ) {
			double n = distro(generator);	// get random n from poisson
			//if ( n == 0 ) { continue; }
			//else { total += min( Nbg, Nsig, n ); }
			total += integral( Nbg, Nsig, n );
		}
			
		chi.push_back( total / samplesize );
		cout << c+1 << " out of " << len << ":	" << total / samplesize << endl;
	}
	
	//cout << "chi length: " << chi.size() << "	m length: " << m.size() << endl;
	// write out
	string savename = "data/limits/chamstats-" + name + "-tPlasmon.dat";
	write2D( savename, m, chi );
}



int main(){
	
	for( int n = 1; n <= 10; n++ ) {

	// thread all 3 at same time
	thread t1(chis, 0, n); usleep(100);	// baby
	thread t2(chis, 1, n); usleep(100);	// baseline
	thread t3(chis, 2, n); usleep(100);	// upgraded
	
	t1.join();
	t2.join();
	t3.join();
	}
	cout << "\n¡¡complete!!" << endl;
	return 0;
}
