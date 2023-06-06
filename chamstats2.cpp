// Tom O'Shea 2023

// script to find the upper limit on phi for a given number of observed events
// using a poisson distribution, for use in IAXO chameleon analysis

#include "cham2.cpp"

using namespace std;

double CL = 0.95;	// confidence level
double dE = 10;	// E range [keV]
int samplesize = 1e1;		// size of random sample
double nModel = 1.; // model parameter


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

double f( double x, double b, double n, double Bm,
		double L, double B, double A, double effO, double effD,
		double effT, double t ) {

	// get flux from model
	//double s = fullFlux(x,nModel,Bm,L,B) * A * effO * effD * effT * t;
	double s = pow(x,2) * wIntegral(nModel,Bm,L,B) * A * effO * effD * effT * t;
	//cout << s << endl;
	double item;

	if ( n == 0 ) { item = ( b + s ); }
	else if ( n == 1 ) { item = ( ( b + s ) - log( b + s ) - 1 ); }
	else { item = ( ( b + s ) - n * ( log( b + s ) + 1 - log(n) ) ); }

	return exp(-item);
}

// integral for getting CL
double integral( double b, double n, double Bm, double L,
		double B, double A, double effO, double effD,
		double effT, double t ) {

	double x = 1e7;
	double x2 = x;
	//double dx = x;
	//cout << x2 << endl;
	double mx = 2;
	//double total = 0.5 * x * ( exp( - f( x, b, s, n ) ) + exp( - f( 0, b, s, n ) ) );
	double total= 0.5 * x * ( f(x,b,n,Bm,L,B,A,effO,effD,effT,t) + f(0,b,n,Bm,L,B,A,effO,effD,effT,t) );
	double norm = total;

	// normalise with intg to inf
	while(true) {
		double dx = x2 * (mx - 1);
		double L1 = f(x2,b,n,Bm,L,B,A,effO,effD,effT,t);
		double L2 = f(x2+dx,b,n,Bm,L,B,A,effO,effD,effT,t);
		if ( isnan(L1 + L2) ) { continue; }
		else { norm += 0.5 * dx * ( L1 + L2 ); }
		if ( L2 + L1 == 0 ) { break; }
		//x2 += dx;
		x2 *= mx;
		//cout << L1+L2 << endl;
	}

	// integrate up to CL
	while ( ( total / norm ) < CL ) {
		double dx = x * (mx - 1);
		double L1 = f(x,b,n,Bm,L,B,A,effO,effD,effT,t);
		double L2 = f(x+dx,b,n,Bm,L,B,A,effO,effD,effT,t);
		//double L2 = L( x*mx, b, s, n );
		
		if ( isnan(L1 + L2) ) { continue; }
		else { total += 0.5 * dx * ( L1 + L2 ); }
		//cout << total * dx / norm << endl;
		x *= mx;
		//x2 += dx;
	}
	cout << x << endl;
	return x;
}


void chis( int detector ) {

	// initialise parameters
	double A, phiBg, a, t, effD, effO, effT, len, L, B;
	string name;
	string ext = "-cham-flux.dat";
	string path = "data/limits/";

	// choose detector

	if ( detector==0 ) {
		// babyIAXO parameters
		name="babyIAXO";// + nModel;
		A = 0.77;	// detector area [m2]
		phiBg = 1e-7 * 1e4 * dE;	// background flux [m-2 s-1]
		a = 0.6 * 1e-4;	// XRay detection area [m2]
		t = 1.5 * 365.25 * 24 * 3600;	// 1.5 years
		effD = 0.7;	// detectior efficiency
		effO = 0.35;	// optical efficiency
		effT = 0.5;	// time efficiency (proportion pointed at sun)
		L = 10;		// detector length [m]
		B = 2;		// B-field [T]
		}

	else if ( detector==1 ) {
		// baseline IAXO parameters
		name="baselineIAXO";// + nModel;
		A = 2.3;	// detector area [m2]
		phiBg = 1e-8 * 1e4 * dE;	// background flux [m-2 s-1]
		a = 1.2 * 1e-4;	// XRay detection area [m2]
		t = 3 * 365.25 * 24 * 3600;	// 3 years
		effD = 0.8;	// detectior efficiency
		effO = 0.7;	// optical efficiency
		effT = 0.5;	// time efficiency (proportion pointed at sun)
		L = 20;		// detector length [m]
		B = 2.5;		// B-field [T]
		}

	else if ( detector==2 ) {
		// upgraded IAXO parameters
		name="upgradedIAXO";// + nModel;
		A = 3.9;	// detector area [m2]
		phiBg = 1e-9 * 1e4 * dE;	// background flux [m-2 s-1]
		a = 1.2 * 1e-4;	// XRay detection area [m2]
		t = 5 * 365.25 * 24 * 3600;	// 5 years
		effD = 0.8;	// detectior efficiency
		effO = 0.7;	// optical efficiency
		effT = 0.5;	// time efficiency (proportion pointed at sun)
		L = 22;		// detector length [m]
		B = 3.5;		// B-field [T]
		}


	// calculate background count
	double Nbg = phiBg * a * t * effT;
	cout << "Detector: "<< name << endl;
	cout << "Expected background count: " << Nbg << endl << endl;

	// get poisson random number
	default_random_engine generator;
	poisson_distribution<int> distro( Nbg );

	vector<double> chi, m;	// initialise chi vector

	// get 95% chi for each Bm value my minimisation
	for ( double Bm = 1e10; Bm < 1e25; Bm*=10 ) {
		
		double total = 0;	// keep total of all runs
		// get sample of random N from poisson
		for ( int i = 0; i < samplesize; i++ ) {
			double n = distro(generator);	// get random n from poisson
			total += integral(Nbg,n,Bm,L,B,A,effO,effD,effT,t);
		}
			
		chi.push_back( total / samplesize );
		m.push_back(Bm);
		cout << name << ":	Bm = " << Bm << "	Bg = " << total/samplesize << endl;
	}
	
	//cout << "chi length: " << chi.size() << "	m length: " << m.size() << endl;
	// write out
	string savename = "data/limits/symstats2-" + name + "-tPlasmon.dat";
	write2D( savename, m, chi );
}


int main(){

	// thread all 3 at same time
	thread t1(chis, 0); usleep(100);	// baby
	//thread t2(chis, 1); usleep(100);	// baseline
	//thread t3(chis, 2); usleep(100);	// upgraded
	
	t1.join();
	//t2.join();
	//t3.join();

	cout << "\n¡¡complete!!" << endl;
	return 0;
}
