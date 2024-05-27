// Tom O'Shea 2024

// Primakoff production of scalars in the sun
// LL version, with full screening

#include "utils.h"

using namespace std;

// constants
double mH = 938.272e6;
double mHe4 = 3.727e9;
double mHe3 = 3.0160293*amu;


// solar model
vector<double> ne = read("data/ne.dat");		// electron number density [eV3]
vector<double> nH = read("data/nH.dat");		// electron number density [eV3]

vector<double> nis;
double mis;
double Zis;

void ions( int sel ) {
	if(sel==0) { ni = ne; mi = me; Zi = 1; }
	else if (sel==1) { nis = nH; mi = mH; Zi = 1; }
	else if (sel==2) { nis = nHe3; mi = mHe3; Zi = 2; }
	else if (sel==3) { nis = nHe4; mi = mHe4; Zi = 2; }
}
