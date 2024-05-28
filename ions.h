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
vector<double> nHe3 = read("data/nHe3.dat");		// electron number density [eV3]
vector<double> nHe4 = read("data/nHe4.dat");		// electron number density [eV3]


vector<double> nis;
double mis;
double Zis;

void ions( int sel ) {
	if(sel==0) { nis = ne; mis = me; Zis = 1; }
	else if (sel==1) { nis = nH; mis = mH; Zis = 1; }
	else if (sel==2) { nis = nHe3; mis = mHe3; Zis = 2; }
	else if (sel==3) { nis = nHe4; mis = mHe4; Zis = 2; }
}
