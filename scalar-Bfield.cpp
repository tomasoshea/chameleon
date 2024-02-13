// Tom O'Shea 2023

// Primakoff production of scalars in the sun

#include "utils.h"

using namespace std;

// constants
double pi = 3.14159265359;
double alpha = 1/137.035999084;
double me = 510998.950;							// e- mass [eV]
double Mpl = 2e27;   							// planck mass [eV]
double E = 2e-3;    							// cham potential energy scale [eV] (2.4e-3 for cosmological const)
double n = 1;									// chameleon potential 1/phi^n
double zeta3 = 1.202056903159594;				// Riemann zeta(3)

double dSolar = 149.5978707e9/m2eV;				// mean earth-sun distance [eV-1]
double rSolar = 6.957e8/m2eV;					// solar radius [eV-1]
double B0 = 3e3*T2eV;							// radiative zone max B [eV2]
double B1 = 50*T2eV;							// tachocline max B [eV2]
double B2 = 3*T2eV;								// outer region max B [T]
double r0 = 0.712*rSolar;						// [eV-1]
double r1 = 0.732*rSolar;						// [eV-1]
double d1 = 0.02*rSolar;						// [eV-1]
double r2 = 0.96*rSolar;						// [eV-1]
double d2 = 0.035*rSolar;						// [eV-1]
double ngamma0 = 1e25*m2eV*m2eV*s2eV;			// photon flux at r0 [eV3]

// solar params
vector<double> ne = read("data/ne.dat");		// electron number density [eV3]
vector<double> T = read("data/T.dat");			// solar temperature [eV]
vector<double> r = read("data/r.dat");			// radial distance [eV-1]
vector<double> rho = read("data/rho.dat");		// solar density [eV4]
vector<double> nH = read("data/nH.dat");		// H number density [eV3]
vector<double> nHe3 = read("data/nHe3.dat");	// He3 number density [eV3]
vector<double> nHe4 = read("data/nHe4.dat");	// He4 number density [eV3]
// get gaunt factors
vector<vector<double>> z1 = readGaunt("data/Z1.dat");	// gaunt factors for Z=1
vector<vector<double>> z2 = readGaunt("data/Z2.dat");	// gaunt factors for Z=2

int core = 0;

// chameleon mass as a function of solar radius and model parameters
// cham model params n (phi-n potential), Bm (matter coupling)
// assume rho dominated by matter density
double mCham2( int c, double Bm ) {
	double E4n = pow(E,4+n);
	return n*(n+1)*E4n*pow( Bm*rho[c]/(n*Mpl*E4n), (n+2)/(n+1) );
}


// T plasmon absorbtion length [eV]
double GammaPhoton( double w, int c, double g1, double g2 ) {

	double p1 = 64 * pow(pi,2) * pow(alpha,3);
	double p2 = 3 * pow(me,2) * pow(w,3);
	double p3 = me * pow(ne[c],2) / (2*pi*T[c]);
	double p4 = 1 - exp(- w / T[c]);
	double p5 = 8 * pi * pow(alpha,2) * ne[c] / (3 * pow(me,2) );

	// sum of ion densities
	double ions = (nH[c] * g1) + g2 * ( (4 * nHe4[c]) + (4 * nHe3[c]) );

	return p1 * pow(p2, -1) * pow(p3, 0.5) * p4 * ions + p5;
}


// return B field
// units eV2
double Bfield( int c ) {
	// B field in solar radiative zone
	if ( r[c] <= r0 ) {
		double l = (10 * (r0/rSolar)) + 1;
		double K = (1+l) * pow( 1 + pow(l,-1) , l ) * B0;
		core=1;
		return K * pow(r[c]/r0 , 2) * pow(1 - pow(r[c]/r0 , 2), l);
		}
	// B-field in tachocline
	else if ( r[c] > (r1 - d1) and r[c] < (r1 + d1)) {
		//cout << "	T" << endl;
		core=0;
		return B1 * ( 1 - pow( (r[c] - r1)/d1 , 2 ) );
		}
	// B-field in outer region
	else if ( r[c] > (r2 - d2) and r[c] < (r2 + d2) ) {
		core=0;
		//cout << "		O" << endl;
		return B2 * ( 1 - pow( (r[c] - r2)/d2 , 2 ) );
		}
	else { return 0; }
}

// solved integral I(a)
double IB(double a) {
	return sqrt(pi/2)*(sqrt(a+sqrt(a*a+4)) - sqrt(2*a));
}

// Gamma boosted by Luca factor for photon scatter
double boostedGamma( int c, double G, double kgamma, double kphi ) {
	double q = sqrt(abs(kgamma*kgamma - kphi*kphi));
	return G /( ngamma0*pow(r0/r[c],2)*pi/4/zeta3*pow(T[c],-3)
			*sqrt(q/2/m2eV) * IB(2*G/q) );
}


// differential scalar production rate d2N/dr/dw times Lambda2
// units eV Bg-2
double integrand( int c, double Bm, double w, double G ) {
	if( T[c]==0 ) { return 0; }					// solves weird behaviour when ne = T = 0
	double mg2 = 4*pi*alpha*ne[c]/me;			// assume mg2 = wp2 [eV2]
	double ms2 = mCham2(c,Bm);					// chameleon mass2 [eV2]
	//double ms2 = Bm*Bm;							// fixed scalar mass2 [eV2]
	if( w*w <= mg2 ) { return 0; }
	if( w*w <= ms2 ) { return 0; }
	double kgamma = sqrt(w*w - mg2);			// photon momentum [eV]
	double kphi = sqrt(w*w - ms2);				// scalar momentum [eV]
	double B = Bfield(c);						// solar B field [eV2]
	//if( core==1 ) { cout<<B<<endl; }
	//double Gboost = boostedGamma(c,G,kgamma,kphi);

	return 1/(pi*Mpl*Mpl) * pow(r[c], 2) *B*B * w*pow(w*w - ms2, 3/2)/( pow(ms2 - mg2, 2) + (w*w*G*G) )
			* G/(exp(w/T[c]) - 1);	// [eV Bg-2]
}


// integral over solar volume, for a given scalar mass and energy
double solarIntg( double w, double Bm ) {
	double total = 0;
	for( int c = 0; c < r.size() - 1; c++ ) {
		// select g(w, T) value from matrix
		int indexT1;
		int indexT2;
		int indexX1;
		int indexX2;
		for( int i = 1; i < 200; i++ ) {
			if( z1[0][i] < T[c] and z1[0][i+1] > T[c] ) { indexT1 = i; }
			if( z2[0][i] < T[c] and z2[0][i+1] > T[c] ) { indexT2 = i; }
		}
		for( int i = 1; i < 500; i++ ) {
			if( (z1[i][0] * T[c]) < w and (z1[i+1][0] * T[c]) > w ) { indexX1 = i; }
			if( (z2[i][0] * T[c]) < w and (z2[i+1][0] * T[c]) > w ) { indexX2 = i; }
		}
		double g1 = z1[ indexT1 ][ indexX1 ];
		double g2 = z2[ indexT2 ][ indexX2 ];
		double G = GammaPhoton(w, c, g1, g2);

		total += 0.5 * (r[c+1] - r[c]) * (integrand(c+1, Bm, w, G) + integrand(c, Bm, w, G));
		if(r[c+1] < r[c]) { cout<<r[c+1]<<"	"<<r[c]<<"	"<<c<<endl; }
	}
	return total;
}


// calculate emission rate profile over solar radius (dN/dr)
// integrated over relevant energies up to 20 keV
void profile() {
	vector<double> radius;
	vector<double> rate;
	double Bm = 1e3;
	double mw = 1.1;
	double dw = 10;		// [eV]
	for( int c = 0; c < r.size(); c++ ) {
		double total = 0;
		// select g(w, T) value from matrix
		int indexT1;
		int indexT2;
		for( int i = 1; i < 200; i++ ) {
			if( z1[0][i] < T[c] and z1[0][i+1] > T[c] ) { indexT1 = i; }
			if( z2[0][i] < T[c] and z2[0][i+1] > T[c] ) { indexT2 = i; }
		}
		for( double w = dw; w < 2e4; w+=dw ) {
			int indexX1;
			int indexX2;
			for( int i = 1; i < 500; i++ ) {
				if( (z1[i][0] * T[c]) < w and (z1[i+1][0] * T[c]) > w ) { indexX1 = i; }
				if( (z2[i][0] * T[c]) < w and (z2[i+1][0] * T[c]) > w ) { indexX2 = i; }
			}
			double g1 = z1[ indexT1 ][ indexX1 ];
			double g2 = z2[ indexT2 ][ indexX2 ];
			double G = GammaPhoton(w, c, g1, g2);

			total += 0.5*dw*( integrand(c,Bm,w*mw,G) + integrand(c,Bm,w,G) );
		}
		//if( r[c] <= r0 ) { cout<<total<<endl;}
		radius.push_back(r[c]);
		rate.push_back(total);
	}
	// write to file
	string name = "data/scalarB_profile_cham-1e3.dat";
	write2D( name , radius, rate );
}


// calculate energy loss as a function of m
// units eV Bm-2
void Eloss() {
	vector<double> mass;
	vector<double> Q;
	double dw = 1e2;
	for( double Bm = 1e0; Bm < 1e8; Bm*=10 ) {
		double total = 0;
		for( double w = dw; w < 2e4; w+=dw ){
			total += 0.5*dw*( (w+dw)*solarIntg(w+dw,Bm) + w*solarIntg(w,Bm) );
		}
		mass.push_back(Bm);
		Q.push_back(total);
		//if((int)(log10(Bm)) % 1 == 0) { cout<<"Bm = 1e"<<(int)(log10(Bm))<<" of 1e8"<<endl; }
		cout<<"Bm = "<<Bm<<endl;
	}
	// write to file
	string name = "data/scalarB_Eloss_cham_1e3--boost.dat";
	write2D( name , mass, Q );
}


// calculate differential particle flux spectrum by intg over solar volume
void spectrum() {
	vector<double> count, energy;
	double Bm = 1e3;		// cham matter coupling, or fixed scalar mass
	double dw = 1e0;
	for( double w = dw; w < 2e4; w+=dw ){
		energy.push_back(w);
		count.push_back( solarIntg(w,Bm) /(4*pi*dSolar*dSolar) );
		if( (int)(w) % (int)(1e3) == 0 ) {
		cout << "w = "<<w/1e3<<"keV (of 20keV)" << endl;
		}
	}
	// write to file
	string name = "data/scalarB_spectrum_cham_1e3.dat";
	write2D( name , energy, count );
}


int main() { 
	// convert Gaunt factor Theta to T in eV
	for( int i = 1; i < 201; i++ ) { z1[0][i] = z1[0][i] * me; }
	for( int i = 1; i < 201; i++ ) { z2[0][i] = z2[0][i] * me; }

	//profile();
	spectrum();
	//Eloss();
	return 0;
}
