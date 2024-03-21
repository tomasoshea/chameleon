// Tom O'Shea 2024

// Primakoff production of scalars in the sun
// V3 - the one using full ee scatter with full disp. rel. but no Raffelt screening=

#include "utils.h"

using namespace std;

// constants
double pi = 3.14159265359;
double alpha = 1/137.035999084;
double me = 510998.950;		// e- mass [eV]
double Mpl = 2e27;   		// planck mass [eV]
double rSolar = 6.957e8/m2eV;					// solar radius [eV-1]
double dSolar = 149.5978707e9/(1.973269804e-7);		// mean earth-sun distance [eV-1]
int nancount = 0;
int ucount = 0;
double rhoLead = 11.34*1e3*kg2eV*pow(m2eV,3);	// density of lead (11.34 g cm-3) [eV4]

// solar params
vector<double> ne = read("data/ne.dat");		// electron number density [eV3]
vector<double> wp = read("data/wp.dat");		// electron number density [eV3]
vector<double> T = read("data/T.dat");			// solar temperature [eV]
vector<double> r = read("data/r.dat");			// radial distance [eV-1]
vector<double> rho = read("data/rho.dat");		// solar density [eV4]
vector<double> nH = read("data/nH.dat");		// H number density [eV3]
vector<double> nHe3 = read("data/nHe3.dat");	// He3 number density [eV3]
vector<double> nHe4 = read("data/nHe4.dat");	// He4 number density [eV3]
// get gaunt factors
vector<vector<double>> z1 = readGaunt("data/Z1.dat");	// gaunt factors for Z=1
vector<vector<double>> z2 = readGaunt("data/Z2.dat");	// gaunt factors for Z=2

// variables
double E = 2.4e-3;    		// cham potential energy scale [eV] (2.4e-3 for cosmological const)
double n = 1;				// chameleon potential 1/phi^n
double Bm = 1e2;			// cham matter coupling
double B_iaxo = 0;			// IAXO B-field [eV2]
double L_iaxo = 0;			// IAXO conversion length [eV-1]



// chameleon mass as a function of solar radius and model parameters
// cham model params n (phi-n potential), Bm (matter coupling)
// assume rho dominated by matter density
double mCham2( double rho ) {
	//if(n==0) { cout<<"ERROR! n = 0"<<endl; return 0; }
	double E4n = pow(E,4+n);
	if(n<0) { return n*(n+1)*E4n*pow( pow( Bm*rho/(n*Mpl*E4n), (n+2) ) , 1/(n+1) ); }
	else { return n*(n+1)*E4n* pow( Bm*rho/(n*Mpl*E4n), (n+2)/(n+1) ); }
}

// T plasmon absorbtion length
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


// integral I(u,v) solved analytically in cross section calc
// dimensionless, with dimensionless arguments
double curlyA( double a, double b ) {
	return a*log( (pow(a+1,2) + b*b)/(pow(a-1,2) + b*b) )
			+ (b*b - a*a + 1)/b *(atan((1-a)/b) + atan(1-a)/b)
			- 2;
}
double curlyAapprox( double u, double v ) {
	return 4*u/(v*v + 2)
			- (u*u - v*v - 1)/v *( atan((1-u)/v) - atan((1+u)/v) )
			- 2;
}

double curlyB( double a ) {
	return 2*log( (a+1)/(a-1) ) - 4;
}

double curlyD( double y, double v) {
	double u = 1/(2*y) + y/2;
	return pow(u-y,2)/v *log((u+1)/(u-1))
			- pow(u+v-y,2)/v*log((u+v+1)/(u+v-1))
			+ 2;
}

// integral I(u,v) solved analytically in cross section calc
// dimensionless, with dimensionless arguments
double curlyI( double u, double v ) {
	return (u*u - 1)/v*log((u-1)/(u+1)) - (pow(u+v,2) - 1)/v*log((u+v-1)/(u+v+1)) - 2;
}
double curlyIapprox( double u, double v ) {		// for u->1
	return u*u/v - (v+2)*log(v/(v+2)) - 2;
}

// differential scalar production rate on earth d2N/dr/dw times Lambda2
// units eV Bm-2
double integrand( int c, double w ) {
	if( T[c]==0 ) { return 0; }				// solves weird behaviour when ne = T = 0
	double mg2 = 4*pi*alpha*ne[c]/me;		// assume mg2 = wp2
	double ms2 = mCham2(c);					// chameleon mass2 [eV2]
	//cout<<ms2<<endl;
	//double ms2 = Bm*Bm;						// fixed scalar mass2 [eV2]
	if( w*w <= mg2 ) { return 0; }
	if( w*w <= ms2 ) { return 0; }
	double K2 = 8*pi*alpha*ne[c]/T[c];		// Debye screening scale ^2 [eV2]
	double kgamma = sqrt(w*w - mg2);		// photon momentum [eV]
	double kphi = sqrt(w*w - ms2);			// scalar momentum [eV]
	double uArg = kgamma/(2*kphi) + kphi/(2*kgamma);	// u for curlyI
	double vArg = K2/(2*kphi*kgamma);		// v for curlyI
	double Iuv = curlyI(uArg,vArg);
	if(uArg < 1.01) { Iuv = curlyIapprox(uArg,vArg); }
//	double aArg = (2*w*w - ms2)/(2*kgamma*kphi);	// u for curlyA
//	double bArg = w*Ggamma/(2*kphi*kgamma);		// v for curlyA
	//cout<<bArg<<endl;
	//double Aab = curlyB(aArg);
	//if( aArg/bArg < 1e3 ) { Aab = curlyA(aArg,bArg); }
	//if(aArg < 1.01) { Aab = curlyIapprox(aArg,bArg); }

	return alpha/(8*Mpl*Mpl*pi) * pow(r[c], 2) * ne[c]/(exp(w/T[c]) - 1) 
			* w*w * kphi/kgamma * Iuv;		// [eV Bg^-2]
}


// differential scalar production rate on earth d2N/dw/dk times Bg2
// pure longitudinal component only
// units Lambda2
double L_integrand( int c, double kgamma ) {
	if( T[c]==0 ) { return 0; }				// solves weird behaviour when ne = T = 0
	double w = wp[c];						// omega is plasma freq
	double ms2 = mCham2(c);					// chameleon mass2 [eV2]
	//double ms2 = Bm*Bm;						// fixed scalar mass2 [eV2]
	if( w*w <= ms2 ) { return 0; }
	double K2 = 8*pi*alpha*ne[c]/T[c];		// Debye screening scale ^2 [eV2]
	//double kgamma = sqrt(w*w - mg2);		// photon momentum [eV]
	double kphi = sqrt(w*w - ms2);			// scalar momentum [eV]
	double yArg = kgamma/kphi;				// y for curlyA
	double vArg = K2/(2*kphi*kgamma);		// v for curlyA
	double Dyuv = curlyD(yArg,vArg);
	//cout << Dyuv << endl;

	//return alpha/(8*Mpl*Mpl*pi) * pow(r[c], 2) * ne[c]/(exp(w/T[c]) - 1) 
	//		* w*w*kphi/kgamma * Dyuv;		// [eV2]
	
	return alpha/(2*Mpl*Mpl*pi) * pow(r[c], 2) * ne[c] * T[c]
			* kphi * Dyuv;		// [eV2]
}

// integral over solar volume, for a given scalar mass and energy
// returns dPhi/dw Bg-2 [eV2]
double solarIntg( double w ) {
	double total = 0;
	for( int c = 0; c < r.size() - 1; c++ ) {
		total += 0.5 * (r[c+1] - r[c]) * (integrand(c+1,w) + integrand(c,w));
	}
	return total;
}

// integral over k_gamma for l-plasmon
double kIntg( int c ) {
	double total = 0;
	double kD = sqrt(8*pi*alpha*ne[c]/T[c]);
	double dk = kD/1000;
	for( double k = dk; k < kD; k+= dk ) {
		//cout << k << endl;
		total += 0.5 * dk * (L_integrand(c, k+dk) + L_integrand(c, k));
	}
	return total;
}


// calculate differential particle flux spectrum by intg over solar volume
// dN/dw, units Bg-2
// (where dN is really dN/dt, sorry)
void spectrum() {
	vector<double> count, energy;
	Bm = 1e2;
	n = 1;
	E = 2.4e-3;
	double dw = 1e0;
	for( double w = dw; w < 2e4; w+=dw ){
		energy.push_back(w);					// eV
		count.push_back( solarIntg(w) );		// Bg-2
		if((int)(w) % (int)(1e3) == 0) { cout<<"w = "<<w/1e3<<"keV of 20keV"<<endl; }
		//if((int)(w) % (int)(100*dw) == 0) { cout<<"w = "<<w<<"eV of 1000eV"<<endl; }
	}
	// write to file
	string name = "data/primakoffV3_spectrum_1e2.dat";
	write2D( name , energy, count );
}

void L_spectrum() {
	vector<double> count, energy;
	Bm = 1e2;
	n = 1;
	E = 2.4e-3;
	double w1, w2 = 0;
	double r1, r2 = rSolar;
	for( int j = wp.size()-1; j >=0; j-- ){
		w1 = wp[j];
		if(w2 > w1) { continue; }
		else{
		r1 = r[j];
		energy.push_back(w1);
		count.push_back( kIntg(j) * abs((r2-r1)/(w2-w1)) );
		r2 = r[j];
		w2 = wp[j];
		}
	}
	// write to file
	string name = "data/primakoffV3_L_spectrum_1e2.dat";
	write2D( name , energy, count );
}


// get combined LT spectrum
void total_spectrum() {
	vector<double> count, energy;
	Bm = 1e2;
	n = 1;
	E = 2.4e-3;
	double w1, w2 = 0;
	double r1, r2 = rSolar;
	double dw = 1e0;
	for( int j = wp.size()-1; j >=0; j-- ){
		w1 = wp[j];
		if(w2 > w1) { continue; }
		else{
		r1 = r[j];
		energy.push_back(w1);
		count.push_back( ( kIntg(j) * abs((r2-r1)/(w2-w1))	// L
						+ solarIntg(w1) ) );					// T
		r2 = r[j];
		w2 = wp[j];
		}
	}
	for( double w = w1+dw; w < 2e4; w+=dw ){
		energy.push_back(w);
		count.push_back( solarIntg(w) );
		if((int)(w) % (int)(1e3) == 0) { cout<<"w = "<<w/1e3<<"keV of 20keV"<<endl; }
		//if((int)(w) % (int)(100*dw) == 0) { cout<<"w = "<<w<<"eV of 1000eV"<<endl; }
	}
	// write to file
	string name = "data/primakoffV3_total-spectrum_cham_1e2_nMinus4.dat";
	write2D( name , energy, count );
}


// chameleon mass (so minimum energy) in lead shielding
// function of chameleon parameters
// units eV
void lead() {
	vector<double> beta, mass;
	double rho = rhoLead;
	n = 1;
	E = 2.4e-3;
	double E4n = pow(E,4+n);
	double mB = 1.1;
	for( Bm = 1e-10; Bm <= 1e1; Bm*=mB ) {
		beta.push_back(Bm);
		mass.push_back( sqrt(mCham2(rhoLead)) );
	}
	// write to file
	string name = "data/lead_n1.dat";
	write2D( name ,beta , mass );
}
		
		

// back conversion prob
// naive version, with low m but a = 0
// (becomes indpt of both w and m)
// dimensionless
double backConversion() {
	return pow(B_iaxo*L_iaxo/2/Mpl, 2);
}


// get IAXO back-converted flux
// input detector params B, L
// units eV3/keV
void IAXO( double Bin, double Lin, string model ) {
	vector<double> beta, flux;
	n = 1;
	E = 2.4e-3;
	B_iaxo = Bin;
	L_iaxo = Lin;
	double Emin = 1e3;
	double Emax = 2e4;
	double P = backConversion();
	double dw = 1e1;
	double w1, w2, r1, r2 = 0;
	for( double Bm = 1e-1; Bm <= 1e18; Bm*=10 ) {
		double total = 0;
		if( Emin < 300 ) {
		for( int j = wp.size()-1; j >= 0; j-- ){
			w1 = wp[j];
			if(w2 >= w1) { continue; }
			else{
			r1 = r[j];
			total += 0.5*(w1-w2)*( (w1+w2)*( kIntg(j) * abs((r2-r1)/(w2-w1)) )
								 + w1*solarIntg(w1) + w2*solarIntg(w2) );
			if(isnan(total)) {cout<<"Bm = "<<Bm<<"	w = "<<w1<<"	j = "<<j<<endl;}
			r2 = r[j];
			w2 = wp[j];
			}
		}
		for( double w = w1+dw; w < 2e4; w+=dw ){
			total += 0.5*dw*( (w+dw)*solarIntg(w+dw) + w*solarIntg(w) );
		}
		}
		else{ for( double w = Emin; w < Emax; w+=dw ){ 
			total += 0.5*dw*(solarIntg(w+dw)+solarIntg(w)); } 
		}
		beta.push_back(Bm);
		flux.push_back(total/4/pi/dSolar/dSolar / ((Emax-Emin)/1e3));
		cout << "Bm = " << Bm << " of 1e18" << endl;
	}
	// write to file
	string name = "data/"+model+"_totalflux.dat";
	write2D( name, beta, flux );
}

int main() { 
	double B_baby = 2*T2eV;	// babyIAXO B-field [eV2]
	double L_baby = 10/m2eV;	// babyIAXO length [eV-1]
	
	double B_base = 2.5*T2eV;	// babyIAXO B-field [eV2]
	double L_base = 20/m2eV;	// babyIAXO length [eV-1]
	
	double B_plus = 3.5*T2eV;	// babyIAXO B-field [eV2]
	double L_plus = 22/m2eV;	// babyIAXO length [eV-1]
	
	double B_cast = 9*T2eV;	// babyIAXO B-field [eV2]
	double L_cast = 9.26/m2eV;	// babyIAXO length [eV-1]
	
	IAXO(B_plus, L_plus, "plusIAXO");
	return 0;
}
