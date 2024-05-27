// Tom O'Shea 2024

// Primakoff production of scalars in the sun
// LL version, with full screening

#include "utils.h"

using namespace std;

// solar model
vector<double> ne = read("data/ne.dat");		// electron number density [eV3]
vector<double> wp = read("data/wp.dat");		// plasma frequency [eV]
vector<double> Ts = read("data/T.dat");			// solar temperature [eV]
vector<double> r = read("data/r.dat");			// radial distance [eV-1]

int c = 0;
double m = 0;
double thres = 1e-6;
bool print = false;


double RePi( double E, double q, double ni, double mi, double T ) {
	double arg1 = sqrt(mi/2/T)*(E/q + q/2/mi);
	double arg2 = sqrt(mi/2/T)*(E/q - q/2/mi);
	double I = -2*ni/q * sqrt(mi/2/T) * ( gsl_sf_dawson(arg1) - gsl_sf_dawson(arg2) );
	return Zis*Zis*I;
}

double ImPi( double E, double q, double ni, double mi, double T ) {
	double I =  -ni*sqrt(2*pi/mi/T)*mi/q
		* exp( -(pow(mi*E/q,2) + pow(q/2,2))/(2*mi*T) ) * sinh(E/2/T);
	return Zis*Zis*I;
}

double F( double E, double q, int cs ) {
	double Ve = 4*pi*alpha/q/q;
	double T = Ts[cs];
	double I = 0
	for( int i = 0; i<=3; i++ ) {
		ions(i);
		double ni = nis[cs];
		double iPi = ImPi(E,q,ni,mis,T);
		double rPi = RePi(E,q,ni,mis,T);
		if(abs(E)<1e-10) { I += (Ve*ni/q*sqrt(2*pi*mi/T)*exp(-q*q/8/mi/T) /pow((1-Ve*rPi),2)); }
		I += ( -2*Ve/(1-exp(-E/T)) * iPi / ( pow((1-Ve*rPi),2) + pow((Ve*iPi),2) ) );
	}
	return I;
}

// simplified integrand with averaged x=>0
double integrand3( double E, double q, double w, double k, int cs ) {
	if(q==0) { return 0; }
	double qi = sqrt(q*q + k*k);	// averaged to x=>0
	if(qi==0) { return 0; }
	double I = pow(2*pi,-3)/4 * pow(k*k - q*q - qi*qi, 2)/q/qi * F(E,q,cs) * F(w-E,qi,cs)/Mpl/Mpl;
	if(print) {cout<<I<<endl;}
	return I;
}

double q_integral3( double w, double E, int cs ) {
	double k = sqrt(w*w-m*m);
	double total = 0;
	double dq = E/10;
	int i = 0;
	double q = 0;
	double iminus = integrand3(E,q,w,k,cs);
	double iplus = 0;
	while(true) {
		iplus = integrand3(E,q+dq,w,k,cs);
		double I = 0.5*dq*(iplus + iminus);
		total+=I;
		//cout<<"total = "<<total<<"	I = "<<I<<endl;
		if(I/total <= thres) {break;}
		if( i>=1000 && total==0 ) { break; }
		iminus = iplus;
		//if( E==w ) {cout<<"total = "<<total<<"	I = "<<I<<endl;}
		q+=dq;
		i++;
	}
	//cout<<total<<endl;
	return total;
}

double E_integral3( double w, int cs ) {
	double total = 0;
	double dE = w/100;
	double E = 0;
	int i = 0;
	double iminus = q_integral3(w,E,cs);
	double iplus = 0;
	while(true) {
		iplus = q_integral3(w,E+dE,cs);
		double I = 0.5*dE*(iplus + iminus);
		total+=I;
		//cout<<"w = "<<w<<"	E = "<<E<<"	total = "<<total<<"	I = "<<I<<endl;
		if(I/total <= thres) {break;}
		if( i>=1000 && total==0 ) { break; }
		iminus = iplus;
		E+=dE;
		i++;
	}
	//cout<<total<<endl;
	return total;
}

double r_integral3( double w ) {
	double total = 0;
	int dc = 10;
	double iminus = E_integral3(w,0);
	double iplus = 0;
	for( int cs=0; cs<r.size()-1; cs+=dc ) {
		iplus = E_integral3(w,cs+dc);
		double I = 0.5*(r[cs+1]-r[cs])*(iplus + iminus);
		total+=I;
		//cout<<"r = "<<r[cs]/rSolar<<" r_solar	total = "<<total<<endl;
		if(I/total <= thres) {break;}
		iminus = iplus;
	}
	return total;
}

void fixedR( int cs ) {
	vector<double> count, energy;
	string name;
	double dw = 1e0;
	for( double w = dw; w < 2e4; w+=dw ){
		energy.push_back(w);					// eV
		count.push_back( E_integral3(w,cs) );	// Bg-2
		if((int)(w) % (int)(1e2) == 0) { cout<<"w = "<<w/1e3<<"keV of 20keV"<<endl; }
		//cout<<"w = "<<w<<"eV of 20keV"<<endl;
	}
	name = "data/LLspectrum_half.dat";
	write2D( name , energy, count );	
}

void fullR() {
	vector<double> count, energy;
	string name;
	double dw = 1e0;
	for( double w = dw; w < 2e4; w+=dw ){
		energy.push_back(w);			
		count.push_back( r_integral3(w) );
		//if((int)(w) % (int)(1e2) == 0) { cout<<"w = "<<w/1e3<<"keV of 20keV"<<endl; }
		cout<<"w = "<<w/1e3<<"keV of 20keV"<<endl;
	}
	name = "data/LLspectrum--1e-6.dat";
	write2D( name , energy, count );	
}


int main() {
	//double I = E_integral3(20,0);
	//double I = q_integral3(20,20,0);
	//double I = r_integral3(1000);
	//cout<<I<<endl;
	//fixedR(998);	// r = 0.5 rSolar
	fullR();
	return 0;
}
