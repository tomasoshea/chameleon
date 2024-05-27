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

double RePi( double E, double q, double ni, double mi, double T ) {
	double arg1 = sqrt(mi/2/T)*(E/q + q/2/mi);
	double arg2 = sqrt(mi/2/T)*(E/q - q/2/mi);
	double daw1, daw2;
	if(arg1>=0) { daw1 = Dawson(arg1); }
	else if(arg1<0) { daw1 = -Dawson(-arg1); }
	if(arg2>=0) { daw2 = Dawson(arg2); }
	else if(arg2<0) { daw2 = -Dawson(-arg2); }
	//cout<<daw2<<endl;
	//if(isnan(daw1) or isnan(daw2)) {return 0;}
	//cout<<daw1<<endl;
	return -2*ni/q * sqrt(mi/2/T) * ( daw1 - daw2 );
}

double ImPi( double E, double q, double ni, double mi, double T ) {
	double I =  -ni*sqrt(2*pi/mi/T)*mi/q
		* exp( -(pow(mi*E/q,2) + pow(q/2,2))/(2*mi*T) ) * sinh(E/2/T);
	//if(isnan(I)) {return 0;}
	return I;
}

double F( double E, double q ) {
	double Ve = 4*pi*alpha/q/q;
	double T = Ts[c];
	double ni = ne[c];
	double iPi = ImPi(E,q,ni,me,T);
	double rPi = RePi(E,q,ni,me,T);
	//cout<<rPi<<endl;
	return -2*Ve/(1-exp(-E/T)) * iPi / ( pow((1-Ve*rPi),2) + pow((Ve*iPi),2) );
}

// simplified integrand with averaged x=>0
double integrand3( double E, double q, double w, double k ) {
	if(E==0) {return 0;}
	if(q==0) {return 0;}
	double qi = sqrt(q*q + k*k);	// averaged to x=>0
	if(qi==0) {return 0;}
	double I = 1/4/pow(2*pi,3) * pow(k*k - q*q - qi*qi, 2)/q/qi * F(E,q) * F(w-E,qi);
	cout<<"E=F(E,q)<<endl;
	return I;
}

double q_integral3( double w, double E ) {
	double k = sqrt(w*w-m*m);
	double total = 0;
	double thres = 1e-10;
	double dq = 1;
	while(true) {
		double q = 0;
		double I = 0.5*dq*(integrand3(E,q+dq,w,k) + integrand3(E,q,w,k));
		total+=I;
		cout<<"total = "<<total<<"	I = "<<I<<endl;
		if(I/total <= thres) {break;}
	}
	return total;
}

double E_integral3( double w ) {
	double total = 0;
	double thres = 1e-99;
	double dE = 1;
	while(true) {
		double E = 0;
		double I = 0.5*dE*(q_integral3(w,E+dE) + q_integral3(w,E));
		total+=I;
		if(I/total <= thres) {break;}
	}
	return total;
}

int main() {
	double I = q_integral3(1000,500);
	cout<<I<<endl;
	return 0;
}
