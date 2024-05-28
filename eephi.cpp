// Tom O'Shea 2024

// Primakoff production of scalars in the sun
// LL version, with full screening

#include "ions.h"	// also contains utils.h
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_integration.h>

using namespace std;

// solar model
vector<double> wp = read("data/wp.dat");		// plasma frequency [eV]
vector<double> Ts = read("data/T.dat");			// solar temperature [eV]
vector<double> r = read("data/r.dat");			// radial distance [eV-1]

int c = 0;
double m = 0;
double ms2 = 0;
double thres = 1e-4;
bool print = false;


double RePi( double E, double q, double ni, double mi, double T ) {
	double arg1 = sqrt(mi/2/T)*(E/q + q/2/mi);
	double arg2 = sqrt(mi/2/T)*(E/q - q/2/mi);
	double daw1, daw2;
	daw1 = gsl_sf_dawson(arg1);
	daw2 = gsl_sf_dawson(arg2);
	double I = -2*ni/q * sqrt(mi/2/T) * ( daw1 - daw2 );
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
	//double I = 0;
	//for( int i = 0; i<=3; i++ ) {
	//	ions(i);
	//	double ni = nis[cs];
	//	double iPi = ImPi(E,q,ni,mis,T);
	//	double rPi = RePi(E,q,ni,mis,T);
	//	if(abs(E)<1e-10) { I += (Ve*ni/q*sqrt(2*pi*mi/T)*exp(-q*q/8/mi/T) /pow((1-Ve*rPi),2)); }
	//	I += ( -2*Ve/(1-exp(-E/T)) * iPi / ( pow((1-Ve*rPi),2) + pow((Ve*iPi),2) ) );
	//}
	//return I;
	ions(0);
	double ni = nis[cs];
	double iPi = ImPi(E,q,ni,mis,T);
	double rPi = RePi(E,q,ni,mis,T);
	//if(abs(E)<1e-10) {return 0;}
	if(abs(E)<1e-10) { return (Ve*ni/q*sqrt(2*pi*mis/T)*exp(-q*q/8/mis/T) /pow((1-Ve*rPi),2)); }
	else { return ( -2*Ve/(1-exp(-E/T)) * iPi / ( pow((1-Ve*rPi),2) + pow((Ve*iPi),2) ) ); }

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

struct args{
	double q;
	double E;
	double w;
	double k;
	int cs;
};

double integrand3_gsl( double x, void * params ) {
	args& recs = *(args*)params;
	double E = recs.E;
	double q = recs.q;
	double w = recs.w;
	double k = recs.k;
	double cs = recs.cs;
	double rc = r[cs];
	double yArg = q/k;
	double uArg = 1/2/yArg + yArg/2;
	if(q==0) { return 0; }
	double qi = sqrt(q*q + k*k - 2*x*q*k);
	if(qi==0) { return 0; }
	double I = pow(k*rc/2/pi/Mpl,2) * q * pow(x-yArg,2)/(uArg-x) * F(E,q,cs) * F(w-E,qi,cs);
	return I;
}

double x_integral3_gsl( double E, void * params ) {
	args& recs = *(args*)params;
	recs.E = E;
	//cout<<recs.E<<endl;
	gsl_integration_workspace * work = gsl_integration_workspace_alloc (1000);
	double result, error;
	gsl_function F;
	F.function = &integrand3_gsl;
	F.params = &recs;
	//gsl_integration_qagiu(&F, 0, 0, 1e-6, 1000, work, &result, &error);
	//gsl_integration_qagi(&F, 0, 1e-2, 10000, work, &result, &error);
	gsl_integration_qag(&F,-1,1,0,1e-3,1000,6,work,&result,&error);
	return result;
}

double E_integral3_gsl( double q, void * params ) {
	args& recs = *(args*)params;
	recs.q = q;
	gsl_integration_workspace * work = gsl_integration_workspace_alloc (1000);
	double result, error;
	gsl_function F;
	F.function = &x_integral3_gsl;
	F.params = &recs;
	//gsl_integration_qagiu(&F, 0, 0, 1e-6, 1000, work, &result, &error);
	gsl_integration_qagi(&F, 0, 1e-2, 1000, work, &result, &error);
	//gsl_integration_qag(&F,-1,1,0,1e-5,1000,6,work,&result,&error);
	return result;
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
	//cout<<r.size()/32<<endl;
	double total = 0;
	int dc = (int) (r.size()/100);
	double iminus = E_integral3(w,0);
	double iplus = 0;
	for( int cs=0; cs<r.size()-2; cs+=dc ) {
		iplus = E_integral3(w,cs+dc);
		double I = 0.5*(r[cs+dc]-r[cs])*(iplus + iminus);
		total+=I;
		//cout<<"r = "<<r[cs]/rSolar<<" r_solar	total = "<<total<<endl;
		iminus = iplus;
	}
	return total;
}


double A( double y, double u ) { return pow(u-y,2)*log((u+1)/(u-1)) - 2*u + 4*y; }


double lowq_integrand( int c, double q ) {
	double Tc = Ts[c];
	double rc = r[c];
	if( Tc==0 ) { return 0; }				// solves weird behaviour when ne = T = 0
	double mg = wp[c];
	//if( 2*mg2 <= ms2 ) { cout<<ms2-(2*mg2)<<endl; return 0; }
	double yArg = q/sqrt(4*mg*mg - ms2);
	double uArg = 1/2/yArg + yArg/2;
	return 1/(4*Mpl*Mpl) * Tc*Tc * rc*rc * (4*mg*mg - ms2) * q * A(yArg,uArg);	// /4/pi/pi
}


//double highq_integrand( int c, double q ) {
//
//}


void fixedR( int cs ) {
	vector<double> count, energy;
	string name;
	double dw = 1e1;
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
	double dw = 1e1;
	for( double w = dw; w < 2e3; w+=dw ){
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
	args qargs;
	double w = 300;
	qargs.k = w;
	qargs.w = w;
	qargs.cs = 0;
	qargs.q = 1000;
	double I = E_integral3_gsl(qargs.q,&qargs);
	//double I = r_integral3(1000);
	cout<<I<<endl;

	//fixedR(998);	// r = 0.5 rSolar
	//fullR();
	return 0;
}
