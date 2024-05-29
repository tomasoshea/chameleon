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


double E_integral3( double q, double w, int cs ) {
	double k = sqrt(w*w - ms2);
	double total = 0;
	double dE = w/100;
	double E = 0;
	int i = 0;
	double iminus = integrand3(E,q,w,k,cs);
	double iplus = 0;
	while(true) {
		iplus = integrand3(E+dE,q,w,k,cs);
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

double q_integral3( double w, int cs, double qmin ) {
	double total = 0;
	double q = qmin;
	double dq = qmin/10;
	if( dq == 0 ) { dq = w/10; }
	int i = 0;
	double iminus = E_integral3(q,w,cs);
	//cout<<iminus<<endl;
	double iplus = 0;
	while(true) {
		//cout<<q<<endl;
		iplus = E_integral3(q+dq,w,cs);
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


double r_integral3( double w, double qmin ) {
	//cout<<r.size()/32<<endl;
	double total = 0;
	int dc = (int) (r.size()/16);
	//int dc = 100;
	double iminus = q_integral3(w,0,qmin);
	//cout<<iminus<<endl;
	double iplus = 0;
	for( int cs=0; cs+dc<r.size()-1; cs+=dc ) {
		//cout<<r[cs]/rSolar<<endl;
		iplus = q_integral3(w,cs+dc,qmin);
		double I = 0.5*(r[cs+dc]-r[cs])*(iplus + iminus);
		total+=I;
		iminus = iplus;
	}
	return total;
}


/*
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
*/


double A( double y, double u ) {
	if(abs(u-1)<1e-5) { return 0; }
	return pow(u-y,2)*log((u+1)/(u-1)) - 2*u + 4*y; }


double lowq_integrand( int c, double q ) {
	double Tc = Ts[c];
	double rc = r[c];
	if( Tc==0 ) { return 0; }				// solves weird behaviour when ne = T = 0
	double mg = wp[c];
	//if( 2*mg2 <= ms2 ) { cout<<ms2-(2*mg2)<<endl; return 0; }
	double Ayu;
	double yArg = q/sqrt(4*mg*mg - ms2);
	//cout<<yArg<<endl;
	double uArg = 1/2/yArg + yArg/2;
	Ayu = A(yArg,uArg);
	if(isnan(Ayu)) { return 0; }
	return 1/(4*Mpl*Mpl) * Tc*Tc * rc*rc * (4*mg*mg - ms2) * q * Ayu;	// /4/pi/pi
}

double lowq_integral( int c, double qmax ) {
	double total = 0;
	double q = 0;
	double dq = qmax/100;
	int i = 0;
	double iminus = lowq_integrand(c,q);
	double iplus = 0;
	while(q<qmax) {
		//cout<<q/qmax<<endl;
		iplus = lowq_integrand(c,q+dq);
		double I = 0.5*dq*(iplus + iminus);
		total+=I;
		iminus = iplus;
		q+=dq;
		i++;
	}
	//cout<<total<<endl;
	return total;
}

void spectrum_2parts(){
	vector<double> count1, count2, count, energy, energy1;
	double w1, w2 = 0;
	double r1, r2 = rSolar;
	double dj = wp.size()/32;
	for( int j = wp.size()-1; j-dj >=0; j-=dj ){
		w1 = wp[j];
		if(w2 > w1) { continue; }
		else{
		//cout<<w1<<endl;
		r1 = r[j];
		double qmax = sqrt(2*me*Ts[j])/20;
		energy.push_back(2*w1);
		energy1.push_back(2*w1);
		double int1 = lowq_integral(j,qmax) * abs((r2-r1)/(w2-w1))/2;
		double int2 = r_integral3(2*w1,0);
		cout<<"w = "<<2*w1<<"	int1 = "<<int1<<"	int2 = "<<int2<<endl;
		count1.push_back( int1 );
		count2.push_back( int2 );
		int1 = 0;
		count.push_back( int1 + int2 );
		r2 = r[j];
		w2 = wp[j];
		}
	}
	for( double w = 2*w2; w<=2e4; w+=w2) {
		energy.push_back(w);
		double int2 = r_integral3(w,0);
		cout<<"w = "<<w<<"	int2 = "<<int2<<endl;
		count.push_back(int2);
		count2.push_back(int2);
	}
	// write to file
	string name = "data/LLspectrum_tot-100.dat";
	string name1 = "data/LLspectrum_1-100.dat";
	string name2 = "data/LLspectrum_2-100.dat";
	write2D( name , energy, count );
	write2D( name1 , energy1, count1 );
	write2D( name2 , energy, count2 );
}

void fixedR( int cs ) {
	vector<double> count, energy;
	string name;
	double dw = 1e1;
	for( double w = dw; w < 5e3; w+=dw ){
		energy.push_back(w);					// eV
		count.push_back( q_integral3(w,cs,1000) );	// Bg-2
		if((int)(w) % (int)(1e2) == 0) { cout<<"w = "<<w/1e3<<"keV of 20keV"<<endl; }
		//cout<<"w = "<<w<<"eV of 20keV"<<endl;
	}
	name = "data/LLspectrum_half.dat";
	write2D( name , energy, count );	
}

//void fullR() {
//	vector<double> count, energy;
//	string name;
//	double dw = 1e1;
//	for( double w = dw; w < 2e3; w+=dw ){
//		energy.push_back(w);			
//		count.push_back( r_integral3(w) );
//		//if((int)(w) % (int)(1e2) == 0) { cout<<"w = "<<w/1e3<<"keV of 20keV"<<endl; }
//		cout<<"w = "<<w/1e3<<"keV of 20keV"<<endl;
//	}
//	name = "data/LLspectrum--1e-6.dat";
//	write2D( name , energy, count );	
//}


int main() {
	//double I = E_integral3(20,0);
	//args qargs;
	//double w = 300;
	//qargs.k = w;
	//qargs.w = w;
	//qargs.cs = 0;
	//qargs.q = 1000;
	//double I = E_integral3_gsl(qargs.q,&qargs);
	//double I = r_integral3(1000);
	//double I = 	q_integral3(300,0,10);
	//double I = lowq_integral(0,sqrt(2*me*Ts[0])/10);
	//cout<<I<<endl;
	//double J = r_integral3(2*wp[0],sqrt(2*me*Ts[0])/10);
	//cout<<J<<endl;
	//fixedR(998);	// r = 0.5 rSolar
	spectrum_2parts();
	return 0;
}
