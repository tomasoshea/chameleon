// Tom O'Shea 2023

// get values to plot meff against solar radius for chameleons of various Bm or n

#include "chameleon.h"	// home of the real code

using namespace std;

int main(){

	//for( double Bm = 1e0 ; Bm < 2e6 ; Bm *= 10 ) { plotMeff( Bm ); }
	//for( double n = 1 ; n < 5 ; n += 1. ) { plotMeff( n ); }

	plotMeff(4.);
	
	return(0);
}