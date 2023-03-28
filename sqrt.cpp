// Tom O'Shea 2023

// complex sqrt

#include "chameleon.h"	// home of the real code

using namespace std;

int main(){

    double a,b;
    csqrt(0,-200,1,1,&a,&b);
    cout << a << " + " << b << "i" << endl;
    
	return(0);
}