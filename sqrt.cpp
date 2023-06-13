// Tom O'Shea 2023

// complex sqrt

#include "chameleon.h"	// home of the real code

using namespace std;

int main( int argc, char** argv ){

 //   if ( argc == 3 ) {
        //double x = argv[1];
       // double y = argv[2];
       double x = 0;
       double y = -8;
       double a,b;
        csqrt(x,y,&a,&b);
        cout << a << " + " << b << "i" << endl;

        return 0;
 //   }

 //   else { cout << "invalid usage! (length " << argc << ")\nusage: " << argv[0] << " <real part> <imaginary part>" << endl; return 69;}
}