// Tom O'Shea 2023

// main chameleon code for oututting flux file for limit setting

#include "cham2.h"

using namespace std;

int main(){

    // set model parameter n
    double n = 100;
    //for( double n = 1; n <= 10; n++ ) {

    // run over each

    // babyIAXO
    double L = 10;  // babyIAXO bore length [m]
    double B = 2 * T2eV;   // babyIAXO B-field [eV2]
    thread t1(calc, n, L, B, "babyIAXO");

    // baseline IAXO
    L = 20;
    B = 2.5 * T2eV;
    thread t2(calc, n, L, B, "baselineIAXO");

    // upgraded IAXO
    L = 22;
    B = 3.5 * T2eV;
    thread t3(calc, n, L, B, "upgradedIAXO");

    //cout << pow( B * (L/m2eV) / mpl , 2 ) << endl;

    t1.join();
    t2.join();
    t3.join();

   // }

    return 0;
}