// Tom O'Shea 2023
// get solar chameleon spectrum (tachocline)

#include "cham2.cpp"

using namespace std;

int main() {

    // choose model params
    double n = 1;
    double Bm = 1e1;

    // initialise vectors
	vector<double> flux, omega;
    for ( double w = 1.; w < 3e3; w++ ) {
		flux.push_back(solarFlux(w,n,Bm));	// [m-2 s-1 eV-1]
		omega.push_back(w);		// [eV]
		//cout << "omega = " << w << endl;
	}

	write2D( "data/chamB1-spectrum.dat", omega, flux );

	return 0;
}