#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <math.h> 
#include <filesystem>
#include <bits/stdc++.h>

using namespace std;

// conversion factors
double s2eV = (6.582119569e-16);	// Hz to eV
double J2eV = (1. / 1.602176634e-19);	// Joules to eV (1 / e)
double m2eV = (1.973269804e-7);	// m-1 to eV
double K2eV = (8.617333262e-5);	// Kelvin to eV
double kg2eV = 5.609588604e35;	// from hbar/c2
double T2eV = 2e-16 * 1e18;		// Tesla to eV2 conversion [eV2/T]


// read in datafiles in csv form
vector<double> read( string name ) {

	//cout << "Reading file " << name << "..." << endl;

	// open file defined in argument
	fstream file;
	file.open(name,ios::in);
	
	char delim('\n');	// define delimiter for file parsing
	
	if (file.is_open()){   // checking whether the file is open
		string temp;	// define temporary storage string
		vector<double> row;	// define vector to store input values and return
		
		while(getline(file, temp, delim)){  // read data from file object and put it into string
			double item = stod(temp);	// convert string to double
			row.push_back(item);	// append double to vector
		}
		
	file.close();   // close the file object.
	
	return row;	// returns vector of values
	
	}
	else{ cout << "some error i guess..." << endl ; vector<double> err = {69.} ; return err ; }
}


// write out 2 column datafile
void write2D( string name, vector<double> data1, vector<double> data2) {

	// delete old file if extant
	if ( remove(name.c_str()) == 0 ) {
		cout << "file " << name << " deleted." << endl;
		}
	
	// create new file
	fstream fout;
	fout.open(name, ios::out | ios::trunc);
	
	// get size of vectors
	int len1 = data1.size();
	int len2 = data2.size();
	
	// check vectors of equal length
	if ( len1 == len2 ) {
		
		// Read the input from vector
		for ( int c = 0; c < len1; c++ ) {
			// Insert the data to file
			fout << data1[c] << "	" << data2[c] << endl;
		}
		
		cout << "file " << name << " created succesfully." << endl;
	}
	
	// if vectors not equal return error
	else { cout << "ERROR - ensure vectors are of equal length" << endl; }
	
	fout.close();
}

// read in gaunt factors from matlab matrix files
vector<vector<double>> readGaunt( string name ) {

	//cout << "Reading file " << name << "..." << endl;

	// open file defined in argument
	fstream file;
	file.open(name,ios::in);
	
	char delim(' ');	// define delimiter for file parsing
	
	if (file.is_open()){   // checking whether the file is open
		string temp;	// define temporary storage string
		vector<vector<double>> g2D; // 2D matrix to store rows
		vector<double> row;	// define vector to store input values and return
		vector<double> x;	// define vector to store input values and return
		
		int c = 0;	// define counter for elements
		int line = 0;	// counter for lines
		
		while( !file.eof() ) {  // loop until end of file
		
			getline(file, temp, delim);	// get one row
			
			// check data is not empty
			if( temp == "\n" ) { continue; }
			
			double item = stod(temp);	// convert string to double

			// add a zero on first line
			if( c == 0 ) { row.push_back(0.); }
			
			row.push_back(item);	// append double to vector
			c++;
			
			// when row is full append to 2D vector and wipe row				
			if( row.size() == 201 ) {
			
				g2D.push_back(row);
				row.clear();
				
				line++;
			}

			if( line == 501 ) { break; }

		}
		
	file.close();   // close the file object.
	
	return g2D;	// returns vector of values
	}
	
	else{ cout << "couldn't find file: " << name << endl ; vector<vector<double>> err = {{69.}} ; return err ; }
}

// write big matrices
void writeREST( string name, vector<double> data, int line ){	// data vector needs values of w for constant r

	//int piece = 0;	// for newlining
	int len = data.size();
	//cout<<len<<endl;

	// delete old file if extant
	if ( line == 0 ) {
	if ( remove(name.c_str()) == 0 ) {
		cout << "file " << name << " deleted." << endl;
		}
	}

	// file pointer
	fstream fout;

	// creates new csv file
	fout.open(name, ios::out | ios::app);

	// Read the input from vector
	int c = 0;
	for ( double item : data ) {
		// Insert the data to file
		//cout<<item<<endl;
		fout << item;
		if(c == len-1) { fout<<endl; break; }
		else { fout << "	"; c++; }
	}
	if( line == 0 ){ cout << "writing to file " << name << "...\n" << endl; }
	
	fout.close();
}
