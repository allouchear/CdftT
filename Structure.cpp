#include<iostream>
#include <string>
#include <fstream>
#include"Structure.h"

using namespace std;

Structure::Structure(const vector<Atom>& A)
{
	_atoms.resize(A.size());
	 _atoms=A;
}

Structure::Structure()
{
	_atoms.resize(0);
	_atoms[0]=Atom();
}

void Structure::cube(ifstream& nameFile)
{
	string bin;
	getline(nameFile, bin);
	getline(nameFile, bin);
	double input;
	nameFile>>input;
	_atoms.resize(input);
	//read origin maybe
	//read Nval
	//read volmetric geometry
	for(int i=0; i<_atoms.size(); i++)
	{	
		nameFile>>input;
		
	}
}
