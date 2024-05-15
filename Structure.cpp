#include<iostream>
#include <string>
#include <fstream>
#include <vector>
#include"Structure.h"

using namespace std;

Structure::Structure(const vector<Atom>& A)
{
	_atoms.resize(A.size());
	_atoms=A;
}

Structure::Structure()
{
	_atoms.resize(1);
	_atoms[0]=Atom();
}

void Structure::read_From_Cube(ifstream& nameFile, int Natoms,const PeriodicTable& Table )
{
	double input;
	_atoms.resize(Natoms);
	cout<<"Natoms="<<Natoms<<endl;
	for(int i=0; i<Natoms; i++)
	{	
		nameFile>>input;
		cout<<"done"<<endl;
		Atom a(Table, int(input));
		cout<<"done"<<endl;
		_atoms[i]=a;
		cout<<"done"<<endl;
		nameFile>>input;
		cout<<"done"<<endl;
		_atoms[i].Set_charge(input);
		cout<<"done"<<endl;
		for(int j=0; j<3; j++)
		{
			nameFile>>input;
			_atoms[i].Set_coordinates(j,input);
		}
	}
}

Structure::Structure(ifstream& nameFile, const int Natoms, const PeriodicTable& Table)
{
	read_From_Cube(nameFile, Natoms, Table);	
}

