#include<iostream>
#include <string>
#include <fstream>
#include <vector>
#include<common/Structure.h>

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
	//cout<<"Natoms="<<Natoms<<endl;
	for(int i=0; i<Natoms; i++)
	{	
		nameFile>>input;
		//cout<<"done"<<endl;
		Atom a(Table, int(input));
		//cout<<"done"<<endl;
		_atoms[i]=a;
		//cout<<"done"<<endl;
		nameFile>>input;
		//cout<<"done"<<endl;
		_atoms[i].Set_charge(input);
		//cout<<"done"<<endl;
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

void Structure::read_from_wfx(WFX& wfx, const PeriodicTable& Table)
{
	int n=0;
	_atoms.resize(wfx.Number_of_Nuclei());
	for(int i=0; i<wfx.Number_of_Nuclei(); i++)
	{
		_atoms[i]=Atom(Table, wfx.Atomic_Number()[i]);
		_atoms[i].Set_charge(wfx.Nuclear_Charges()[i]);
		n=3*i;
		for(int j=0; j<3; j++)
		{
			_atoms[i].Set_coordinates(j, wfx.Nuclear_Cartesian_Coordinates()[n]);
			n++;
		}
	}
}

Structure::Structure(WFX& wfx, const PeriodicTable& Table)
{
	read_from_wfx(wfx, Table);
}