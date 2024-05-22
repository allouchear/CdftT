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
	for(int i=0; i<Natoms; i++)
	{	
		nameFile>>input;
		Atom a(Table, int(input));
		_atoms[i]=a;
		nameFile>>input;
		_atoms[i].Set_charge(input);
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

Structure Structure::operator+(const Structure& S) const
{
	if(S._atoms.size()>_atoms.size())
	{
		return S;
	}
	else
	{
		return *this;
	}
}

Structure Structure::add(const Structure& S)
{
	try
	{
		for(int j=0; j<S.number_of_atoms(); j++)
		{	
			bool b=true;
			for(int k=0; k<number_of_atoms(); k++)
			{
				if( _atoms[k].get_distance(S._atoms[j])==0 )
				{
					if(_atoms[k].symbol()==S._atoms[j].symbol())
					{
						throw string("::add(const Structure& S): can't add two identical atoms at same coordinates");
					}
					else
					{
						b=false;
						break;
					}
				}
			}
			if(b)
			{
					cout<<"added"<<endl;
					_atoms.push_back(S._atoms[j]);
			}
		}
		return *this;
	}
	catch(string error)
	{
		cout<<error<<endl;
		exit(1);
	}
}




















