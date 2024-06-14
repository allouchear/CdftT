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
	_atoms=vector<Atom> ();
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

<<<<<<< HEAD
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




















=======
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
>>>>>>> dbuffat
