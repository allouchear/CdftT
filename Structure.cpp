#include<iostream>
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