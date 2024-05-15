using namespace std;
#include "Domain.h"
#include "Constants.h"
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <iostream>


Domain::Domain()
{
	_Natoms=0;
	_Nval=1;
	_N1=0;
	_N2=0;
	_N3=0;
	_T.resize(3, vector<double>(3));
	
	for(int i=0;i<3;i++)
	{
		_O[i]=0;
		for(int j=0; j<3;j++)
		{
			_T[i][j]=0;
		}
	}
}

void Domain::read_From_Cube(ifstream& nameFile)
{
	string bin;
	for(int i=0;i<12;i++) {nameFile>>bin;}
	nameFile>>_Natoms;
	cout<<_Natoms<<endl;
	_O.resize(3);
	for(int i=0;i<3;i++)
	{
		nameFile>>_O[i];
		cout<<_O[i];
	}
	nameFile>>_Nval;
	nameFile>>_N1;
	_T.resize(3, vector<double>(3));
	if(_N1<0)
	{
		for(int i=0; i<3; i++)
		{
			nameFile>>_T[1][i];
		}
		nameFile>>_N2;
		for(int i=0; i<3; i++)
		{
			nameFile>>_T[2][i];
		}
		nameFile>>_N3;
		for(int i=0; i<3; i++)
		{
			nameFile>>_T[3][i];
		}
		nameFile>>_Natoms;
		_N1=abs(_N1);
	}
	
	else
	{
		for(int i=0; i<3; i++)
		{
			nameFile>>_T[1][i];
			_T[1][i]=_T[1][i]*ANGTOBOHR;
		}
		nameFile>>_N2;
		for(int i=0; i<3; i++)
		{
			nameFile>>_T[2][i];
			_T[1][i]=_T[1][i]*ANGTOBOHR;
		}
		nameFile>>_N3;
		for(int i=0; i<3; i++)
		{
			nameFile>>_T[3][i];
			_T[1][i]=_T[1][i]*ANGTOBOHR;
		}
		nameFile>>_Natoms;
	}
}

Domain::Domain(ifstream& nameFile)
{
	read_From_Cube(nameFile);	
}

int Domain::Natoms() const
{
	return _Natoms;
}

int Domain::Nval() const
{
	return _Nval;
}

int Domain::N1() const
{
	return _N1;
}

int Domain::N2() const
{
	return _N2;
}

int Domain::N3() const
{
	return _N3;
}

vector<double> Domain::O() const
{
	return _O;
}

vector<vector<double>> Domain::T() const
{
	return _T;
}


