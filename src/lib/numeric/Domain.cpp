using namespace std;
#include <numeric/Domain.h>
#include <common/Constants.h>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <iostream>


Domain::Domain()
{
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
	for(int i=0;i<3;i++)
	{
		nameFile>>_O[i];
	}
	nameFile>>_Nval;
	nameFile>>_N1;
	_T.resize(3, vector<double>(3));
	if(_N1<0)
	{
		for(int i=0; i<3; i++)
		{
			nameFile>>_T[0][i];
		}
		nameFile>>_N2;
		for(int i=0; i<3; i++)
		{
			nameFile>>_T[1][i];
		}
		nameFile>>_N3;
		for(int i=0; i<3; i++)
		{
			nameFile>>_T[2][i];
		}
		_N1=abs(_N1);
	}
	
	else
	{
		for(int i=0; i<3; i++)
		{
			nameFile>>_T[0][i];
			_T[0][i]=_T[0][i]*ANGTOBOHR;
		}
		nameFile>>_N2;
		for(int i=0; i<3; i++)
		{
			nameFile>>_T[1][i];
			_T[1][i]=_T[1][i]*ANGTOBOHR;
		}
		nameFile>>_N3;
		for(int i=0; i<3; i++)
		{
			nameFile>>_T[2][i];
			_T[2][i]=_T[2][i]*ANGTOBOHR;
		}
	}
}

Domain::Domain(ifstream& nameFile)
{
	read_From_Cube(nameFile);	
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

double* Domain::O()
{
	return _O;
}

vector<vector<double>> Domain::T() const
{
	return _T;
}

bool Domain::operator==(const Domain& D) const
{
	if(D._Nval==_Nval and D._N1==_N1 and D._N2==_N2 and D._N3==_N3 and D._T==_T)
		return true;
	else
		return false;
}

bool Domain::operator!=(const Domain& D) const
{
	if( D._Nval!=_Nval or D._N1!=_N1 or D._N2!=_N2 or D._N3!=_N3 or D._T!=_T)
		return true;
	else
		return false;
}
