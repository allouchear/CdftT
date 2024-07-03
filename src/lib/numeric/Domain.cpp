using namespace std;
#include <numeric/Domain.h>
#include <common/Constants.h>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <iostream>
#include <sstream>


Domain::Domain()
{
	_Nval=1;
	_N1=0;
	_N2=0;
	_N3=0;
	_T.resize(3, vector<double>(3));
	_inv_T.resize(3, vector<double>(3));
	for(int i=0;i<3;i++)
	{
		_O[i]=0;
		for(int j=0; j<3;j++)
		{
			_T[i][j]=0;
		}
	}
	inverse_T();
	_dx=0;
	_dy=0;
	_dz=0;
	_dv=0;
}

void Domain::read_From_Cube(ifstream& nameFile)
{
	string s;
	getline(nameFile,s);
	stringstream ss(s);
	// Tokenize the input string by comma delimiter

	for(int i=0;i<3;i++)
	{
		ss>>_O[i];
	}
	_Nval = 1;
	ss>>_Nval;
	if(ss.fail()) _Nval=1;	
	nameFile>>_N1;
	_T.resize(3, vector<double>(3));
	if(_N1>0)
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
		for(int i=0;i<3;i++)
		{
			_O[i]=_O[i]*ANGTOBOHR;
		}
	}
	for(int i=0; i<3;i++)
	{
		_dx += _T[0][i]*_T[0][i];
		_dy += _T[1][i]*_T[1][i];
		_dz += _T[2][i]*_T[2][i];
	}
	_inv_T.resize(3, vector<double>(3));
	inverse_T();
	_dx = sqrt(_dx);
	_dy = sqrt(_dy);
	_dz = sqrt(_dz);
	_dv = _dx*_dy*_dz;
}

Domain::Domain(ifstream& nameFile)
{
	read_From_Cube(nameFile);
}

Domain::Domain(int i, int n, int m, int l, double* O)
{
	
	{
	set_Nval(i);
	set_N1(n);
	set_N2(m);
	set_N3(l);
	_T.resize(3, vector<double>(3));
	if(O==NULL)
	{
		for(int i=0;i<3;i++)
		{
			_O[i]=0;
			for(int j=0; j<3;j++)
			{
				_T[i][j]=0;
			}
		}
	}
	else
	{
		for(int i=0;i<3;i++)
		{
			_O[i]=O[i];
			for(int j=0; j<3;j++)
			{
				_T[i][j]=0;
			}
		}
	}
	for(int i=0; i<3;i++)
	{
		_dx += _T[0][i]*_T[0][i];
		_dy += _T[1][i]*_T[1][i];
		_dz += _T[2][i]*_T[2][i];
	}
	_inv_T.resize(3, vector<double>(3));
	inverse_T();
	_dx = sqrt(_dx);
	_dy = sqrt(_dy);
	_dz = sqrt(_dz);
	_dv = _dx*_dy*_dz;
	}
}

int Domain::Nval() const
{
	return _Nval;
}

void Domain::set_Nval(int N)
{
	_Nval=N;
}

int Domain::N1() const
{
	return _N1;
}

void Domain::set_N1(int N)
{
	_N1=N;
}

int Domain::N2() const
{
	return _N2;
}

void Domain::set_N2(int N)
{
	_N2=N;
}

int Domain::N3() const
{
	return _N3;
}

void Domain::set_N3(int N)
{
	_N3=N;
}

double* Domain::O()
{
	return _O;
}

vector<vector<double>> Domain::T() const
{
	return _T;
}

double Domain::Tij(int i , int j) const
{
	return _T[i][j];
}

void Domain::set_T(double v, int i, int j)
{
	_T[i][j]=v;
}

double Domain::dx() const
{
	return _dx;
}

double Domain::dy() const
{
	return _dy;
}

double Domain::dz() const
{
	return _dz;
}
double Domain::x(int i, int j, int k) const
{
	return _O[0]+_T[0][0]*i+_T[0][1]*j+_T[0][2]*k;
}
double Domain::y(int i, int j, int k) const
{
	return _O[1]+_T[1][0]*i+_T[1][1]*j+_T[1][2]*k;
}
double Domain::z(int i, int j, int k) const
{
	return _O[2]+_T[2][0]*i+_T[2][1]*j+_T[2][2]*k;
}

int Domain::i(double x, double y, double z) const
{
	return floor((x-_O[0])*_inv_T[0][0]+(y-_O[1])*_inv_T[0][1]+(z-_O[2])*_inv_T[0][2]);
}

int Domain::j(double x, double y, double z) const
{
	return floor((x-_O[0])*_inv_T[1][0]+(y-_O[1])*_inv_T[1][1]+(z-_O[2])*_inv_T[1][2]);
}

int Domain::k(double x, double y, double z) const
{
	return floor((x-_O[0])*_inv_T[2][0]+(y-_O[1])*_inv_T[2][1]+(z-_O[2])*_inv_T[2][2]);
}

double Domain::dv() const
{
	return _dv;
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

void Domain::inverse_T()
{
	double t4,t6,t8,t10,t12,t14,t17;
	
	t4 = _T[0][0]*_T[1][1];     
	t6 = _T[0][0]*_T[1][2];
	t8 = _T[0][1]*_T[1][0];
	t10 = _T[0][2]*_T[1][0];
	t12 = _T[0][1]*_T[2][0];
	t14 = _T[0][2]*_T[2][0];
	t17 = (t4*_T[2][2]-t6*_T[2][1]-t8*_T[2][2]+t10*_T[2][1]+t12*_T[1][2]-t14*_T[1][1]);
	if(abs(t17)<1e-12)
	{
		for(int i=0;i<3;i++)
		{
			for(int j=0; j<3;j++)
			{
				if(i==j)
				{
					_inv_T[i][j]=1;	
				}
				else
				{
					_inv_T[i][j]=0;
				}
			}
		}
	}
	else
	{
		t17 = 1/t17;
		_inv_T[0][0] = (_T[1][1]*_T[2][2]-_T[1][2]*_T[2][1])*t17;
		_inv_T[0][1] = -(_T[0][1]*_T[2][2]-_T[0][2]*_T[2][1])*t17;
		_inv_T[0][2] = -(-_T[0][1]*_T[1][2]+_T[0][2]*_T[1][1])*t17;
		_inv_T[1][0] = -(_T[1][0]*_T[2][2]-_T[1][2]*_T[2][0])*t17;
		_inv_T[1][1] = (_T[0][0]*_T[2][2]-t14)*t17;
		_inv_T[1][2] = -(t6-t10)*t17;
		_inv_T[2][0] = -(-_T[1][0]*_T[2][1]+_T[1][1]*_T[2][0])*t17;
		_inv_T[2][1] = -(_T[0][0]*_T[2][1]-t12)*t17;
		_inv_T[2][2] = (t4-t8)*t17;
	}
}
