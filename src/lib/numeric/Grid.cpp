using namespace std;
#include <numeric/Grid.h>
#include <common/Constants.h>
#include <cmath>


Grid::Grid()
{
	Domain d;
	_dom=d;
	Structure S;
	_str=S;
	_V.resize(_dom.N1(),vector<vector<vector<double>>>(_dom.N2(),vector<vector<double>>(_dom.N3(), vector<double>(_dom.Nval(), 0))));
}

void Grid::read_From_Cube(ifstream& nameFile, const PeriodicTable& Table)
{
	int Natoms;
	string bin;
	getline(nameFile, bin);
	getline(nameFile, bin);
	nameFile>>Natoms;
	Domain d(nameFile);
	_dom=d;
	Structure S(nameFile, Natoms, Table);
	_str=S;
	_V.resize(_dom.N1());
	for(int i=0; i<_dom.N1();i++)
	{	
		_V[i].resize(_dom.N2());
		for(int j=0; j<_dom.N2();j++)
		{
			_V[i][j].resize(_dom.N3());
			for(int k=0; k<_dom.N3();k++)
			{
				_V[i][j][k].resize(_dom.Nval());
				for(int l=0; l<_dom.Nval();l++)
				{
				nameFile>>_V[i][j][k][l];
				}
			}
		}
	}
}

Grid::Grid(ifstream& nameFile, const PeriodicTable& Table)
{
	read_From_Cube(nameFile, Table);	
}

vector<vector<vector<vector<double>>>> Grid::V() const
{
	return _V;
}

Domain Grid::dom() const
{
	return _dom;
}

Structure Grid::str() const
{
	return _str;
}

void Grid::set_dom(const Domain& d)
{
	_dom=d;
}

void Grid::set_str(const Structure& S)
{
	_str=S;
}

void Grid::set_V(const vector<vector<vector<vector<double>>>>& U)
{
	_V=U;
}

Grid Grid::operator+(const Grid& g)
{
	try
	{
		if(g._dom==_dom)
		{
			Grid sum;			
			sum.set_dom(_dom);
			sum.set_str(_str+g._str);
			sum._V.resize(g._dom.N1());
			for(int i=0; i<g._dom.N1();i++)
			{	
				sum._V[i].resize(g._dom.N2());
				for(int j=0; j<g._dom.N2();j++)
				{
					sum._V[i][j].resize(g._dom.N3());
					for(int k=0; k<g._dom.N3();k++)
					{
						sum._V[i][j][k].resize(g._dom.Nval());
						for(int l=0; l<g._dom.Nval();l++)
						{
						sum._V[i][j][k][l]=_V[i][j][k][l]+g._V[i][j][k][l];
						}
					}
				}
			}
			return sum;
		}
		else if(g.dom().Nval()!=dom().Nval())
		{
			throw string("Attribute (grid).(domain)._Nval aren't equal in both grids. This is required for addition");
		}
		else if (g.dom().N1()!=dom().N1())
		{
			throw string("Attributes (grid).(domain)._N1 aren't equal in both grids. This is required for addition");
		}
		else if (g.dom().N2()!=dom().N2())
		{
			throw string("Attributes (grid).(domain)._N2 aren't equal in both grids. This is required for addition");
		}
		else if(g.dom().N3()!=dom().N3())
		{
			throw string("Attributes (grid).(domain)._N3 aren't equal in both grids. This is required for addition");
		}
		else
		{
			throw string("Attributes (grid).(domain)._T aren't equal in both grids. This is required for addition.");
		}
	}
	catch (string Error)
	{
		cout<<Error<<endl;
		exit(1);
	}
}

Grid Grid::add(const Grid& g)
{
	try
	{
		if(g._dom==_dom)
		{
			set_dom(_dom);
			_str.add(g._str);
			_V.resize(g._dom.N1());
			for(int i=0; i<g._dom.N1();i++)
			{	
				_V[i].resize(g._dom.N2());
				for(int j=0; j<g._dom.N2();j++)
				{
					_V[i][j].resize(g._dom.N3());
					for(int k=0; k<g._dom.N3();k++)
					{
						_V[i][j][k].resize(g._dom.Nval());
						for(int l=0; l<g._dom.Nval();l++)
						{
						_V[i][j][k][l]=_V[i][j][k][l]+g._V[i][j][k][l];
						}
					}
				}
			}
			return *this;
		}
		else if(g.dom().Nval()!=dom().Nval())
		{
			throw string("Attribute (grid).(domain)._Nval aren't equal in both grids. This is required for addition");
		}
		else if (g.dom().N1()!=dom().N1())
		{
			throw string("Attributes (grid).(domain)._N1 aren't equal in both grids. This is required for addition");
		}
		else if (g.dom().N2()!=dom().N2())
		{
			throw string("Attributes (grid).(domain)._N2 aren't equal in both grids. This is required for addition");
		}
		else if(g.dom().N3()!=dom().N3())
		{
			throw string("Attributes (grid).(domain)._N3 aren't equal in both grids. This is required for addition");
		}
		else
		{
			throw string("Attributes (grid).(domain)._T aren't equal in both grids. This is required for addition.");
		}
	}
	catch (string Error)
	{
		cout<<Error<<endl;
		exit(1);
	}
}

Grid Grid::operator*(const Grid& g)
{
	try
	{
		if(g._dom==_dom)
		{
			Grid prod;
			prod.set_dom(_dom);
			prod.set_str(_str);
			prod._V.resize(g._dom.N1());
			for(int i=0; i<g._dom.N1();i++)
			{	
				prod._V[i].resize(g._dom.N2());
				for(int j=0; j<g._dom.N2();j++)
				{
					prod._V[i][j].resize(g._dom.N3());
					for(int k=0; k<g._dom.N3();k++)
					{
						prod._V[i][j][k].resize(g._dom.Nval());
						for(int l=0; l<g._dom.Nval();l++)
						{
						prod._V[i][j][k][l]=_V[i][j][k][l]*g._V[i][j][k][l];
						}
					}
				}
			}
			return prod;
		}
		else if(g.dom().Nval()!=dom().Nval())
		{
			throw string("Attribute (grid).(domain)._Nval aren't equal in both grids. This is required for product");
		}
		else if (g.dom().N1()!=dom().N1())
		{
			throw string("Attributes (grid).(domain)._N1 aren't equal in both grids. This is required for product");
		}
		else if (g.dom().N2()!=dom().N2())
		{
			throw string("Attributes (grid).(domain)._N2 aren't equal in both grids. This is required for product");
		}
		else if(g.dom().N3()!=dom().N3())
		{
			throw string("Attributes (grid).(domain)._N3 aren't equal in both grids. This is required for product");
		}
		else
		{
			throw string("Attributes (grid).(domain)._T aren't equal in both grids. This is required for product");
		}
	}
	catch (string Error)
	{
		cout<<Error<<endl;
		exit(1);
	}
}

Grid Grid::operator-(const Grid& g)
{
	try
	{
		if(g._dom==_dom)
		{
			Grid diff;
			diff.set_dom(_dom);
			diff.set_str(_str);
			diff._V.resize(g._dom.N1());
			for(int i=0; i<g._dom.N1();i++)
			{	
				diff._V[i].resize(g._dom.N2());
				for(int j=0; j<g._dom.N2();j++)
				{
					diff._V[i][j].resize(g._dom.N3());
					for(int k=0; k<g._dom.N3();k++)
					{
						diff._V[i][j][k].resize(g._dom.Nval());
						for(int l=0; l<g._dom.Nval();l++)
						{
						diff._V[i][j][k][l]=_V[i][j][k][l]-g._V[i][j][k][l];
						}
					}
				}
			}
			return diff;
		}
		else if(g.dom().Nval()!=dom().Nval())
		{
			throw string("Attribute (grid).(domain)._Nval aren't equal in both grids. This is required for difference");
		}
		else if (g.dom().N1()!=dom().N1())
		{
			throw string("Attributes (grid).(domain)._N1 aren't equal in both grids. This is required for difference");
		}
		else if (g.dom().N2()!=dom().N2())
		{
			throw string("Attributes (grid).(domain)._N2 aren't equal in both grids. This is required for difference");
		}
		else if(g.dom().N3()!=dom().N3())
		{
			throw string("Attributes (grid).(domain)._N3 aren't equal in both grids. This is required for difference");
		}
		else
		{
			throw string("Attributes (grid).(domain)._T aren't equal in both grids. This is required for difference");
		}
	}
	catch (string Error)
	{
		cout<<Error<<endl;
		exit(1);
	}
}

double Grid::sum()
{
	double sum=0;
	for(int i=0; i<_dom.N1();i++)
	{	
		for(int j=0; j<_dom.N2();j++)
		{	
			for(int k=0; k<_dom.N3();k++)
			{
				for(int l=0; l<_dom.Nval();l++)
				{
					sum+=_V[i][j][k][l];
				}
			}
		}
	}
	return sum;	
}

double Grid::integrate_over_dom()
{
	double I=0;
	double dx=0;
	double dy=0;
	double dz=0;
	for(int i=0; i<3;i++)
	{
		dx += _dom.T()[0][i]*_dom.T()[0][i];
		dy += _dom.T()[1][i]*_dom.T()[1][i];
		dz += _dom.T()[2][i]*_dom.T()[2][i];
	}
	I = sum()*sqrt(dx*dy*dz);
	return I;
}
/*
Grid Grid::resize(int n, int m, int p)
{
	_V.erase(_dom.N1()-n, _dom.N1());
	for(int i=0; i<_dom.N1();i++)
	{	
		_V[i].erase(_dom.N2()-m, _dom.N2());
		for(int j=0; j<_dom.N2();j++)
		{	
			_V[i][j].erase(_dom.N3()-p, _dom.N3());
		}
	}
	return *this;
}
*/
Grid Grid::resize_zeros(int n, int m, int p)
{
	for(int i=_dom.N1()-n; i<_dom.N1();i++)
	{	
		for(int j=_dom.N2()-m; j<_dom.N2();j++)
		{	
			for(int k=_dom.N3()-p; k<_dom.N3();k++)
			{
				for(int l=0; l<_dom.Nval();l++)
				{
					_V[i][j][k][l]= 0;
				}
			}
		}
	}
	return *this;
}

void Grid::set_V_Func(vector<double> fpar, function<double(vector<double>,double, double, double, const Grid& g )> f, const Grid& g)
{
	set_dom(g.dom());
	int N=_dom.N1();
	int M=_dom.N2();
	int P=_dom.N3();
	if(int(fpar[1])==1)
	{
		N= N-fpar[0];
		M= M-fpar[0];
		P= P-fpar[0];
	}
	double x = _dom.O()[ 0 ];
	double y = _dom.O()[ 1 ];
	double z = _dom.O()[ 2 ];
	_V.resize(g._dom.N1());
	for(int i=0; i<N;i++)
	{	
		_V[i].resize(g._dom.N2());
		for(int j=0; j<M;j++)
		{	
			_V[i][j].resize(g._dom.N3());
			for(int k=0; k<P;k++)
			{
				_V[i][j][k].resize(g._dom.Nval());
				for(int l=0; l<_dom.Nval();l++)
				{
					x = _dom.O()[ 0 ] + i*_dom.T()[0][0] + j*_dom.T()[0][1] +  k*_dom.T()[0][2]; 
					y = _dom.O()[ 1 ] + i*_dom.T()[1][0] + j*_dom.T()[1][1] +  k*_dom.T()[1][2]; 
					z = _dom.O()[ 2 ] + i*_dom.T()[2][0] + j*_dom.T()[2][1] +  k*_dom.T()[2][2];
					
					_V[i][j][k][l]= f(fpar,x,y,z, g);
				}
			}
		}
	}
	if(int(fpar[1])==1)
	{	
		resize_zeros(int(fpar[0]),int(fpar[0]),int(fpar[0]));
	}	
}

