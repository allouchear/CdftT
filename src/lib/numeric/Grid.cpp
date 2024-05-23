using namespace std;
#include <numeric/Grid.h>
#include <common/Constants.h>
#include <cmath>
#ifdef ENABLE_OMP
#include <omp.h>
#endif


void Grid::reset()
{
	if (_dom.Nval()<1 || _dom.N1()<1 || _dom.N2()<1 || _dom.N3()<1 )
		_V = vector<vector<vector<vector<double>>>>();
	else 
	{
		vector<double> V(_dom.Nval(),0);
		vector<vector<double>> M(_dom.N3(), V);
		vector<vector<vector<double>>> A(_dom.N2(), M);
		_V = vector<vector<vector<vector<double>>>>(_dom.N1(), A);
	}
}

Grid::Grid()
{
	Domain d;
	_dom=d;
	Structure S;
	_str=S;
	reset();
}
Grid::Grid(const Domain& d)
{
	_dom=d;
	Structure S;
	_str=S;
	reset();
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
	cout<<"done struct"<<endl;
	reset();
	cout<<"done reset"<<endl;
	for(int i=0; i<_dom.N1();i++)
	{
		for(int j=0; j<_dom.N2();j++)
		{
			for(int k=0; k<_dom.N3();k++)
			{
				for(int l=0; l<_dom.Nval();l++)
				{
				nameFile>>_V[i][j][k][l];
				}
			}
		}
	}
	cout<<"done read"<<endl;
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
	reset();
	
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
			Grid sum(_dom);
			sum.set_str(_str+g._str);
#ifdef ENABLE_OMP
#pragma omp parallel for
#endif
			for(int i=0; i<g._dom.N1();i++)
			{
				for(int j=0; j<g._dom.N2();j++)
				{
					for(int k=0; k<g._dom.N3();k++)
					{
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
#ifdef ENABLE_OMP
#pragma omp parallel for
#endif
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
			Grid prod(_dom);
			prod.set_str(_str);
#ifdef ENABLE_OMP
#pragma omp parallel for
#endif
			for(int i=0; i<g._dom.N1();i++)
			{
				for(int j=0; j<g._dom.N2();j++)
				{
					for(int k=0; k<g._dom.N3();k++)
					{
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
			Grid diff(_dom);
			diff.set_str(_str);
#ifdef ENABLE_OMP
#pragma omp parallel for
#endif
			for(int i=0; i<g._dom.N1();i++)
			{
				for(int j=0; j<g._dom.N2();j++)
				{
					for(int k=0; k<g._dom.N3();k++)
					{
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
#ifdef ENABLE_OMP
#pragma omp parallel for reduction (+:sum)
#endif
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

double Grid::integrate_Over_Dom()
{
	return sum()*_dom.dv();
}

Grid Grid::coulomb_Grid(double q, vector<double> R)
{
	Grid V(_dom);
	double x = 0;
	double y = 0;
	double z = 0;
	double v = 0;
#ifdef ENABLE_OMP
#pragma omp parallel for
#endif
	for(int i=0; i<_dom.N1();i++)
	{	
		for(int j=0; j<_dom.N2();j++)
		{		
			for(int k=0; k<_dom.N3();k++)
			{	
				for(int l=0; l<_dom.Nval();l++)
				{
					x = _dom.O()[ 0 ] + i*_dom.T()[0][0] + j*_dom.T()[0][1] +  k*_dom.T()[0][2]; 
					y = _dom.O()[ 1 ] + i*_dom.T()[1][0] + j*_dom.T()[1][1] +  k*_dom.T()[1][2]; 
					z = _dom.O()[ 2 ] + i*_dom.T()[2][0] + j*_dom.T()[2][1] +  k*_dom.T()[2][2];
					v = q/sqrt( (x-R[0])*(x-R[0])+(y-R[1])*(y-R[1])+(z-R[2])*(z-R[2]));
					if(v<1e-10)
					{
						v=0;
					}
					V._V[i][j][k][l]= v;
				}
			}
		}
	}
	return V;	
}

void Grid::coefs_Laplacian(int nBound, vector<double>& fcx, vector<double>& fcy, vector<double>& fcz, double& cc)
{
	vector<double> coefs(nBound+1);	
		if(nBound==1)
		{
			vector<double> c = {-2.0, 1.0};
			for(int i=0;i<=nBound;i++)
					coefs[i] = c[i];
		}
		else if(nBound==2)
		{
			double denom = 12.0;
			vector<double> c = {-30.0, 16.0, -1.0};
			for(int i=0;i<=2;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==3)
		{
			double denom = 180.0;
			vector<double> c = {-490.0, 270.0,-27.0, 2.0};
			for(int i=0;i<=3;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==4)
		{
			double denom = 5040.0;
			vector<double> c = {-14350.0, 8064.0, -1008.0, 128.0, -9.0};
			for(int i=0;i<=4;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==5)
		{
			double denom = 25200.0;
			vector<double> c = {-73766.0, 42000.0, -6000.0, 1000.0, -125.0, 8.0};
			for(int i=0;i<=5;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==6)
		{
			double denom = 831600.0;
			vector<double> c = {-2480478.0,1425600.0,-222750.0,44000.0,-7425.0,864.0,-50.0};
			for(int i=0;i<=6;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==7)
		{
			double denom = 75675600.0;
			vector<double> c = {-228812298.0,132432300.0,-22072050.0,4904900.0,-1003275.0, 160524.0,-17150.0,900.0};
			for(int i=0;i<=7;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==8)
		{
			double denom = 302702400.0;
			vector<double> c = {-924708642.0,538137600.0,-94174080.0,22830080.0,-5350800.0,1053696.0,-156800.0,15360.0,-735.0};
			for(int i=0;i<=8;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==9)
		{
			double denom = 15437822400.0;
			vector<double> c = {-47541321542.0,+27788080320.0, -5052378240.0,+1309875840.0,-340063920.0,+77728896.0,-14394240.0,+1982880.0,-178605.0,+7840.0};
			for(int i=0;i<=9;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==10)
		{
			double denom = 293318625600.0;
			vector<double> c = {-909151481810.0,+533306592000.0, -99994986000.0,+27349056000.0,-7691922000.0,+1969132032.0,-427329000.0,+73872000.0,-9426375.0,+784000.0,-31752.0};
			for(int i=0;i<=10;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==11)
		{
			double denom = 3226504881600.0;
			vector<double> c = {-10053996959110.0,+5915258949600.0,-1137549798000.0,+325014228000.0,-97504268400.0,+27301195152.0,-6691469400.0,+1365606000.0,-220114125.0,+26087600.0,-2012472.0,+75600.0};
			for(int i=0;i<=11;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==12)
		{
			double denom = 74209612276800.0;
			vector<double> c = {-232272619118930.0,+137002361126400.0,-26911178078400.0,+7973682393600.0,-2522922944850.0,+759845028096.0,-205205061600.0,+47609337600.0,-9112724775.0,+1371462400.0,-151484256.0,+10886400.0,-381150.0};
			for(int i=0;i<=12;i++)
				coefs[i] = c[i]/denom;
		}
		else
		{
			exit(1);
		}
		cc = 1/(_dom.dx()*_dom.dx()) + 1/(_dom.dy()*_dom.dy()) + 1/(_dom.dz()*_dom.dz());
		cc *= coefs[0];

		for(int i=0;i<=nBound;i++)
		{
			fcx[i] =  coefs[i]/(_dom.dx()*_dom.dx());
			fcy[i] =  coefs[i]/(_dom.dy()*_dom.dy());
			fcz[i] =  coefs[i]/(_dom.dz()*_dom.dz());
		}
}

Grid Grid::laplacian(int nBound)
{
	try
	{
		if(nBound<1 or nBound>12)
		{
			throw string("nBound oustide bounds");
		}
		else
		{
			vector<double> fcx(nBound);
			vector<double> fcy(nBound); 
			vector<double> fcz(nBound);
			double cc=0;
			Grid g(_dom);
			cout<<"end create grid"<<endl;
			coefs_Laplacian(nBound, fcx, fcy, fcz, cc);
			cout<<"coefs done"<<endl;
#ifdef ENABLE_OMP
#pragma omp parallel for
#endif
			for(int i=nBound;i<g._dom.N1()-nBound;i++)
			{
				for(int j=nBound;j<g._dom.N2()-nBound;j++)
				{	
					for(int k=nBound;k<g._dom.N3()-nBound;k++)
					{
						for(int l=0; l<g._dom.Nval();l++)
						{
							double v = cc*_V[i][j][k][l];
							for(int n=1;n<=nBound;n++)
							{
								v += fcx[n] *(_V[i-n][j][k][l]+_V[i+n][j][k][l]);
								v += fcy[n] *(_V[i][j-n][k][l]+_V[i][j+n][k][l]);
								v += fcz[n] *(_V[i][j][k-n][l]+_V[i][j][k+n][l]);
							}
							g._V[i][j][k][l]=v;
						}
					}
				}
			}
			return g;
		}
	}
	catch(string error)
	{
		cout<<error<<endl;
		exit(1);
	}				
}

/*
void Grid::coefs_Gradient(int nBound, vector<double>& fcx, vector<double>& fcy, vector<double>& fcz, double& cc)
{
	vector<double> coefs(nBound+1);	
		if(nBound==1)
		{
			vector<double> c = {-2.0, 1.0};
			for(int i=0;i<=nBound;i++)
					coefs[i] = c[i];
		}
		else if(nBound==2)
		{
			double denom = 12.0;
			vector<double> c = {-30.0, 16.0, -1.0};
			for(int i=0;i<=2;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==3)
		{
			double denom = 180.0;
			vector<double> c = {-490.0, 270.0,-27.0, 2.0};
			for(int i=0;i<=3;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==4)
		{
			double denom = 5040.0;
			vector<double> c = {-14350.0, 8064.0, -1008.0, 128.0, -9.0};
			for(int i=0;i<=4;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==5)
		{
			double denom = 25200.0;
			vector<double> c = {-73766.0, 42000.0, -6000.0, 1000.0, -125.0, 8.0};
			for(int i=0;i<=5;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==6)
		{
			double denom = 831600.0;
			vector<double> c = {-2480478.0,1425600.0,-222750.0,44000.0,-7425.0,864.0,-50.0};
			for(int i=0;i<=6;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==7)
		{
			double denom = 75675600.0;
			vector<double> c = {-228812298.0,132432300.0,-22072050.0,4904900.0,-1003275.0, 160524.0,-17150.0,900.0};
			for(int i=0;i<=7;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==8)
		{
			double denom = 302702400.0;
			vector<double> c = {-924708642.0,538137600.0,-94174080.0,22830080.0,-5350800.0,1053696.0,-156800.0,15360.0,-735.0};
			for(int i=0;i<=8;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==9)
		{
			double denom = 15437822400.0;
			vector<double> c = {-47541321542.0,+27788080320.0, -5052378240.0,+1309875840.0,-340063920.0,+77728896.0,-14394240.0,+1982880.0,-178605.0,+7840.0};
			for(int i=0;i<=9;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==10)
		{
			double denom = 293318625600.0;
			vector<double> c = {-909151481810.0,+533306592000.0, -99994986000.0,+27349056000.0,-7691922000.0,+1969132032.0,-427329000.0,+73872000.0,-9426375.0,+784000.0,-31752.0};
			for(int i=0;i<=10;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==11)
		{
			double denom = 3226504881600.0;
			vector<double> c = {-10053996959110.0,+5915258949600.0,-1137549798000.0,+325014228000.0,-97504268400.0,+27301195152.0,-6691469400.0,+1365606000.0,-220114125.0,+26087600.0,-2012472.0,+75600.0};
			for(int i=0;i<=11;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==12)
		{
			double denom = 74209612276800.0;
			vector<double> c = {-232272619118930.0,+137002361126400.0,-26911178078400.0,+7973682393600.0,-2522922944850.0,+759845028096.0,-205205061600.0,+47609337600.0,-9112724775.0,+1371462400.0,-151484256.0,+10886400.0,-381150.0};
			for(int i=0;i<=12;i++)
				coefs[i] = c[i]/denom;
		}
		else
		{
			exit(1);
		}
		cc = 1/(_dom.dx()*_dom.dx()) + 1/(_dom.dy()*_dom.dy()) + 1/(_dom.dz()*_dom.dz());
		cc *= coefs[0];

		for(int i=0;i<=nBound;i++)
		{
			fcx[i] =  coefs[i]/_dom.dx();
			fcy[i] =  coefs[i]/_dom.dy()*_dom.dy();
			fcz[i] =  coefs[i]/_dom.dz()*_dom.dz();
		}
}


*/














