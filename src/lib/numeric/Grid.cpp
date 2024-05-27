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
//cout<<"Number of cores = "<<omp_get_num_procs()<<endl;
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


void Grid::coefs_Gradient(int nBound, vector<double>& fcx, vector<double>& fcy, vector<double>& fcz, double& cc)
{
	vector<double> coefs(nBound+1);	
		if(nBound==1)
		{
			double denom = 2.0;
			vector<double> c = {-1.0};
			for(int i=0;i<nBound;i++)
					coefs[i] = c[i]/denom;
		}
		else if(nBound==2)
		{
			double denom = 12.0;
			vector<double> c = { 1.0, -8.0};
			for(int i=0;i<2;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==3)
		{
			double denom = 60.0;
			vector<double> c ={ -1.0, +9.0, -45.0};
			for(int i=0;i<3;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==4)
		{
			double denom = 840.0;
			vector<double> c = { 3.0, -32.0, +168.0, -672.0};
			for(int i=0;i<4;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==5)
		{
			double denom = 2520.0;
			vector<double> c = { -2.0, +25.0, -150.0,+600.0, -2100.0};
			for(int i=0;i<5;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==6)
		{
			double denom = 27720.0;
			vector<double> c = { 5.0, -72.0, +495.0, -2200.0, +7425.0, -23760.0};
			for(int i=0;i<6;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==7)
		{
			double denom = 360360.0;
			vector<double> c = { -15.0, +245.0, -1911.0, +9555.0, -35035.0, +105105.0, -315315.0};
			for(int i=0;i<7;i++)
				coefs[i] = c[i]/denom;
		}
		else if (nBound==8)
		{
			double denom = 720720.0;
			vector<double> c = { 7.0, -128.0, +1120.0, -6272.0, +25480.0, -81536.0, +224224.0, -640640.0};
			for(int i=0;i<8;i++)
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


Grid Grid::gradient(int nBound)
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
			_dom.set_Nval(4);
			Grid g(_dom);
			cout<<"end create grid"<<endl;
			coefs_Gradient(nBound, fcx, fcy, fcz, cc);
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
								if(l==0)
								{
									v= v/cc;
									break;
									
								}
								if(l==1)
								{
								v += fcx[n] *(_V[i-n][j][k][l]+_V[i+n][j][k][l]);
								}
								if(l==2)
								{
								v += fcy[n] *(_V[i][j-n][k][l]+_V[i][j+n][k][l]);
								}
								if(l==3)
								{
								v += fcz[n] *(_V[i][j][k-n][l]+_V[i][j][k+n][l]);
								}
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
void bicub_Coef(vector<double> z, vector<double> dzdx, vector<double> dzdy, vector<double> d2zdxdy, double dx, double dy, vector<vector<double>> c)
{
	vector<vector<int>> Ainv[16][16]={
		{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},{-3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0}{2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0},{0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},{0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1},{0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1},{-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0},{9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2},{-6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2},{2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0},{-6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1},{4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1}};
	int i,j,k,l;
	double s, dxdy;
	vector<double> c1D(16),x(16);
	dxdy= dx*dy;
	for (i=0;i<4;i++)
	{
		x[i]=z[i];
		x[i+4]=dzdx[i]*dx;
		x[i+8]=dzdy[i]*dy;
		x[i+12]=d2zdxdy[i]*dxdy;
	}
	for (i=0;i<16;i++)
	{
		s=0.0;
		for (k=0;k<16;k++) 
		{
			s += Ainv[i][k]*x[k];
		}
		c1D[i]=s;
	}
	l=0;
	for (i=0;i<4;i++)
	{
		for (j=0;j<4;j++)
		{
		c[i][j]=c1D[l++];
		}
	}
}
*/



Grid Grid::finer_Grid()
{
	Grid g(Domain(_dom.Nval() ,2*_dom.N1()-1,2*_dom.N2()-1,2*_dom.N3()-1, _dom.O()));
	g._str=_str;
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			g._dom.set_T(_dom.Tij(i,j)/2, i, j);
		}
	}
	//coarse grid points to fine grid

#ifdef ENABLE_OMP
#pragma omp parallel for 
#endif
	for(int i=0;i<_dom.N1();i++)
	{	
		for(int j=0;j<_dom.N2();j++)
		{	
			for(int k=0;k<_dom.N3();k++)
			{
				g._V[2*i][2*j][2*k]=_V[i][j][k];
			}	
		}
	}
		
	//center cube points
#ifdef ENABLE_OMP
#pragma omp parallel for 
#endif
	for(int i=1;i<g._dom.N1()-1;i+=2)
	{	
		for(int j=1;j<g._dom.N2()-1;j+=2)
		{	
			for(int k=1;k<g._dom.N3()-1;k+=2)
			{
				for(int l=0; l<g._dom.Nval();l++)
				{	
					g._V[i][j][k][l] = 0.125*g._V[i-1][j-1][k-1][l] + 0.125*g._V[i-1][j-1][k+1][l] + 0.125*g._V[i-1][j+1][k-1][l] + 0.125*g._V[i-1][j+1][k+1][l] + 0.125*g._V[i+1][j-1][k-1][l] + 0.125*g._V[i+1][j-1][k+1][l] + 0.125*g._V[i+1][j+1][k-1][l] + 0.125*g._V[i+1][j+1][k+1][l];
				}	
			}
		}
	}
	
				
	//cubes face points
#ifdef ENABLE_OMP
#pragma omp parallel for 
#endif
	for(int i=0;i<g._dom.N1()-1;i+=2)
	{	
		for(int j=0;j<g._dom.N2()-1;j+=2)
		{	
			for(int k=1;k<g._dom.N3()-1;k+=2)
			{
				for(int l=0; l<g._dom.Nval();l++)
				{
					g._V[i][j][k][l]= 0.5*g._V[i][j][k-1][l]+0.5*g._V[i][j][k+1][l];
				}
			}	
		}
	}

#ifdef ENABLE_OMP
#pragma omp parallel for 
#endif
	for(int i=0;i<g._dom.N1()-1;i+=2)
	{	
		for(int j=1;j<_dom.N2()-1;j+=2)
		{	
			for(int k=0;k<_dom.N3()-1;k+=2)
			{
				for(int l=0; l<g._dom.Nval();l++)
				{
					g._V[i][j][k][l] = 0.5*g._V[i][j-1][k][l]+0.5*g._V[i][j+1][k][l];
				}
			}	
		}
	}

#ifdef ENABLE_OMP
#pragma omp parallel for 
#endif
	for(int i=1;i<g._dom.N1()-1;i+=2)
	{	
		for(int j=0;j<g._dom.N2()-1;j+=2)
		{	
			for(int k=0;k<g._dom.N3()-1;k+=2)
			{
				for(int l=0; l<g._dom.Nval();l++)
				{
					g._V[i][j][k][l] = 0.5*g._V[i-1][j][k][l]+0.5*g._V[i+1][j][k][l];
				}
			}	
		}
	}

#ifdef ENABLE_OMP
#pragma omp parallel for 
#endif
	for(int i=0;i<g._dom.N1()-1;i+=2)
	{	
		for(int j=1;j<g._dom.N2()-1;j+=2)
		{	
			for(int k=1;k<g._dom.N3()-1;k+=2)
			{
				for(int l=0; l<g._dom.Nval()-1;l++)
				{
					g._V[i][j][k][l] = 0.25*g._V[i][j-1][k-1][l] + 0.25*g._V[i+1][j-1][k+1][l] + 0.25*g._V[i][j+1][k-1][l] + 0.25*g._V[i][j+1][k+1][l];
				}
			}	
		}
	}
	
#ifdef ENABLE_OMP
#pragma omp parallel for 
#endif
	for(int i=1;i<g._dom.N1()-1;i+=2)
	{	
		for(int j=0;j<g._dom.N2()-1;j+=2)
		{	
			for(int k=1;k<g._dom.N3()-1;k+=2)
			{
				for(int l=0; l<g._dom.Nval();l++)
				{
					g._V[i][j][k][l] = 0.25*g._V[i-1][j][k-1][l] + 0.25*g._V[i-1][j][k+1][l] + 0.25*g._V[i+1][j][k-1][l] + 0.25*g._V[i+1][j][k+1][l];
				}
			}	
		}
	}

#ifdef ENABLE_OMP
#pragma omp parallel for 
#endif
	for(int i=1;i<g._dom.N1()-1;i+=2)
	{	
		for(int j=1;j<g._dom.N2()-1;j+=2)
		{	
			for(int k=0;k<g._dom.N3()-1;k+=2)
			{
				for(int l=0; l<g._dom.Nval();l++)
				{
					g._V[i][j][k][l] = 0.25*g._V[i-1][j-1][k][l] + 0.25*g._V[i+1][j-1][k][l] + 0.25*g._V[i-1][j+1][k][l] + 0.25*g._V[i+1][j+1][k][l];
				}
			}	
		}
	}
	return g;
}

Grid Grid::coarser_Grid()
{
	int N[3]={_dom.N1(),_dom.N2(), _dom.N3()};
	for(int i=0; i<3; i++)
	{
		if(N[i]%2==1)
		{
			N[i]=N[i]-1;
		}
		N[i] =N[i]/2;
	}
	Grid g(Domain(_dom.Nval(), N[0], N[1], N[2], _dom.O()));
	double scale = 1.0 / 64.0;
	

	/*printf("Begin restriction\n");*/

	int iXBegin = 1;
	int iXEnd = g._dom.N1();
	int iYBegin = 1;
	int iYEnd = g._dom.N2();
	int iZBegin = 1;
	int iZEnd = g._dom.N3();

#ifdef ENABLE_OMP
#pragma omp parallel for
#endif
	for(int i = iXBegin ; i <= iXEnd-1 ; i++)
	{
		cout<<"i "<<i<<endl;
		int x0, xp, xm;
		x0 = 2 * i;
		xp = x0 + 1;
		xm = x0 - 1;
		for(int j = iYBegin ; j <= iYEnd-1 ; j++)
		{
			cout<<"j "<<j<<endl;
			int y0, yp, ym;
			y0 = 2 * j;
			yp = y0 + 1;
			ym = y0 - 1;
			for(int k= iZBegin ; k <= iZEnd-1 ; k++)
			{
				cout<<"k "<<k<<endl;
				int z0, zp, zm;
				z0 = 2 * k;
				zp = z0 + 1;
				zm = z0 - 1;
				for(int l=0;l<_dom.Nval(); l++)
				{ 	
					cout<<"l "<<l<<endl;
					double face, corner, edge;
					
					face = _V [xm][y0][z0][l] +_V [xp][y0][z0][l] +_V [x0][ym][z0][l] +_V [x0][yp][z0][l] + _V [x0][y0][zm][l] +_V [x0][y0][zp][l];

					corner =  _V [xm] [ym] [zm][l] +_V [xm] [ym] [zp][l] +_V [xm] [yp] [zm][l] +_V [xm] [yp] [zp][l]+_V [xp] [ym] [zm][l] +_V [xp] [ym] [zp][l] +_V [xp] [yp] [zm][l] +_V [xp] [yp] [zp][l];
					
					edge = _V [xm] [y0] [zm][l] +_V [xm] [ym] [z0][l] +_V [xm] [yp] [z0][l] +_V [xm] [y0] [zp][l] +						_V [x0] [ym] [zm][l] +_V [x0] [yp] [zm][l] +_V [x0] [ym] [zp][l] +_V [x0] [yp] [zp][l] +_V [xp] [y0] [zm][l] +_V [xp] [ym] [z0][l] +_V [xp] [yp] [z0][l] +_V [xp] [y0] [zp][l];
	
					g._V [i][j][k][l]=scale * (8.0 * _V [x0][y0][z0][l] +4.0 * face +2.0 * edge +corner);
          			}
			}
		}
	}
	return g;
}

void Grid::save(ofstream& nameFile)
{
	nameFile<<_str.number_of_atoms();
	for(int i=0;i<3;i++)
	{
		nameFile<<_dom.O()[i];
	}
	nameFile<<_dom.Nval();
	nameFile<<endl<<_dom.N1();
	for(int i=0;i<3;i++)
	{
		nameFile<<_dom.Tij(1,i);
	}
	nameFile<<endl;
	nameFile<<_dom.N2();
	for(int i=0;i<3;i++)
	{
		nameFile<<_dom.Tij(2,i);
	}
	nameFile<<endl<<_dom.N3();
	for(int i=0;i<3;i++)
	{
		nameFile<<_dom.Tij(3,i);
	}
	nameFile<<endl;
	for(int i=0;i<_str.number_of_atoms();i++)
	{
		nameFile<<_str.atom(i).atomic_number();
		nameFile<<_str.atom(i).charge();
		for(int j=0; j<3;j++)
		{
			nameFile<<_str.atom(i).coordinates()[j];
		}
		nameFile<<endl;
	}
	for(int i=0; i<_dom.N1();i++)
	{	
		for(int j=0; j<_dom.N2();j++)
		{	
			for(int k=0; k<_dom.N3();k++)
			{
				for(int l=0; l<_dom.Nval();l++)
				{
					nameFile<<_V[i][j][k][l];
				}
			}
		}
	}
	
	
}








