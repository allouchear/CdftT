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
	reset();
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
#pragma omp parallel for shared(sum,g,_V)
#endif
#ifdef ENABLE_ACC
#pragma acc kernels loop
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
#ifdef ENABLE_ACC
#pragma acc kernels loop
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
#ifdef ENABLE_ACC
#pragma acc kernels loop
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
#ifdef ENABLE_ACC
#pragma acc kernels loop
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
#ifdef ENABLE_ACC
#pragma acc kernels loop reduction(+:sum)
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
#ifdef ENABLE_ACC
#pragma acc kernels loop
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
			g.set_str(_str);
			cout<<"end create grid"<<endl;
			coefs_Laplacian(nBound, fcx, fcy, fcz, cc);
			cout<<"coefs done"<<endl;
#ifdef ENABLE_OMP
//cout<<"Number of cores = "<<omp_get_num_procs()<<endl;
#pragma omp parallel for shared(g,_V,fcx,fcy,fcz,cc,nBound)
#endif
#ifdef ENABLE_ACC
#pragma acc kernels loop
#endif
			//for(int i=nBound;i<g._dom.N1()-nBound;i++) // DEBUG
			for(int ii=0;ii<g._dom.N1()-2*nBound;ii++)
			{
				int i=ii+nBound;
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
		if(nBound<1 or nBound>8)
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
			g._str=_str;
			coefs_Gradient(nBound, fcx, fcy, fcz, cc);
#ifdef ENABLE_OMP
#pragma omp parallel for
#endif
#ifdef ENABLE_ACC
#pragma acc kernels loop
#endif
			for(int i=nBound;i<g._dom.N1()-nBound;i++) // DEBUG
			//for(int ii=0;ii<g._dom.N1()-2*nBound;ii++)
			{
				//int i=ii+nBound;
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
#ifdef ENABLE_ACC
#pragma acc kernels loop
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
#ifdef ENABLE_ACC
#pragma acc kernels loop
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
#ifdef ENABLE_ACC
#pragma acc kernels loop
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
#ifdef ENABLE_ACC
#pragma acc kernels loop
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
#ifdef ENABLE_ACC
#pragma acc kernels loop
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
#ifdef ENABLE_ACC
#pragma acc kernels loop
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
#ifdef ENABLE_ACC
#pragma acc kernels loop
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
#ifdef ENABLE_ACC
#pragma acc kernels loop
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
	g._str=_str;
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
#ifdef ENABLE_ACC
#pragma acc kernels loop
#endif
	for(int i = iXBegin ; i <= iXEnd-1 ; i++)
	{
		int x0, xp, xm;
		x0 = 2 * i;
		xp = x0 + 1;
		xm = x0 - 1;
		for(int j = iYBegin ; j <= iYEnd-1 ; j++)
		{
			int y0, yp, ym;
			y0 = 2 * j;
			yp = y0 + 1;
			ym = y0 - 1;
			for(int k= iZBegin ; k <= iZEnd-1 ; k++)
			{
				int z0, zp, zm;
				z0 = 2 * k;
				zp = z0 + 1;
				zm = z0 - 1;
				for(int l=0;l<_dom.Nval(); l++)
				{ 	
					double face, corner, edge;
					
					face = _V [xm][y0][z0][l] +_V [xp][y0][z0][l] +_V [x0][ym][z0][l] +_V [x0][yp][z0][l] + _V [x0][y0][zm][l] +_V [x0][y0][zp][l];

					corner =  _V [xm] [ym] [zm][l] +_V [xm] [ym] [zp][l] +_V [xm] [yp] [zm][l] +_V [xm] [yp] [zp][l]+_V [xp] [ym] [zm][l] +_V [xp] [ym] [zp][l] +_V [xp] [yp] [zm][l] +_V [xp] [yp] [zp][l];
					
					edge = _V [xm] [y0] [zm][l] +_V [xm] [ym] [z0][l] +_V [xm] [yp] [z0][l] +_V [xm] [y0] [zp][l] +_V [x0] [ym] [zm][l] +_V [x0] [yp] [zm][l] +_V [x0] [ym] [zp][l] +_V [x0] [yp] [zp][l] +_V [xp] [y0] [zm][l] +_V [xp] [ym] [z0][l] +_V [xp] [yp] [z0][l] +_V [xp] [y0] [zp][l];
	
					g._V [i][j][k][l]=scale * (8.0 * _V [x0][y0][z0][l] +4.0 * face +2.0 * edge +corner);
          			}
			}
		}
	}
	return g;
}

void Grid::save(ofstream& nameFile)
{
	nameFile.precision(14);
	nameFile<<scientific;
	nameFile<<_str.number_of_atoms()<<" ";
	for(int i=0;i<3;i++)
	{
		nameFile<<_dom.O()[i]<<" ";
	}
	nameFile<<_dom.Nval()<<" ";
	nameFile<<endl<<_dom.N1()<<" ";
	if(_dom.N1()<0)
	{
		for(int i=0;i<3;i++)
		{
			nameFile<<_dom.Tij(0,i)<<" ";
		}
		nameFile<<endl;
		nameFile<<_dom.N2()<<" ";
		for(int i=0;i<3;i++)
		{
			nameFile<<_dom.Tij(1,i)<<" ";
		}
		nameFile<<endl<<_dom.N3()<<" ";
		for(int i=0;i<3;i++)
		{
			nameFile<<_dom.Tij(2,i)<<" ";
		}
	}
	else
	{
		for(int i=0;i<3;i++)
		{
			nameFile<<_dom.Tij(0,i)*BOHRTOANG<<" ";
		}
		nameFile<<endl;
		nameFile<<_dom.N2()<<" ";
		for(int i=0;i<3;i++)
		{
			nameFile<<_dom.Tij(1,i)*BOHRTOANG<<" ";
		}
		nameFile<<endl<<_dom.N3()<<" ";
		for(int i=0;i<3;i++)
		{
			nameFile<<_dom.Tij(2,i)*BOHRTOANG<<" ";
		}
	}
	nameFile<<endl;
	for(int i=1;i<=_str.number_of_atoms();i++)
	{
		nameFile<<_str.atom(i).atomic_number()<<" ";
		nameFile<<_str.atom(i).charge()<<" ";
		for(int j=1; j<=3;j++)
		{
			nameFile<<_str.atom(i).coordinates()[j]<<" ";
		}
		nameFile<<endl;
	}
	for(int i=0; i<_dom.N1();i++)
	{	
		
		for(int j=0; j<_dom.N2();j++)
		{	
			int R=0;
			for(int k=0; k<_dom.N3();k++)
			{
				for(int l=0; l<_dom.Nval();l++)
				{
					nameFile<<_V[i][j][k][l]<<" ";
					R++;
					if(R%6==0)
					{
						nameFile<<endl;
					}
				}
			}
			if(R%6!=0)
			{
				nameFile<<endl;
			}
		}
	}
}
vector<double> Grid::find_max_neighbour(int i, int j, int k,int xm,int xp, int ym,int yp, int zm, int zp, double Current, vector<vector<double>>& traj, const Grid& g, int& repeat,vector<double>& v)
{	
	repeat++;
	cout<<"xm "<<xm<<endl;
	cout<<"xp "<<xp<<endl;
	cout<<"ym "<<ym<<endl;
	cout<<"yp "<<yp<<endl;
	cout<<"zm "<<zm<<endl;
	cout<<"zp "<<zp<<endl;
	cout<<Current<<endl;
	if(-_V [xm][j][k][1]>Current)
	{
		cout<<"yo are here"<<endl;
		Current= -_V [xm][j][k][1];
		v[0]=xm;
		v[1]=j;
		v[2]=k;
		cout<<"1"<<endl;
	}
	if(_V [xp][j][k][1]>Current)
	{
		Current=_V [xp][j][k][1];
		v[0]=xp;
		v[1]=j;
		v[2]=k;
		cout<<"2"<<endl;
	}
	if(-_V [i][ym][k][2]>Current)
	{
		Current= -_V [i][ym][k][2];
		v[0]=i;
		v[1]=ym;
		v[2]=k;
		cout<<"3"<<endl;
	}
	if(_V [i][yp][k][2] >Current)
	{
		Current=_V [i][yp][k][2];
		v[0]=i;
		v[1]=yp;
		v[2]=k;
		cout<<"4"<<endl;
	}
	if(-_V [i][j][zm][3]>Current)
	{
		Current= -_V [i][j][zm][3];
		v[0]=i;
		v[1]=j;
		v[2]=zm;
		cout<<"5"<<endl;
	}
	if(_V [i][j][zp][3]>Current)
	{
		Current=_V [i][j][zp][3];
		v[0]=i;
		v[1]=j;
		v[2]=zp;
		cout<<"6"<<endl;
	}
	if((-_V [xm] [ym] [zm][1]-_V [xm] [ym] [zm][2]-_V [xm] [ym] [zm][3])>Current)
	{
		Current= -_V [xm] [ym] [zm][1]-_V [xm] [ym] [zm][2]-_V [xm] [ym] [zm][3];
		v[0]=xm;
		v[1]=ym;
		v[2]=zm;
		cout<<"7"<<endl;
	}
	if((-_V [xm] [ym] [zp][1]-_V [xm] [ym] [zp][2]+_V [xm] [ym] [zp][3])>Current)
	{
		Current= -_V [xm] [ym] [zp][1]-_V [xm] [ym] [zp][2]+_V [xm] [ym] [zp][3];
		v[0]=xm;
		v[1]=ym;
		v[2]=zp;
		cout<<"8"<<endl;
	}
	if((-_V [xm] [yp] [zm][1]+_V [xm] [yp] [zm][2]-_V [xm] [yp] [zm][3])>Current)
	{
		Current= -_V [xm] [yp] [zm][1]+_V [xm] [yp] [zm][2]-_V [xm] [yp] [zm][3];
		v[0]=xm;
		v[1]=yp;
		v[2]=zm;
		cout<<"9"<<endl;
	}
	if((-_V [xm] [yp] [zp][1]+_V [xm] [yp] [zp][2]+_V [xm] [yp] [zp][3])>Current)
	{
		Current= -_V [xm] [yp] [zp][1]+_V [xm] [yp] [zp][2]+_V [xm] [yp] [zp][3];
		v[0]=xm;
		v[1]=yp;
		v[2]=zp;
		cout<<"10"<<endl;
	}
	if((_V [xp] [ym] [zm][1]-_V [xp] [ym] [zm][2]-_V [xp] [ym] [zm][3])>Current)
	{
		Current=_V [xp] [ym] [zm][1]-_V [xp] [ym] [zm][2]-_V [xp] [ym] [zm][3];
		v[0]=xp;
		v[1]=ym;
		v[2]=zm;
		cout<<"11"<<endl;
	}
	if((_V [xp] [ym] [zp][1]-_V [xp] [ym] [zp][2]+_V [xp] [ym] [zp][3])>Current)
	{
		Current=_V [xp] [ym] [zp][1]-_V [xp] [ym] [zp][2]+_V [xp] [ym] [zp][3];
		v[0]=xp;
		v[1]=ym;
		v[2]=zp;
		cout<<"12"<<endl;
	}
	if((_V [xp] [yp] [zm][1]+_V [xp] [yp] [zm][2]-_V [xp] [yp] [zm][3])>Current)
	{
		Current=_V [xp] [yp] [zm][1]+_V [xp] [yp] [zm][2]-_V [xp] [yp] [zm][3];
		v[0]=xp;
		v[1]=yp;
		v[2]=zm;
		cout<<"13"<<endl;
	}
	if((_V [xp] [yp] [zp][1]+_V [xp] [yp] [zp][2]+_V [xp] [yp] [zp][3])>Current)
	{
		Current=_V [xp] [yp] [zp][1]+_V [xp] [yp] [zp][2]+_V [xp] [yp] [zp][3];
		v[0]=xp;
		v[1]=yp;
		v[2]=zp;
		cout<<"14"<<endl;
	}
	if((-_V [xm] [j] [zm][1]-_V [xm] [j] [zm][3])>Current)
	{
		Current= -_V [xm] [j] [zm][1]-_V [xm] [j] [zm][3];
		v[0]=xm;
		v[1]=j;
		v[2]=zm;
		cout<<"15"<<endl;
	}
	if((-_V [xm] [ym] [k][1]-_V [xm] [ym] [k][2])>Current)
	{
		Current= -_V [xm] [ym] [k][1]-_V [xm] [ym] [k][2];
		v[0]=xm;
		v[1]=ym;
		v[2]=k;
		cout<<"16"<<endl;
	}
	if((-_V [xm] [yp] [k][1]+_V [xm] [yp] [k][2])>Current)
	{
		Current= -_V [xm] [yp] [k][1]+_V [xm] [yp] [k][2];
		v[0]=xm;
		v[1]=yp;
		v[2]=k;
		cout<<"17"<<endl;
	}
	if((-_V [xm] [j] [zp][1]+_V [xm] [j] [zp][3])>Current)
	{
		Current= -_V [xm] [j] [zp][1]+_V [xm] [j] [zp][3];
		v[0]=xm;
		v[1]=j;
		v[2]=zp;
		cout<<"18"<<endl;
	}
	if((-_V [i] [ym] [zm][2]-_V [i] [ym] [zm][3])>Current)
	{
		Current= -_V [i] [ym] [zm][2]-_V [i] [ym] [zm][3];
		v[0]=i;
		v[1]=ym;
		v[2]=zm;
		cout<<"19"<<endl;
	}
	if((_V [i] [yp] [zm][2]-_V [i] [yp] [zm][3])>Current)
	{
		Current=_V [i] [yp] [zm][2]-_V [i] [yp] [zm][3];
		v[0]=i;
		v[1]=yp;
		v[2]=zm;
		cout<<"20"<<endl;
	}
	if((-_V [i] [ym] [zp][2]+_V [i] [ym] [zp][3])>Current)
	{
		Current=-_V [i] [ym] [zp][2]+_V [i] [ym] [zp][3];
		v[0]=i;
		v[1]=ym;
		v[2]=zp;
		cout<<"21"<<endl;
	}
	if((_V [i] [yp] [zp][2]+_V [i] [yp] [zp][3])>Current)
	{
		Current=_V [i] [yp] [zp][2]+_V [i] [yp] [zp][3];
		v[0]=i;
		v[1]=yp;
		v[2]=zp;
		cout<<"22"<<endl;
	}
	if((_V [xp] [j] [zm][1]-_V [xp] [j] [zm][3])>Current)
	{
		Current=_V [xp] [j] [zm][1]-_V [xp] [j] [zm][3];
		v[0]=xp;
		v[1]=j;
		v[2]=zm;
		cout<<"23"<<endl;
	}
	if((_V [xp] [ym] [k][1]-_V [xp] [ym] [k][2])>Current)
	{
		Current=_V [xp] [ym] [k][1]-_V [xp] [ym] [k][2];
		v[0]=xp;
		v[1]=ym;
		v[2]=k;
		cout<<"24"<<endl;
	}
	if((_V [xp] [yp] [k][1]+_V [xp] [yp] [k][2])>Current)
	{
		Current=_V[xp] [yp] [k][1]+_V [xp] [yp] [k][2];
		v[0]=xp;
		v[1]=yp;
		v[2]=k;
		cout<<"25"<<endl;
	}
	
	if((_V [xp] [j] [zp][1]+_V [xp] [j] [zp][3])>Current)
	{
		Current=_V [xp] [j] [zp][1]+_V [xp] [j] [zp][3];
		v[0]=xp;
		v[1]=j;
		v[2]=zp;
		cout<<"26"<<endl;
	}
	
	//test if point is a max
	bool HiMax=true;
	if(v[0]!=i)
	{
		HiMax=false;
	}
	if(v[1]!=j)
	{
		HiMax=false;
	}
	if(v[2]!=k)
	{
		HiMax=false;
	}
	
	v[3]=Current;
	
	//if point is max: exit function
	if(HiMax)
	{
		cout<<"end recion"<<endl;
		return v;
	}
	else
	{
		//test if new point is indexed
		if(g._V[v[0]][v[1]][v[2]][0]!=0)
		{
			return v;
		}
		//repeat function for new point
		else
		{
			cout<<"are we doing it?"<<endl;
			traj.push_back(v);
			cout<<"weve done it"<<endl;
			for(int i=0;i<int(v.size());i++)
			{
				cout<<v[i]<<endl;
			}
			find_max_neighbour(v[0],v[1],v[2],v[0]-1,v[0]+1,v[1]-1,v[1]+1,v[2]-1,v[2]+1,v[3],traj,g,repeat,v);
		}
	}
	
	//have to index points
	//if its been indexed have to break loop
}
/*
double Grid::atom_attract_diff(int j,const vector<double>& attract)
{
	double v;
	for(int i=0;i<3;i++)
	{
		v+= (_str.atom()[j+1].coords()[i]-attract[i])*(_str.atom()[j+1].coords()[i]-attract[i]);
	}
	v=sqrt(v);
	return v;
}
*/
Grid Grid::aim_On_Grid(int nBound)
{
	try
	{
		if(_dom.Nval()==4)
		{
			cout<<"begin aim"<<endl;
			Grid g(Domain(1,_dom.N1(), _dom.N2(), _dom.N3(), _dom.O()));
			for(int i=0;i<3;i++)
			{
				for(int j=0;j<3;j++)
				{
					g._dom.set_T(_dom.Tij(i,j)/2, i, j);
				}
			}
			g._str=_str;
			vector<vector<double>> attractors;
			/*vector<int> attractIndex(_str.number_of_atoms);
			for(int i=0; i<_str.number_of_atoms;i++)
			{
				attractIndex[i]=i+1;
			}*/
			
			for(int i=nBound;i<_dom.N1()-nBound;i++)
			{
				cout<<"i "<<i<<endl;
				int xp, xm;
				xp = i + 1;
				xm = i - 1;
				for(int j=nBound;j<_dom.N2()-nBound;j++)
				{
					cout<<"j "<<j<<endl;
					int yp,ym;
					yp = j + 1;
					ym = j - 1;
					for(int k=nBound;k<_dom.N3()-nBound;k++)
					{
						cout<<"k "<<k<<endl;
						int zp,zm;
						zp = k + 1;
						zm = k - 1;
							
						if(g._V[i][j][k][0]!=0)
						{
							break;
						}
						else
						{	
							double w=_V[i][j][k][1]*_V[i][j][k][1]+_V[i][j][k][2]*_V[i][j][k][2]+_V[i][j][k][3]*_V[i][j][k][3];
							w=sqrt(w);
							vector<vector<double>> trajectory={{double(i),double(j),double(k),w}};
							cout<<"Begin recursion"<<endl;
							int repeat=0;
							vector<double> v=trajectory[0];
							while(uplus-u>1e-10)
							{
								u=uplus;
								uplus = find_max_neighbour(i,j,k,xm,xp,ym,yp,zm,zp,w,trajectory,g,repeat,v)[3];
							}
							cout<<"end recu"<<endl;
							//if point has been already indexed exit loop: new point
							if(g._V[int(u[0])][int(u[1])][int(u[2])][0]!=0)
							{
								break;
							}
							
							//Do we know max?
							bool WhoIsMax=true;
							//compare max to list
							cout<<"done"<<endl;
							for(int m=0;m<int(attractors.size());m++)
							{
								if(attractors[m]==u)
								{
									int I=0;
									WhoIsMax=false;
									I++;
									for(int p=0;p<int(trajectory.size());p++)
									{
										g._V[int(trajectory[p][0])][int(trajectory[p][1])][int(trajectory[p][2])][0]=double(I);
									}
								}
							}
							cout<<"done"<<endl;
							//dont know max
							if(WhoIsMax)
							{
								//add to known attractors
								attractors.push_back(u);
								//index all points from trajectory
								for(int p=0;p<int(trajectory.size());p++)
								{
									g._V[int(trajectory[p][0])][int(trajectory[p][1])][int(trajectory[p][2])][0]=double(trajectory.size());
								}
							}
						}					
					}
				}
			}
			/*for(int m=0; m<int(attractors.size());m++)
			{
				for(int a=0;a<4;a++)
				{
					cout<<attractors[m][a]<<endl;
				}
			}*/
			return g;
		}
		else
		{
			throw string("Gradient grid must be used for AIM");
		}
	}
	
	catch(string error)
	{
		cout<<error<<endl;
		exit(1);
	}
}







