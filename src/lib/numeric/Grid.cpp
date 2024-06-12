using namespace std;
#include <numeric/Grid.h>
#include <numeric/GridCP.h>
#include <common/Constants.h>
#include <cmath>
#include <list>
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
	cout<<"End read domain "<<endl;
	_dom=d;
	Structure S(nameFile, Natoms, Table);
	cout<<"End read struct "<<endl;
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
	cout<<"End read "<<endl;
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

void Grid::reset_Boundary(int nBound)
{
	int NV=(_dom.Nval()>1)?1:0;
	/* left */
	for(int i=0;i<nBound;i++)
	{
		for(int j=0;j<_dom.N2();j++)
		{
			for(int k=0;k<_dom.N3();k++)
			{	
				for(int l=NV; l<_dom.Nval();l++)
				{
				_V[i][j][k][l] = _V[nBound][j][k][l];
				}
			}
		}
	}
	
	/* right */
	for(int i=_dom.N1()-nBound;i<_dom.N1();i++)
	{	
		for(int j=0;j<_dom.N2();j++)
		{
			for(int k=0;k<_dom.N3();k++)
			{	
				for(int l=NV; l<_dom.Nval();l++)
				{
					_V[i][j][k][l] = _V[_dom.N1()-nBound-1][j][k][l];
				}
			}
		}
	}
	/* front */
	for(int j=0;j<nBound;j++)
	{
		for(int i=0;i<_dom.N1();i++)
		{	
			for(int k=0;k<_dom.N2();k++)
			{	
				for(int l=NV; l<_dom.Nval();l++)
				{
					_V[i][j][k][l] = _V[i][nBound][k][l];
				}
			}
		}
	}
	
	/* back */
	for(int j=_dom.N2()-nBound;j<_dom.N2();j++)
	{
		for(int i=0;i<_dom.N2();i++)
		{	
			for(int k=0;k<_dom.N3();k++)
			{
				for(int l=NV; l<_dom.Nval();l++)
				{
					_V[i][j][k][l] = _V[i][_dom.N2()-nBound-1][k][l];
				}
			}
		}
	}
	/* top */
	for(int k=0;k<nBound;k++)
	{	
		for(int j=0;j<_dom.N2();j++)
		{
			for(int i=0;i<_dom.N1();i++)
			{
				for(int l=NV; l<_dom.Nval();l++)
				{	
					_V[i][j][k][l] = _V[i][j][nBound][l];
				}
			}
		}
	}


	/* bottom */
	for(int k=_dom.N3()-nBound;k<_dom.N3();k++)
	{	
		for(int j=0;j<_dom.N2();j++)
		{
			for(int i=0;i<_dom.N1();i++)
			{
				for(int l=NV; l<_dom.Nval();l++)
				{
					_V[i][j][k][l] = _V[i][j][_dom.N3()-nBound-1][l];
				}
			}
		}
	}
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

void Grid::coefs_Laplacian(int nBound, vector<double>& fcx, vector<double>& fcy, vector<double>& fcz, double& cc) const
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

Grid Grid::laplacian(int nBound) const
{
	try
	{
		if(nBound<1 or nBound>12)
		{
			throw string("nBound oustide acceptable precision values");
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
						double v = cc*_V[i][j][k][0];
						for(int n=1; n<=nBound ; n++)
						{
							v += fcx[n] *(_V[i-n][j][k][0]+_V[i+n][j][k][0]);
							v += fcy[n] *(_V[i][j-n][k][0]+_V[i][j+n][k][0]);
							v += fcz[n] *(_V[i][j][k-n][0]+_V[i][j][k+n][0]);
						}
						g._V[i][j][k][0]=v;
					}
				}
			}
			g.reset_Boundary(nBound);
			return g;
		}
	}
	catch(string error)
	{
		cout<<error<<endl;
		exit(1);
	}				
}

void Grid::coefs_Gradient(int nBound, vector<double>& fcx, vector<double>& fcy, vector<double>& fcz) const
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

		for(int i=0;i<nBound;i++)
		{
			fcx[i] =  coefs[i]/_dom.dx();
			fcy[i] =  coefs[i]/_dom.dy();
			fcz[i] =  coefs[i]/_dom.dz();
		}
}


Grid Grid::gradient(int nBound) const
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
			Domain dg=_dom;
			dg.set_Nval(4);
			Grid g(dg);
			g._str=_str;
			coefs_Gradient(nBound, fcx, fcy, fcz);
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
						g._V[i][j][k][0]=_V[i][j][k][0];
						double gx=0;
						double gy=0;
						double gz=0;
						for(int n=-nBound, kn=0 ; kn<nBound ; n++, kn++)
						{
							gx += fcx[kn] * (_V[i+n][j][k][0]-_V[i-n][j][k][0]);
							gy += fcy[kn] * (_V[i][j+n][k][0]-_V[i][j-n][k][0]);
							gz += fcz[kn] * (_V[i][j][k+n][0]-_V[i][j][k-n][0]);
						}
						g._V[i][j][k][1]=gx;
						g._V[i][j][k][2]=gy;
						g._V[i][j][k][3]=gz;
					}
				}
			}
			g.reset_Boundary(nBound);
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
	if(_dom.Nval()==1)
	{
		nameFile<<"Grid generated by CdftT "<<endl;
		nameFile<<"Density"<<endl;
	}
	else if(_dom.Nval()==4)
	{
		nameFile<<"Grid generated by CdftT "<<endl;
		nameFile<<"Gradient"<<endl;
	}
	else
	{
		nameFile<<"Grid generated by CdftT "<<endl;
		nameFile<<"Orbitals"<<endl;
	}
	nameFile<<scientific;
	nameFile<<_str.number_of_atoms()<<" ";
	for(int i=0;i<3;i++)
	{
		nameFile<<_dom.O()[i]<<" ";
	}
	nameFile<<_dom.Nval()<<" ";
	nameFile<<endl<<_dom.N1()<<" ";
	if(_dom.N1()>0)
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
	for(int i=0;i<_str.number_of_atoms();i++)
	{
		nameFile<<_str.atom(i).atomic_number()<<" ";
		nameFile<<_str.atom(i).charge()<<" ";
		for(int j=0; j<3;j++)
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
double Grid::value(int i, int j, int k) const
{
	return _V[i][j][k][0];
}
double Grid::value(int i, int j, int k, int l) const
{
	return _V[i][j][k][l];
}

void Grid::next(int i, int j, int k, double& current, vector<vector<int>>& trajectory)
{
	vector<int> v=trajectory.back();
	int I[3];
	int i2 =i-1;
	int i1 =i+1;
	I[0] = i2;
	I[1] = i;
	I[2] = i1;
		
	int J[3];
	int j1 = j+1;
	int j2 = j-1;
	J[0] = j2;
	J[1] = j;
	J[2] = j1;
			
	int K[3];
	int k1 = k+1;
	int k2 = k-1;
	K[0] = k2;
	K[1] = k;
	K[2] = k1;
	if(i2<0) i2 = i;
	if(i1>_dom.N1()-1) i1 = i;
	if(j2<0) j2 = j;
	if(j1>_dom.N2()-1) j1 = j;
	if(k2<0) k2 = k;
	if(k1>_dom.N3()-1) k1 = k;
	for(int ic=0;ic<3;ic++)
	{
		for(int jc=0;jc<3;jc++)
		{	
				
			for(int kc=0;kc<3;kc++)
			{
				if(ic==1 and jc==1 and kc==1) continue;
				vector<double> ds(3);
				ds[0]=(I[ic]-I[1])*_dom.dx();
				ds[1]=(J[jc]-J[1])*_dom.dy();
				ds[2]=(K[kc]-K[1])*_dom.dz();
				double grad =0;
				double normds=sqrt(ds[0]*ds[0]+ds[1]*ds[1]+ds[2]*ds[2]);
				for(int m=1;m<=3;m++)
				{
					grad +=_V[I[ic]][J[jc]][K[kc]][m]*ds[m-1];
				}
				grad=grad/normds;
				if(grad>current)
				{
					v[0]=I[ic];
					v[1]=J[jc];
					v[2]=K[kc];
					current = grad;
				}
			}
		}
	}
	
	trajectory.push_back(v);
}

vector<double> Grid::atom_attract_diff(const vector<vector<int>>& attract)
{
	vector<double> v(_str.number_of_atoms());
	double distance=0;
	double d1=100;
	vector<double> ds(3);
	ds[0]=_dom.dx();
	ds[1]=_dom.dy();
	ds[2]=_dom.dz();
	for(int j=0; j<_str.number_of_atoms();j++)
	{
		for(int n=0;n<int(attract.size()); n++)
		{
			for(int i=0;i<3;i++)
			{
				distance += double((_str.atoms()[j].coordinates()[i]-attract[n][i])*(_str.atoms()[j].coordinates()[i]-attract[n][i]))*ds[i]*ds[i];
			}
			distance=sqrt(distance);
			if(distance<d1)
			{
				d1=distance;
			}	
		}
		v[j]=d1;
		
	}
	return v;
}


void Grid::addSurroundingEqualPoints(int i,int j,int k, vector<vector<int>>& equals, double& current)
{
	vector<int> v(3);
	int I[3];
	int i2 =i-1;
	int i1 =i+1;
	I[0] = i2;
	I[1] = i;
	I[2] = i1;
	
	int J[3];
	int j1 = j+1;
	int j2 = j-1;
	J[0] = j2;
	J[1] = j;
	J[2] = j1;
	
	int K[3];
	int k1 = k+1;
	int k2 = k-1;
	K[0] = k2;
	K[1] = k;
	K[2] = k1;
	for(int ic=0;ic<3;ic++)
	{	
		for(int jc=0;jc<3;jc++)
		{
			for(int kc=0;kc<3;kc++)
			{
				if(ic==1 and jc==1 and kc==1) continue;
				if(i2<0) i2 = i;
				if(i1>_dom.N1()-1) i1 = i;
				if(j2<0) j2 = j;
				if(j1>_dom.N2()-1) j1 = j;
				if(k2<0) k2 = k;
				if(k1>_dom.N3()-1) k1 = k;
				vector<double> ds(3);
				ds[0]=(I[ic]-I[1])*_dom.dx();
				ds[1]=(J[jc]-J[1])*_dom.dy();
				ds[2]=(K[kc]-K[1])*_dom.dz();
				double grad=0;
				double normds=sqrt((ds[0]*ds[0])+(ds[1]*ds[1])+(ds[2]*ds[2]));
				for(int m=1;m<3;m++)
				{
					grad +=_V[I[1]][J[1]][K[1]][m]*ds[m-1];
				}
				grad=grad/normds;
				if(grad==current)	
				{
					v[0]=I[ic];
					v[1]=J[jc];
					v[2]=K[kc];
					equals.push_back(v);
				}	
			}
		}
	}
}
Grid Grid::aim_On_Grid(int nBound)
{
	try
	{
		if(_dom.Nval()==4)
		{
			Grid g(Domain(1,_dom.N1(), _dom.N2(), _dom.N3(), _dom.O()));
			for(int i=0;i<3;i++)
			{
				for(int j=0;j<3;j++)
				{
					g._dom.set_T(_dom.Tij(i,j), i, j);
				}
			}
			g._str=_str;
			//make list instead
			vector<vector<int>> attractors(0, vector<int>(3));
			cout<<_V[36][29][32][2]<<endl;
			cout<<_V[36][30][32][2]<<endl;
			for(int i=nBound;i<_dom.N1()-nBound;i++)
			{
				for(int j=nBound;j<_dom.N2()-nBound;j++)
				{
					for(int k=nBound;k<_dom.N3()-nBound;k++)
					{
						if(g._V[i][j][k][0]>PRECISION)
						{
							continue;
						}
						else
						{	
							double current=0;
							bool KnownPt=false;
							vector<vector<int>> trajectory={{i,j,k}};
							vector<vector<int>> equals(0, vector<int>(3));
							do
							{
								current=0;
								int I=trajectory.back()[0];
								int J=trajectory.back()[1];
								int K=trajectory.back()[2];
								next(I,J,K,current,trajectory);
								
								
								cout<<trajectory.back()[0]<<endl;
								cout<<trajectory.back()[1]<<endl;
								cout<<trajectory.back()[2]<<endl;
								if(g._V[trajectory.back()[0]][trajectory.back()[1]][trajectory.back()[2]][0]>PRECISION)
								{
									KnownPt=true;
									break;
								}
							} while(current>PRECISION);
							for(int p=0;p<int(trajectory.size()-1);p++)
							{
								g._V[trajectory[p][0]][trajectory[p][1]][trajectory[p][2]][0]=g._V[trajectory.back()[0]][trajectory.back()[1]][trajectory.back()[2]][0];
							}
							addSurroundingEqualPoints(trajectory.back()[0],trajectory.back()[1],trajectory.back()[2],equals,current);
							for(int p=0;p<int(equals.size());p++)
							{
								g._V[equals[p][0]][equals[p][1]][equals[p][2]][0]=g._V[trajectory.back()[0]][trajectory.back()[1]][trajectory.back()[2]][0];	
							}
							//if point has been already indexed exit loop: new point
							if(KnownPt)
							{
								cout<<"gone"<<endl;
								continue;	
							}
							
							//Do we know max?
							bool WhoIsMax=true;
							//compare max to list
							cout<<i<<endl;
							
							for(int m=0;m<int(attractors.size());m++)
							{
								int I=1;
								
								if(attractors[m]==trajectory.back())
								{
									WhoIsMax=false;
									for(int p=0;p<int(trajectory.size());p++)
									{
										g._V[trajectory[p][0]][trajectory[p][1]][trajectory[p][2]][0]=double(I);
									}
								}
								I++;
							}
							//dont know max
							//add to known attractors
							if(WhoIsMax)
							{	
								if(int(trajectory.size())<(_dom.N1()*_dom.N2()*_dom.N3())/1000)
								{
									
									for(int p=0;p<int(equals.size());p++)
									{
										g._V[equals[p][0]][equals[p][1]][equals[p][2]][0]=double(attractors.size());		
									}
									for(int p=0;p<int(trajectory.size());p++)
									{
										g._V[trajectory[p][0]][trajectory[p][1]][trajectory[p][2]][0]=double(attractors.size());
									}
								}
								attractors.push_back(trajectory.back());
								//index all points from trajectory
								cout<<"traj size "<<trajectory.size()<<endl;
								cout<<-_V[36][29][32][2]<<endl;
								cout<<"attract size "<<attractors.size()<<endl;
								for(int p=0;p<int(equals.size());p++)
								{
									g._V[equals[p][0]][equals[p][1]][equals[p][2]][0]=g._V[trajectory.back()[0]][trajectory.back()[1]][trajectory.back()[2]][0];		
								}
								for(int p=0;p<int(trajectory.size());p++)
								{
									g._V[trajectory[p][0]][trajectory[p][1]][trajectory[p][2]][0]=double(attractors.size());
								}
							}
						
						}					
					}
				}
			}
			/*
			vector<vector<int>> vpos(0,vector<int>(3));
			for(int p=0;p<int(attractors.size()); p++)
			{
				for(int m=0;m<int(attractors.size());m++)
				{
					int I[2];
					I[0]=-1;
					I[1]=1;
					for(int i=0;i<2;i++)
					{
						if(attractors[p][0]==attractors[m][0] +I[i] and attractors[p][1]==attractors[m][1] and attractors[p][2]==attractors[m][2])
						{
							vpos.push_back(attractors[m]);
							for(int l=nBound;l<_dom.N1()-nBound;l++)
							{
								for(int j=nBound;j<_dom.N2()-nBound;j++)
								{
									for(int k=nBound;k<_dom.N3()-nBound;k++)
									{
										if(g._V[l][j][k][0]==m)
										{
											cout<<"okay"<<endl;
											g._V[l][j][k][0]=p;
											
										}
									}
								}
							}
						}
						if( attractors[p][0]==attractors[m][0] and attractors[p][1]==attractors[m][1] +I[i] and attractors[p][1]==attractors[m][1])
						{
							vpos.push_back(attractors[m]);
							for(int l=nBound;l<_dom.N1()-nBound;l++)
							{
								for(int j=nBound;j<_dom.N2()-nBound;j++)
								{
									for(int k=nBound;k<_dom.N3()-nBound;k++)
									{
										if(g._V[l][j][k][0]==m)
										{
											cout<<"okay"<<endl;
											g._V[l][j][k][0]=p;
											
										}
									}
								}
							}
						}
						if(attractors[p][0]==attractors[m][0] and attractors[p][1]==attractors[m][1] and attractors[p][2]==attractors[m][2] +I[i])
						{
							vpos.push_back(attractors[m]);
							for(int l=nBound;l<_dom.N1()-nBound;l++)
							{
								for(int j=nBound;j<_dom.N2()-nBound;j++)
								{
									for(int k=nBound;k<_dom.N3()-nBound;k++)
									{
										if(g._V[l][j][k][0]==m)
										{
											cout<<"okay"<<endl;
											g._V[l][j][k][0]=p;
											
										}
									}
								}
							}
						}
					}
				}
			}
			for(int k=0;k<int(vpos.size());k++)
			{
				for(auto element:attractors)
				{	
					auto it=find(attractors.begin(),attractors.end(), vpos[k]);
					if(it!=attractors.end()) attractors.erase(it);
				}
			}*/
			for(int m=0; m<int(attractors.size());m++)
			{
				vector<double> ds(3);
				ds[0]=_dom.dx();
				ds[1]=_dom.dy();
				ds[2]=_dom.dz();
				cout<<"attractor "<<m<<endl;
				for(int a=0;a<3;a++)
				{
					cout<<(attractors[m][a]*ds[a] + _dom.O()[a])*BOHRTOANG<<endl;
				}
			}
			for(int m=0; m<_str.number_of_atoms();m++)
			{
				cout<<atom_attract_diff(attractors)[m]<<endl;
				cout<<atom_attract_diff(attractors).size()<<endl;
			}
			
			
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

void Grid::next_Density(int i, int j, int k, double& rhocenter, vector<vector<int>>& trajectory)
{
	//coordinates on grid of the latest point
	vector<int> v=trajectory.back();
	
	// Initialize indices for neighboring points
	int I[3] = { max(0, i-1), i, min(i+1, _dom.N1()-1) };
	int J[3] = { max(0, j-1), j, min(j+1, _dom.N2()-1) };
	int K[3] = { max(0, k-1), k, min(k+1, _dom.N3()-1) };
	
	double drhomax=0;
	vector<double> ds(3);
	
	// Iterate through neighboring cells
	for(int ic=0;ic<3;ic++)
	{
		for(int jc=0;jc<3;jc++)
		{	
				
			for(int kc=0;kc<3;kc++)
			{
				// Skip the center point
				if(ic==1 and jc==1 and kc==1) continue;
				
				// Calculate distances
				ds[0]=(I[ic]-I[1])*_dom.dx();
				ds[1]=(J[jc]-J[1])*_dom.dy();
				ds[2]=(K[kc]-K[1])*_dom.dz();
				
				// Calculate the norm of the distance vector
				double normds=sqrt(ds[0]*ds[0]+ds[1]*ds[1]+ds[2]*ds[2]);
				
				// Get the density at the neighboring point
				double rho=_V[I[ic]][J[jc]][K[kc]][0];
				
				// Calculate the density gradient
				double drho = (rho-rhocenter)/normds;
				
				// Check if this is the maximum density gradient
				if(drho-drhomax>PRECISION)
				{
					v[0]=I[ic];
					v[1]=J[jc];
					v[2]=K[kc];
					drhomax = drho;
				}
			}
		}
	}
	
	//Update density for the new point
	rhocenter=_V[v[0]][v[1]][v[2]][0];
	
	//Add the new point to the ascent trajectory
	trajectory.push_back(v);
}

void Grid::addSurroundingDensity(int i,int j,int k, vector<vector<int>>& equals, double& rhocenter)
{
	//coordinates on grid of the latest point
	vector<int> v(3);
	
	// Initialize indices for neighboring points
	int I[3] = { max(0, i-1), i, min(i+1, _dom.N1()-1) };
	int J[3] = { max(0, j-1), j, min(j+1, _dom.N2()-1) };
	int K[3] = { max(0, k-1), k, min(k+1, _dom.N3()-1) };
	
	for(int ic=0;ic<3;ic++)
	{	
		for(int jc=0;jc<3;jc++)
		{
			for(int kc=0;kc<3;kc++)
			{
				// Skip the center point
				if(ic==1 and jc==1 and kc==1) 
				{
					continue;
				}
				
				// Calculate distances
				vector<double> ds(3);
				ds[0]=(I[ic]-I[1])*_dom.dx();
				ds[1]=(J[jc]-J[1])*_dom.dy();
				ds[2]=(K[kc]-K[1])*_dom.dz();
				double normds=sqrt((ds[0]*ds[0])+(ds[1]*ds[1])+(ds[2]*ds[2]));
				
				// Get the density at the neighboring point
				double rho=_V[I[ic]][J[jc]][K[kc]][0];
				
				// Calculate the density gradient
				rho =(rho-rhocenter)/normds;
				
				// Check if the gradient is zero
				if(rho<PRECISION)
				{
					v[0]=I[ic];
					v[1]=J[jc];
					v[2]=K[kc];
					
					// If it is, add point to list
					equals.push_back(v);
				}	
			}
		}
	}
}

Grid Grid::aim_On_Grid_Density()
{

	try
	{
		if(_dom.Nval()==1)
		{
			// Initialise new grid for indices
			Grid g(Domain(1,_dom.N1(), _dom.N2(), _dom.N3(), _dom.O()));
			for(int i=0;i<3;i++)
			{
				for(int j=0;j<3;j++)
				{
					g._dom.set_T(_dom.Tij(i,j), i, j);
				}
			}
			g._str=_str;
			//make list instead
			// Create a list of the attractors' coordinates
			vector<vector<int>> attractors(0, vector<int>(3));
			for(int i=0;i<_dom.N1();i++)
			{
				for(int j=0;j<_dom.N2();j++)
				{
					for(int k=0;k<_dom.N3();k++)
					{
						// Check if the point has already been indexed
						if(g._V[i][j][k][0]>PRECISION)
						{
							cout<<"already indexed"<<endl;
							continue;
						}
						else
						{	
							// Initialise density at point i j k and the known point flag
							double rhocenter=_V[i][j][k][0];
							bool KnownPt=false;
							
							// Initialise ascent trajectory and equal density neighbours
							vector<vector<int>> trajectory={{i,j,k}};
							vector<vector<int>> equals(0, vector<int>(3));
							
							
							// Start the loop to find the path of highest density gradient
							do
							{
								//Initialize indices to the point
								int I=trajectory.back()[0];
								int J=trajectory.back()[1];
								int K=trajectory.back()[2];

								next_Density(I,J,K,rhocenter,trajectory);
								
								// Update indices to the new point
								int newI = trajectory.back()[0];
								int newJ = trajectory.back()[1];
								int newK = trajectory.back()[2];
								cout<<newI<<endl;
								cout<<newJ<<endl;
								cout<<newK<<endl;
								cout<<rhocenter<<endl;
								
								// Check if a new point was added to the trajectory
								if(I==newI and J==newJ and K==newK)
								{
									cout<<"no new point found"<<endl;
									break;
								}
								
								//check if the new point is already indexed
								if(g._V[I][J][K][0]>PRECISION)
								{
									KnownPt=true;
									break;
								}
							} while(true);
							
							// Iterate through the trajectory and set the index value of each point in the trajectory to the index of the associated maximum (last point)
							for(int p=0;p<int(trajectory.size()-1);p++)
							{
								g._V[trajectory[p][0]][trajectory[p][1]][trajectory[p][2]][0]=g._V[trajectory.back()[0]][trajectory.back()[1]][trajectory.back()[2]][0];
							}
							
							
							// Add points surrounding the maximum with equal density gradients to equals for indexing
							addSurroundingDensity(trajectory.back()[0],trajectory.back()[1],trajectory.back()[2],equals,rhocenter);
							
							
							// Iterate through the equal points and set the index value of each point in the the vector to the index of the associated maximum 
							for(int p=0;p<int(equals.size());p++)
							{
								g._V[equals[p][0]][equals[p][1]][equals[p][2]][0]=g._V[trajectory.back()[0]][trajectory.back()[1]][trajectory.back()[2]][0];	
							}
							
							
							// If the point has been already indexed move on to a new point
							if(KnownPt)
							{
								cout<<"gone"<<endl;
								continue;	
							}
							
							// Known attractor flag
							bool WhoIsMax=true;

							// Compare new attractor to known attractors. If the attractor is known index all points from trajectory and change flag value
							for(int m=0;m<int(attractors.size());m++)
							{
								int I=1;
								
								if(attractors[m]==trajectory.back())
								{
									WhoIsMax=false;
									for(int p=0;p<int(trajectory.size());p++)
									{
										g._V[trajectory[p][0]][trajectory[p][1]][trajectory[p][2]][0]=double(I);
									}
								}
								I++;
							}
							
							// If the attractor is unknown
							if(WhoIsMax)
							{	
							
								//if the number of points leading to the attractor is insufficient, index trajectory and equals to previous attractor and disregard current
								if(int(trajectory.size())<(_dom.N1()*_dom.N2()*_dom.N3())/1000)
								{
									
									for(int p=0;p<int(equals.size());p++)
									{
										g._V[equals[p][0]][equals[p][1]][equals[p][2]][0]=double(attractors.size());		
									}
									for(int p=0;p<int(trajectory.size());p++)
									{
										g._V[trajectory[p][0]][trajectory[p][1]][trajectory[p][2]][0]=double(attractors.size());
									}
									
								}
								
								// Add the attractor to the list of attractors
								attractors.push_back(trajectory.back());
								
								
								cout<<"traj size "<<trajectory.size()<<endl;
								cout<<"attract size "<<attractors.size()<<endl;
								
								// Index all points from trajectory and equals
								for(int p=0;p<int(equals.size());p++)
								{
									g._V[equals[p][0]][equals[p][1]][equals[p][2]][0]=g._V[trajectory.back()[0]][trajectory.back()[1]][trajectory.back()[2]][0];		
								}
								for(int p=0;p<int(trajectory.size());p++)
								{
									g._V[trajectory[p][0]][trajectory[p][1]][trajectory[p][2]][0]=double(attractors.size());
								}
							}
						
						}					
					}
				}
			}
			/*
			vector<vector<int>> vpos(0,vector<int>(3));
			for(int p=0;p<int(attractors.size()); p++)
			{
				for(int m=0;m<int(attractors.size());m++)
				{
					int I[2];
					I[0]=-1;
					I[1]=1;
					for(int i=0;i<2;i++)
					{
						if(attractors[p][0]==attractors[m][0] +I[i] and attractors[p][1]==attractors[m][1] and attractors[p][2]==attractors[m][2])
						{
							vpos.push_back(attractors[m]);
							for(int l=nBound;l<_dom.N1()-nBound;l++)
							{
								for(int j=nBound;j<_dom.N2()-nBound;j++)
								{
									for(int k=nBound;k<_dom.N3()-nBound;k++)
									{
										if(g._V[l][j][k][0]==m)
										{
											cout<<"okay"<<endl;
											g._V[l][j][k][0]=p;
											
										}
									}
								}
							}
						}
						if( attractors[p][0]==attractors[m][0] and attractors[p][1]==attractors[m][1] +I[i] and attractors[p][1]==attractors[m][1])
						{
							vpos.push_back(attractors[m]);
							for(int l=nBound;l<_dom.N1()-nBound;l++)
							{
								for(int j=nBound;j<_dom.N2()-nBound;j++)
								{
									for(int k=nBound;k<_dom.N3()-nBound;k++)
									{
										if(g._V[l][j][k][0]==m)
										{
											cout<<"okay"<<endl;
											g._V[l][j][k][0]=p;
											
										}
									}
								}
							}
						}
						if(attractors[p][0]==attractors[m][0] and attractors[p][1]==attractors[m][1] and attractors[p][2]==attractors[m][2] +I[i])
						{
							vpos.push_back(attractors[m]);
							for(int l=nBound;l<_dom.N1()-nBound;l++)
							{
								for(int j=nBound;j<_dom.N2()-nBound;j++)
								{
									for(int k=nBound;k<_dom.N3()-nBound;k++)
									{
										if(g._V[l][j][k][0]==m)
										{
											cout<<"okay"<<endl;
											g._V[l][j][k][0]=p;
											
										}
									}
								}
							}
						}
					}
				}
			}
			for(int k=0;k<int(vpos.size());k++)
			{
				for(auto element:attractors)
				{	
					auto it=find(attractors.begin(),attractors.end(), vpos[k]);
					if(it!=attractors.end()) attractors.erase(it);
				}
			}*/
			for(int m=0; m<int(attractors.size());m++)
			{
				vector<double> ds(3);
				ds[0]=_dom.dx();
				ds[1]=_dom.dy();
				ds[2]=_dom.dz();
				cout<<"attractor "<<m<<endl;
				for(int a=0;a<3;a++)
				{
					cout<<(attractors[m][a]*ds[a]+_dom.O()[a])*BOHRTOANG<<endl;
					cout<<attractors[m][a]<<endl;
				}
			}
			for(int m=0; m<_str.number_of_atoms();m++)
			{
				cout<<atom_attract_diff(attractors)[m]<<endl;
				cout<<atom_attract_diff(attractors).size()<<endl;
			}
			//int nbrofindices=0;
			double value=0;
			cout<<"V0000"<<g._V[0][0][0][0]<<endl;
			for(int i=0;i<_dom.N1();i++)
			{
				for(int j=0;j<_dom.N2();j++)
				{
					for(int k=0;k<_dom.N3();k++)
					{
						if(g._V[i][j][k][0]>value)
						{
							value =g._V[i][j][k][0];
						}
					}
				}
			}
			cout<<"number of of idices "<<int(value)<<endl;
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
