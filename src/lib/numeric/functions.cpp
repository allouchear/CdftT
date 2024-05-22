using namespace std;
#include <numeric/functions.h>
#include <cmath>
#include <vector>


double coulomb(vector<double> fpar, double x, double y, double z, const Grid& g)
{
	double normO=0;
	double r=sqrt(x*x+y*y+z*z);
	for(int i=0; i<3;i++)
	{
		normO+= g.dom().O()[i]*g.dom().O()[i];
	}
	normO= sqrt(normO);
	if(r==normO)
		return 0;
	else
		return fpar[0]/( r- normO);
}

double laplacian(vector<double> fpar, double x, double y, double z, const Grid& g)
{
	double dx=0;
	double dy=0;
	double dz=0;
	for(int i=0; i<3;i++)
	{
		dx += g.dom().T()[0][i]*g.dom().T()[0][i];
		dy += g.dom().T()[1][i]*g.dom().T()[1][i];
		dz += g.dom().T()[2][i]*g.dom().T()[2][i];
	}
	dx=sqrt(dx);
	dy=sqrt(dy);
	dz=sqrt(dz);
	double L=0;
	switch(int(fpar[0]))
	{	
		case 1:
		{
			vector<double> c = {-2.0, 1.0};
			L = c[0]*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0]+c[1]*g.V()[int(x/dx)+1][int(y/dy)][int(z/dz)][0] + c[1]*g.V()[int(x/dx)-1][int(y/dy)][int(z/dz)][0];
			
			break;
		}
		case 2:
		{
			double denom = 12.0;
			vector<double> c = {-30.0, 16.0, -1.0};
			for(int i=0;i<=2;i++)
				L = (c[i]/(denom*dx))*g.V()[int(x/dx)+i][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dx))*g.V()[int(x/dx)-i][int(y/dy)][int(z/dz)][0]-(c[i]/(denom*dx))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)+i][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)-i][int(z/dz)][0] - (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)+i][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz) - i][0]-(c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0];
			break;
		}
		case 3:
		{
			double denom = 180.0;
			vector<double> c = {-490.0, 270.0,-27.0, 2.0};
			for(int i=0;i<=3;i++)
				L = (c[i]/(denom*dx))*g.V()[int(x/dx)+i][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dx))*g.V()[int(x/dx)-i][int(y/dy)][int(z/dz)][0]-(c[i]/(denom*dx))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)+i][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)-i][int(z/dz)][0] - (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)+i][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz) - i][0]-(c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0];
			break;
		}
		case 4:
		{
			double denom = 5040.0;
			vector<double> c = {-14350.0, 8064.0, -1008.0, 128.0, -9.0};
			for(int i=0;i<=4;i++)
				L = (c[i]/(denom*dx))*g.V()[int(x/dx)+i][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dx))*g.V()[int(x/dx)-i][int(y/dy)][int(z/dz)][0]-(c[i]/(denom*dx))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)+i][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)-i][int(z/dz)][0] - (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)+i][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz) - i][0]-(c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0];
			break;
		}
		case 5:
		{
			double denom = 25200.0;
			vector<double> c = {-73766.0, 42000.0, -6000.0, 1000.0, -125.0, 8.0};
			for(int i=0;i<=5;i++)
				L = (c[i]/(denom*dx))*g.V()[int(x/dx)+i][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dx))*g.V()[int(x/dx)-i][int(y/dy)][int(z/dz)][0]-(c[i]/(denom*dx))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)+i][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)-i][int(z/dz)][0] - (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)+i][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz) - i][0]-(c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0];
			break;
		}
		case 6:
		{
			double denom = 831600.0;
			vector<double> c = {-2480478.0,1425600.0,-222750.0,44000.0,-7425.0,864.0,-50.0};
			for(int i=0;i<=6;i++)
				L = (c[i]/(denom*dx))*g.V()[int(x/dx)+i][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dx))*g.V()[int(x/dx)-i][int(y/dy)][int(z/dz)][0]-(c[i]/(denom*dx))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)+i][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)-i][int(z/dz)][0] - (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)+i][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz) - i][0]-(c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0];
			break;
		}
		case 7:
		{
			double denom = 75675600.0;
			vector<double> c = {-228812298.0,132432300.0,-22072050.0,4904900.0,-1003275.0, 160524.0,-17150.0,900.0};
			for(int i=0;i<=7;i++)
				L = (c[i]/(denom*dx))*g.V()[int(x/dx)+i][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dx))*g.V()[int(x/dx)-i][int(y/dy)][int(z/dz)][0]-(c[i]/(denom*dx))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)+i][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)-i][int(z/dz)][0] - (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)+i][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz) - i][0]-(c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0];
			break;
		}
		case 8:
		{
			double denom = 302702400.0;
			vector<double> c = {-924708642.0,538137600.0,-94174080.0,22830080.0,-5350800.0,1053696.0,-156800.0,15360.0,-735.0};
			for(int i=0;i<=8;i++)
				L = (c[i]/(denom*dx))*g.V()[int(x/dx)+i][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dx))*g.V()[int(x/dx)-i][int(y/dy)][int(z/dz)][0]-(c[i]/(denom*dx))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)+i][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)-i][int(z/dz)][0] - (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)+i][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz) - i][0]-(c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0];
			break;
		}
		case 9:
		{
			double denom = 15437822400.0;
			vector<double> c = {-47541321542.0,+27788080320.0, -5052378240.0,+1309875840.0,-340063920.0,+77728896.0,-14394240.0,+1982880.0,-178605.0,+7840.0};
			for(int i=0;i<=9;i++)
				L = (c[i]/(denom*dx))*g.V()[int(x/dx)+i][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dx))*g.V()[int(x/dx)-i][int(y/dy)][int(z/dz)][0]-(c[i]/(denom*dx))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)+i][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)-i][int(z/dz)][0] - (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)+i][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz) - i][0]-(c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0];
			break;
		}
		case 10:
		{
			double denom = 293318625600.0;
			vector<double> c = {-909151481810.0,+533306592000.0, -99994986000.0,+27349056000.0,-7691922000.0,+1969132032.0,-427329000.0,+73872000.0,-9426375.0,+784000.0,-31752.0};
			for(int i=0;i<=10;i++)
				L = (c[i]/(denom*dx))*g.V()[int(x/dx)+i][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dx))*g.V()[int(x/dx)-i][int(y/dy)][int(z/dz)][0]-(c[i]/(denom*dx))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)+i][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)-i][int(z/dz)][0] - (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)+i][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz) - i][0]-(c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0];
			break;
		}
		case 11:
		{
			double denom = 3226504881600.0;
			vector<double> c = {-10053996959110.0,+5915258949600.0,-1137549798000.0,+325014228000.0,-97504268400.0,+27301195152.0,-6691469400.0,+1365606000.0,-220114125.0,+26087600.0,-2012472.0,+75600.0};
			for(int i=0;i<=11;i++)
				L = (c[i]/(denom*dx))*g.V()[int(x/dx)+i][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dx))*g.V()[int(x/dx)-i][int(y/dy)][int(z/dz)][0]-(c[i]/(denom*dx))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)+i][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)-i][int(z/dz)][0] - (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)+i][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz) - i][0]-(c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0];
			break;
		}
		case 12:
		{
			double denom = 74209612276800.0;
			vector<double> c = {-232272619118930.0,+137002361126400.0,-26911178078400.0,+7973682393600.0,-2522922944850.0,+759845028096.0,-205205061600.0,+47609337600.0,-9112724775.0,+1371462400.0,-151484256.0,+10886400.0,-381150.0};
			for(int i=0;i<=12;i++)
				L = (c[i]/(denom*dx))*g.V()[int(x/dx)+i][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dx))*g.V()[int(x/dx)-i][int(y/dy)][int(z/dz)][0]-(c[i]/(denom*dx))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)+i][int(z/dz)][0] + (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)-i][int(z/dz)][0] - (c[i]/(denom*dy))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)+i][0] + (c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz) - i][0]-(c[i]/(denom*dz))*g.V()[int(x/dx)][int(y/dy)][int(z/dz)][0];
			break;
		}
	}
	//call reset boundaries function
	
	return L;
}
/*
double grad(vector<double>, double x, double y, double z, const Grid& g)
{
	
}*/
