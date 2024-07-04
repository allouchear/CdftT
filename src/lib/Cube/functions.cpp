using namespace std;
#include <Cube/functions.h>
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

