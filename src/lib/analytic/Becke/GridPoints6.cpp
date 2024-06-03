#include<iostream>
#include<analytic/Becke/GridPoints.h>

using namespace std;

void GridPoints::GridPoints6()
{
	_Lebedev_Npts = 6;
	_Lebedev_Lmax = 3;
	_Lebedev_L2max = 1;
	_LebedevGridPoints = {
		// (   theta,            phi,            weight) 
    	{  1.570796326795,  0.000000000000,  0.166666666667},
    	{  1.570796326795,  3.141592653590,  0.166666666667},
    	{  1.570796326795,  1.570796326795,  0.166666666667},
    	{  1.570796326795,  4.712388980385,  0.166666666667},
    	{  0.000000000000,  3.141592653590,  0.166666666667},
    	{  3.141592653590,  3.141592653590,  0.166666666667}
	};
}