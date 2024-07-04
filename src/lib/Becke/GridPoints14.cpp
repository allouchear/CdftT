#include<iostream>
#include <Becke/GridPoints.h>

using namespace std;

void GridPoints::GridPoints14()
{
	_Npts = 14;
	_Lmax = 5;
	_L2max = 2;
	_LebedevGridPoints = {
		// (   theta,            phi,            weight) 
    	{  1.570796326795,  0.000000000000,  0.066666666667},
        {  1.570796326795,  3.141592653590,  0.066666666667},
        {  1.570796326795,  1.570796326795,  0.066666666667},
        {  1.570796326795,  4.712388980385,  0.066666666667},
        {  0.000000000000,  3.141592653590,  0.066666666667},
        {  3.141592653590,  3.141592653590,  0.066666666667},
        {  0.955316618125,  0.785398163397,  0.075000000000},
        {  0.955316618125,  2.356194490192,  0.075000000000},
        {  0.955316618125,  5.497787143782,  0.075000000000},
        {  0.955316618125,  3.926990816987,  0.075000000000},
        {  2.186276035465,  0.785398163397,  0.075000000000},
        {  2.186276035465,  2.356194490192,  0.075000000000},
        {  2.186276035465,  5.497787143782,  0.075000000000},
        {  2.186276035465,  3.926990816987,  0.075000000000}
	};
}