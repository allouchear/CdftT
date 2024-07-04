#include<iostream>
#include <Becke/GridPoints.h>

using namespace std;

GridPoints::GridPoints()
{
	_Lebedev_Npts={6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 5810};
	_Lebedev_Lmax={3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53, 59, 65, 71, 77, 83, 131};
	_Lebedev_L2max={1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 20, 23, 26, 29, 32, 35, 38, 41, 65};
	_LebedevGridPoints=vector<vector<double>> ();
}

GridPoints::GridPoints(int Lmax)
{
	if(Lmax==3)
		GridPoints6();
	else if(Lmax==5)
		GridPoints14();
	else if(Lmax==7)
		GridPoints26();
	else if(Lmax==9)
		GridPoints38();
	else if(Lmax==11)
		GridPoints50();
	else if(Lmax==13)
		GridPoints74();
	else if(Lmax==15)
		GridPoints86();
	else if(Lmax==17)
		GridPoints110();
	else if(Lmax==19)
		GridPoints146();
	else if(Lmax==21)
		GridPoints170();
	else if(Lmax==23)
		GridPoints194();
	else if(Lmax==25)
		GridPoints230();
	else if(Lmax==27)
		GridPoints266();
	else if(Lmax==29)
		GridPoints302();
	else if(Lmax==31)
		GridPoints350();
	else if(Lmax==35)
		GridPoints434();
	else if(Lmax==41)
		GridPoints590();
	else if(Lmax==47)
		GridPoints770();
	else if(Lmax==53)
		GridPoints974();
	else if(Lmax==59)
		GridPoints1202();
	else if(Lmax==65)
		GridPoints1454();
	else if(Lmax==71)
		GridPoints1730();
	else if(Lmax==77)
		GridPoints2030();
	else if(Lmax==83)
		GridPoints2354();
	else if(Lmax==131)
		GridPoints5810();
}