#include<iostream>
#include<analytic/Becke/GridPoints.h>

using namespace std;

GridPoints::GridPoints()
{
	_Lebedev_Npts=0;
	_Lebedev_Lmax=0;
	_Lebedev_L2max=0;
	_LebedevGridPoints=vector<vector<double>> ();
}

GridPoints::GridPoints(int i)
{
	if(i==6)
		GridPoints6();
	else if(i==14)
		GridPoints14();
	else if(i==26)
		GridPoints26();
	else if(i==38)
		GridPoints38();
	else if(i==50)
		GridPoints50();
	else if(i==74)
		GridPoints74();
	else if(i==86)
		GridPoints86();
	else if(i==110)
		GridPoints110();
	else if(i==146)
		GridPoints146();
	else if(i==170)
		GridPoints170();
	else if(i==194)
		GridPoints194();
	else if(i==230)
		GridPoints230();
	else if(i==266)
		GridPoints266();
	else if(i==302)
		GridPoints302();
	else if(i==350)
		GridPoints350();
	else if(i==434)
		GridPoints434();
	else if(i==590)
		GridPoints590();
	else if(i==770)
		GridPoints770();
	else if(i==974)
		GridPoints974();
	else if(i==1202)
		GridPoints1202();
	else if(i==1454)
		GridPoints1454();
	else if(i==1730)
		GridPoints1730();
	else if(i==2030)
		GridPoints2030();
	else if(i==2354)
		GridPoints2354();
	else if(i==5810)
		GridPoints5810();
}