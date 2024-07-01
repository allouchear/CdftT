#include<iostream>
#include <cmath>
#include <numeric/functions.h>
#include <numeric/Grid.h>
#include <numeric/GridCP.h>
#include <common/Element.h>
#include <sys/time.h>
#include "Timer.h"
#include <common/Descriptors.h>

using namespace std;

int main()
{
	PeriodicTable Table;
	Timer timer;

	ifstream h2o("./h2o80.gcube");
	Grid grid(h2o, Table);

	GridCP gridcp;
	cout<<"Begin Sign"<<endl;
	timer.init();
	gridcp.buildBasins(grid,0);
	gridcp.computeIntegrals(grid);
	gridcp.printCriticalPoints();
	vector<double> charges=gridcp.computeAIMCharges(grid);
	cout<<"Time in ms "<<timer.get()<<endl;
	return 0;
}
