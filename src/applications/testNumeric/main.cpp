#include<iostream>
#include <cmath>
#include <numeric/functions.h>
#include <numeric/Grid.h>
#include <numeric/GridCP.h>
#include <common/Element.h>
#include <sys/time.h>
#include "Timer.h"
#include <common/Descriptors.h>
#include <analytic/Becke/Becke.h>

using namespace std;

int main()
{
	PeriodicTable Table;
	Factorial fact(100);
	Binomial bino (100, fact);
	Timer timer;
	ifstream h2o80("./h2o80.gcube");
	ifstream h2o256("./h2o_256.gcube");
	ifstream h2o400("./h2o_400.gcube");

	ifstream wf;
	wf.open("h2o.wfx");
	WFX wfx_h2o(wf);
	wf.close();
	Becke beck_h2o (wfx_h2o, bino, Table);
	vector<double> qw = beck_h2o.PartialChargeAndEnergy();
	for(size_t i=0;i<qw.size();i++)
	{
		cout<<"wfx Qk = "<<qw[i]<<endl;
	}

	Grid grid80(h2o80, Table);
	Becke B80(grid80);
	vector<double> Q80=B80.PartialChargeAndEnergy(grid80);
	cout<<"grid80 80 "<<endl;
	for(size_t i=0;i<Q80.size();i++)
	{
		cout<<"Qk = "<<Q80[i]<<endl;
	}

	Grid grid256(h2o256, Table);
	Becke B256(grid256);
	vector<double> Q256=B256.PartialChargeAndEnergy(grid256);
	cout<<"grid256 256 "<<endl;
	for(size_t i=0;i<Q80.size();i++)
	{
		cout<<"Qk = "<<Q256[i]<<endl;
	}

	Grid grid4(h2o400, Table);
	Becke B4(grid4);
	vector<double> Q4=B4.PartialChargeAndEnergy(grid4);
	cout<<"grid4 400 "<<endl;
	for(size_t i=0;i<Q80.size();i++)
	{
		cout<<"Qk = "<<Q4[i]<<endl;
	}


	/*
	GridCP gridcp;
	cout<<"Begin Sign"<<endl;
	timer.init();
	gridcp.buildBasins(grid,0);
	gridcp.computeIntegrals(grid);
	gridcp.printCriticalPoints();
	vector<double> charges=gridcp.computeAIMCharges(grid);
	cout<<"Time in ms "<<timer.get()<<endl;
	return 0;
	*/
}
