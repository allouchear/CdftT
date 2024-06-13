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
	Timer timer;
	//ifstream f("ch4256.cube");
	//ifstream f("propane400.cube");
	//ifstream f("propaneGrad400.cube");
	//ifstream f("h2o80.cube");
	//ifstream f("h2o80.gcube");
	//ifstream f("lap.gcube");
	//ifstream f("h2o256.cube");
	//ifstream f("h2oGrad256.cube");
	//ifstream f("h2oGrad80.cube");
	//ifstream f("propaneGrad256.cube");
	//ifstream f("h2o192.cube");
	//ifstream f("/home/tmaamaatuai/tmp/CH4/ch4256.cube");
	//ifstream f("/home/tmaamaatuai/tmp/h2o/h2odens.cube");
	//ifstream f("/home/tmaamaatuai/tmp/test256.cube");
	//ifstream f("/home/tmaamaatuai/tmp/test500.cube");
	
	ifstream f("h2o_0_Grad80.cube");
	ifstream e("h2o_m_Grad80.cube");
	ifstream d("h2o_p_Grad80.cube");
	
	//ifstream f("Ala_0_256.cube");
	//ifstream e("Ala_m_256.cube");
	//ifstream d("Ala_p_256.cube");
	
	if(f.is_open()==true) cout<<"fichier créé"<<endl;
	PeriodicTable Table;
	cout<<"tp créé"<< endl;
	/*
	Atom A(Table, 64);
	Atom B;
	*/
	
	Grid AIM0(f, Table);
	Grid AIMM(e, Table);
	Grid AIMP(d, Table);
	
	/*
	cout<<g.dom().dx()<<endl;
	cout<<g.dom().O()[0]<<endl;
	cout<<g.str().atoms()[0].coordinates()[0]<<endl;
	//cout<<g.V()[127][92][100][0]<<endl;
	//cout<<g.V()[128][92][100][0]<<endl;
	Grid lap=g.laplacian(2);
	ofstream S("savelap.cube");
	lap.save(S);
	*/
	
	/*
	GridCP gridcp0;
	cout<<"Begin AIM on grid"<<endl;
	timer.init();
	gridcp0.buildAttractors(AIM0,0);
	gridcp0.printCriticalPoints();
	vector<double> charges0=gridcp0.computeAIMCharges(AIM0);
	cout<<"Time in ms "<<timer.get()<<endl;
	
	GridCP gridcpm;
	cout<<"Begin AIM near grid without refinement"<<endl;
	
	timer.init();
	gridcpm.buildAttractors(AIMM,0);
	gridcpm.printCriticalPoints();
	vector<double> chargesM=gridcpm.computeAIMCharges(AIMM);
	cout<<"Time in ms "<<timer.get()<<endl;
	
	GridCP gridcpp;
	cout<<"Begin AIM near grid with refinement"<<endl;
	timer.init();
	gridcpp.buildAttractors(AIMP,0);
	//gridcp.printCriticalPoints();
	vector<double> chargesP=gridcpp.computeAIMCharges(AIMP);
	cout<<"Time in ms "<<timer.get()<<endl;
	
	for(int i=0; i<int(chargesP.size());i++)
	{
		cout<<"Q0 "<<charges0[i]<<endl;
		cout<<"Qm "<<chargesM[i]<<endl;
		cout<<"Qp "<<chargesP[i]<<endl;
	}
	
	Descriptors D(gridcp0.str(), charges0, chargesM, chargesP, 0.46, 0.11);
	cout<<D;
	
	vector<double> Q0={-0.86200, 0.43100,0.43100};
	vector<double> QM={-1.06949,0.03474,0.03474};
	vector<double> QP={-0.018,0.50930,0.50930};
	Descriptors E(gridcp0.str(), Q0, QM, QP, 0.46, 0.11);
	cout<<E;
	*/
	cout<<"Begin Descriptor calculation"<<endl;
	timer.init();
	Descriptors D(AIM0,AIMM,AIMP,0.309,0.057,2);
	cout<<D;
	cout<<"Time in ms "<<timer.get()<<endl;
	
	/*
	Descriptors D(f,e,d,0.46,0.11,2);
	cout<<D;
	*/
	
	/*GridCP gridcp;
	cout<<"Begin AIM on grid"<<endl;
	timer.init();
	gridcp.buildAttractors(g,0);
	gridcp.printCriticalPoints();
	gridcp.computeAIMCharges(g);
	cout<<"Time in ms "<<timer.get()<<endl;
	cout<<"Begin AIM near grid without refinement"<<endl;
	timer.init();
	gridcp.buildAttractors(g,1);
	gridcp.printCriticalPoints();
	vector<double> v=gridcp.computeAIMCharges(g);
	cout<<"Time in ms "<<timer.get()<<endl;
	cout<<"Begin AIM near grid with refinement"<<endl;
	timer.init();
	gridcp.buildAttractors(g,2);
	gridcp.printCriticalPoints();
	vector<double> u=gridcp.computeAIMCharges(g);
	cout<<"Time in ms "<<timer.get()<<endl;
	*/
	return 0;
}
