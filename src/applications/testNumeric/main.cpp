#include<iostream>
#include <cmath>
#include <numeric/functions.h>
#include <numeric/Grid.h>
#include <common/Element.h>
#include <sys/time.h>
#include "Timer.h"

using namespace std;

int main()
{
	Timer timer;
	//ifstream f("/home/tmaamaatuai/tmp/propane/Propane256.cube");
	ifstream f("/home/tmaamaatuai/tmp/CH4/ch4256.cube");
	//ifstream f("/home/tmaamaatuai/tmp/h2o/h2odens.cube");
	//ifstream f("/home/tmaamaatuai/tmp/test256.cube");
	//ifstream f("/home/tmaamaatuai/tmp/test500.cube");
	if(f.is_open()) cout<<"fichier créé"<<endl;
	PeriodicTable Table;
	cout<<"tp créé"<< endl;
	/*
	Atom A(Table, 64);
	Atom B;
	*/
	
	Grid g(f, Table);
	cout<<g.dom().dx()<<endl;
	cout<<g.dom().O()[0]<<endl;
	cout<<g.str().atoms()[0].coordinates()[0]<<endl;
	cout<<g.V()[127][92][100][0]<<endl;
	cout<<g.V()[128][92][100][0]<<endl;
	
	
	cout<<"Begin gradient"<<endl;
	timer.init();
	Grid lap = g.gradient(4);
	cout<<"Time in ms "<<timer.get()<<endl;
	
	/*
	cout<<"Begin AIM"<<endl;
	timer.init();
	Grid AIM=g.aim_On_Grid_Density();
	cout<<"Time in ms "<<timer.get()<<endl;
	*/
	ofstream out("/home/tmaamaatuai/tmp/gradsave.cube");
	cout<<"Begin write"<<endl;
	timer.init();
	lap.save(out);
	cout<<"Time in ms "<<timer.get()<<endl;
	
	
	/*
	for(int i=0;i<10;i++)
	{
		cout<<"i="<<i<<endl;
		Grid lap = g.gradient(4);
		cout<<"Sum laplacian= "<<lap.sum()<<endl;
	}*/
	/*
	for(int i=0;i<10;i++)
	{
		cout<<"i="<<i<<endl;
		Grid lap = g.gradient(4);
		cout<<"Sum gradient= "<<lap.sum()<<endl;
	}
	*/
	/*
	timer.init();
	cout<<"Begin finer grid function"<<endl;
	Grid G=g.finer_Grid();
	cout<<"Time in ms "<<timer.get()<<endl;
	*/

/*
	timer.init();
	cout<<"Begin coarser grid function"<<endl;
	Grid G=g.coarser_Grid();
	cout<<"Time in ms "<<timer.get()<<endl;

	timer.init();
	cout<<"Begin Coulomb grid function"<<endl;
	Grid h = G.coulomb_Grid(12, {0,0,0});
	cout<<"Time in ms "<<timer.get()<<endl;

	cout<<"grid done"<<endl;
	cout<<"h done"<<endl;
	cout<<"set V done"<<endl;
	timer.init();
	cout<<"Begin product grid function"<<endl;
	Grid prod=g*h;
	cout<<"Time in ms "<<timer.get()<<endl;

	cout<<"product done"<<endl;
	double I=prod.integrate_Over_Dom();
	cout<<"sum done"<<endl;
	cout<<I<<endl;*/
	return 0;
}
