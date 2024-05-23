#include<iostream>
#include <cmath>
#include <numeric/functions.h>
#include <numeric/Grid.h>
#include <common/Element.h>

using namespace std;

int main()
{
	//ifstream f("/home/tmaamaatuai/tmp/test.cube");
	//ifstream f("/home/tmaamaatuai/tmp/test256.cube");
	ifstream f("/home/tmaamaatuai/tmp/test500.cube");
	if(f.is_open()==true) cout<<"fichier créé"<<endl;
	PeriodicTable Table;
	cout<<"tp créé"<< endl;
	/*
	Atom A(Table, 64);
	Atom B;
	*/
	Grid g(f, Table);
	/*for(int i=0;i<10;i++)
	{
		cout<<"i="<<i<<endl;
		Grid lap = g.laplacian(4);
		cout<<"Sum laplacian"<<lap.sum()<<endl;
	}
	*/
	
	Grid h = g.coulomb_Grid(12, {0,0,0});
	cout<<"grid done"<<endl;
	cout<<"h done"<<endl;
	cout<<"set V done"<<endl;
	Grid prod=g*h;
	cout<<"product done"<<endl;
	double I=prod.integrate_Over_Dom();
	cout<<"sum done"<<endl;
	cout<<I<<endl;
	return 0;
}
