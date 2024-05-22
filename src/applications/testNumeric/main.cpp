#include<iostream>
#include <cmath>
#include <numeric/functions.h>
#include <numeric/Grid.h>
#include <common/Element.h>

using namespace std;

int main()
{
	ifstream f("/home/tmaamaatuai/tmp/test.cube");
	if(f.is_open()==true) cout<<"fichier créé"<<endl;
	PeriodicTable Table;
	cout<<"tp créé"<< endl;
	Atom A(Table, 64);
	Atom B;
	Grid g(f, Table);
	vector<double> par={1};
	Grid h;
	cout<<"grid done"<<endl;
	h.resize_zeros(1,2,3);
	cout<<"h done"<<endl;
	h.set_V_Func(par, &laplacian, g);
	cout<<"set V done"<<endl;
	Grid prod=g*h*g;
	cout<<"product done"<<endl;
	double I=prod.integrate_over_dom();
	cout<<"sum done"<<endl;
	cout<<I<<endl;
	return 0;
}
