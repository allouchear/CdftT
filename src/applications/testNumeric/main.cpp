#include<iostream>
#include <numeric/Grid.h>
#include <common/Element.h>

using namespace std;

int main()
{
	ifstream f("/home/tmaamaatuai/tmp/test.cube");
	if(f.is_open()==true) cout<<"fichier créé"<<endl;
	Element e;
	cout<<"element créé"<<endl;
	PeriodicTable Table;
	cout<<"tp créé"<< endl;
	Atom A(Table, 64);
	Atom B;
	cout<<"atom créé"<<endl;
	Domain d(f);
	cout<<d.Natoms()<<endl;
	Structure S(f, d.Natoms(), Table);
	cout<<"Structure créé"<<endl;
	Grid g(f, Table);
	for(int i=0; i<g.dom().N1();i++)
	{	
		for(int j=0; j<g.dom().N2();i++)
		{
			for(int k=0; k<g.dom().N3();i++)
			{
				for(int l=0; l<g.dom().Nval();i++)
				{
				cout<<g.V()[i][j][k][l]<< " ";
				}
			}
		}
	}
	return 0;
}
