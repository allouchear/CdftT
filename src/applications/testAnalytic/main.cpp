#include<iostream>
#include<vector>
#include<analytic/MathFunction.h>
#include<analytic/WFX.h>

using namespace std;

int main()
{
	/*
	Factorial Table(100);
	cout<<"test pre-initialisation"<<endl;
	Binomial Bin(100, 50, Table);
	cout<<"Test post-initialisation"<<endl;

	cout<<Bin.binomial(18,7, Table)<<endl;
	cout<<Table.factorial(18)/Table.factorial(7)/Table.factorial(11)<<endl;
	cout<<endl;
	cout<<Bin.binomial(20,18, Table)<<endl;
	cout<<Table.factorial(20)/Table.factorial(18)/Table.factorial(2);
	*/
	
	ifstream f;
	f.open("h2o.wfx");
	WFX test (f);
	f.close();

	ofstream g;
	g.open("test.wfx");
	test.write_file_wfx(g);
	g.close();
	
	ifstream h;
	h.open("format.wfx");
	WFX test2;
	test2.read_file_wfx(h);
	h.close();

	ofstream i;
	i.open("test2.wfx");
	test2.write_file_wfx(i);
	i.close();

	/*cout<<"Test 1"<<endl;
	cout<<test.Electronic_Spin_Multiplicity()<<endl;
	cout<<"Test 2"<<endl;
	cout<<test.Nuclear_Charges()[1]<<endl;
	cout<<"Test 3"<<endl;
	cout<<test.Molecular_Orbital_Primitive_Coefficients()[3].Coefficients()[15]<<endl;
	cout<<"Test 4"<<endl;*/

	return 0;
}
