#include<iostream>
#include<vector>
#include<analytic/MathFunction.h>
#include<analytic/GTF.h>
#include<analytic/CGTF.h>
#include<analytic/LCAO.h>
#include<analytic/WFX.h>
#include"Timer.h"
#include<ctime>

using namespace std;

int main()
{
	Timer time;

	Factorial Table(100);
	cout<<"Test pre-initialisation"<<endl;
	Binomial Bin(100, Table);
	cout<<"Test post-initialisation"<<endl;

	cout<<Bin.binomial(18,7)<<endl;
	cout<<Table.factorial(18)/Table.factorial(7)/Table.factorial(11)<<endl;
	cout<<endl;
	cout<<Bin.binomial(20,18)<<endl;
	cout<<Table.factorial(20)/Table.factorial(18)/Table.factorial(2)<<endl;
	cout<<endl;
	cout<<Table.double_factorial(10)<<endl<<endl;
	
	
	ifstream f;
	f.open("h2o.wfx");
	WFX wfx_test (f);
	f.close();

	ofstream g;
	g.open("test.wfx");
	wfx_test.write_file_wfx(g);
	g.close();
	
	int i,j;

	GTF gtf_test;
	vector<vector<CGTF>> vcgtf_test(wfx_test.Number_of_Occupied_Molecular_Orbital(), vector<CGTF> (wfx_test.Number_of_Primitives()));
	vector<LCAO> vlcao_test(wfx_test.Number_of_Occupied_Molecular_Orbital());

	vector<vector<double>> Coord_test(wfx_test.Number_of_Nuclei(), vector<double> (0));
	for(i=0; i<wfx_test.Number_of_Nuclei(); i++)
		for(j=i*3; j<3*(1+i); j++)
			Coord_test[i].push_back(wfx_test.Nuclear_Cartesian_Coordinates()[j]);

	for(i=0; i<wfx_test.Number_of_Occupied_Molecular_Orbital(); i++)
	{
		gtf_test=GTF();

		for(j=0; j<wfx_test.Number_of_Primitives(); j++)
		{			
			gtf_test.push_back(wfx_test.Primitive_Exponents()[j], 1.0, Coord_test[wfx_test.Primitive_Centers()[j]-1], setLxyz(wfx_test.Primitive_Types()[j]), Bin);
			vcgtf_test[i][j].push_back(gtf_test);		
			vlcao_test[i].push_back(vcgtf_test[i][j], wfx_test.Molecular_Orbital_Primitive_Coefficients()[i].Coefficients()[j]);
		}
	}

	//for(i=0; i<wfx_test.Number_of_Occupied_Molecular_Orbital(); i++)
	//	cout<<vlcao_test[i]<<endl;

	cout<<"Test Post Overlap ok"<<endl;

	for(i=0; i<wfx_test.Number_of_Occupied_Molecular_Orbital(); i++)
		for(j=0; j<wfx_test.Number_of_Occupied_Molecular_Orbital(); j++)
			cout<<"OverlapLCAO <"<<i<<"|"<<j<<"> = "<<vlcao_test[i].overlapLCAO(vlcao_test[j])<<endl;

	cout<<"Test Overlap ok"<<endl;
	/*
	ifstream h;
	h.open("format.wfx");
	WFX test2;
	test2.read_file_wfx(h);
	h.close();

	ofstream i;
	i.open("test2.wfx");
	test2.write_file_wfx(i);
	i.close();

	cout<<"Test 1"<<endl;
	cout<<test.Electronic_Spin_Multiplicity()<<endl;
	cout<<"Test 2"<<endl;
	cout<<test.Nuclear_Charges()[1]<<endl;
	cout<<"Test 3"<<endl;
	cout<<test.Molecular_Orbital_Primitive_Coefficients()[3].Coefficients()[15]<<endl;
	cout<<"Test 4"<<endl;*/

	cout<<"Temps d'execution : "<<time.get()<<" ms"<<endl;

	return 0;
}