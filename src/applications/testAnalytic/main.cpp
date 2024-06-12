#include<iostream>
#include<vector>
#include<analytic/Utils/MathFunction.h>
#include<analytic/Basis/GTF.h>
#include<analytic/Basis/CGTF.h>
#include<analytic/Orbitals/LCAO.h>
#include<analytic/Orbitals/Orbitals.h>
#include<analytic/Utils/WFX.h>
#include<analytic/Becke/Becke.h>
#include"Timer.h"
#include<ctime>

using namespace std;

int main()
{
	Timer time;

	Factorial Fact(100);
	//cout<<"Test pre-initialisation"<<endl;
	Binomial Bin(100, Fact);
	//cout<<"Test post-initialisation"<<endl;
	PeriodicTable Table;
/*
	cout<<Bin.binomial(18,7)<<endl;
	cout<<Table.factorial(18)/Table.factorial(7)/Table.factorial(11)<<endl;
	cout<<endl;
	cout<<Bin.binomial(20,18)<<endl;
	cout<<Table.factorial(20)/Table.factorial(18)/Table.factorial(2)<<endl;
	cout<<endl;
	cout<<Table.double_factorial(10)<<endl<<endl;
*/
	
	ifstream f;
	f.open("h2o.wfx");
	WFX wfx_h2o (f);
	f.close();

	cout<<endl;

	ifstream g;
	g.open("h2ominus.wfx");
	WFX wfx_h2ominus (g);
	g.close();

	cout<<endl;

	ifstream h;
	h.open("h2oplus.wfx");
	WFX wfx_h2oplus (h);
	h.close();

	cout<<endl;

	/*ofstream g;
	g.open("test.wfx");
	wfx_test.write_file_wfx(g);
	g.close();*/
	
	Becke h2o (wfx_h2o, Bin, Table);
	Becke h2ominus (wfx_h2ominus, Bin, Table);
	Becke h2oplus (wfx_h2oplus, Bin, Table);

	vector<vector<double>> all (3, vector<double> ());

	all[0]=h2o.PartialChargeAndEnergy(3, 23, 1);
	all[1]=h2ominus.PartialChargeAndEnergy(3, 23, 1);
	all[2]=h2oplus.PartialChargeAndEnergy(3, 23, 1);

	Descriptors test(wfx_h2o, Table);
	test.set_mu_fk_data(all);
	test.compute_all();
	cout<<test<<endl;

	//Orbitals Orbtest(wfx_h2o, Bin, Table);



	//Orbtest.PrintDescriptors();

/*
	GTF A,B,C;
	vector<int> l1(3), l2(3);
	l1 = {0, 0 ,2};
	l2 = {0, 1, 1};
	vector<double> coord1(3), coord2(3);
	coord1 = {-6.778167794095e-34,  5.422534235276e-33,  2.199641208994e-01};
	coord2 = {0.000000000000e+00,  1.419184319548e+00, -8.798564835975e-01};
	A=B=GTF(8.588500000000e+03, 1.0, coord1, l1, Bin);
	C=GTF(1.297230000000e+03, 1.0, coord2, l2, Bin);

	A=B=GTF(8.588500000000e+00, 1.0, coord1, l1, Bin);
	C=GTF(1.297230000000e+00, 1.0, coord2, l2, Bin);
*/
/*
	LCAO A,B,C;
	A=Orbtest.vlcao()[0];
	B=Orbtest.vlcao()[1];
	C=Orbtest.vlcao()[2];
	cout<<setprecision(10);

	double OLAA, OLCC, OLAB, OLAC, OLBC;
	OLAA=Becketest.overlapLCAO(A,A);
	OLCC=Becketest.overlapLCAO(C,C);
	OLAB=Becketest.overlapLCAO(A,B);
	OLAC=Becketest.overlapLCAO(A,C);
	OLBC=Becketest.overlapLCAO(B,C);

	cout<<endl<<right<<setw(30)<<"*** TEST BECKE ***"<<setw(30)<<"*** ANALYTIQUE ***"<<setw(30)<<"*** DIFF ***"<<endl;

	cout<<left<<setw(10)<<"<A|A> = "<<right<<setw(20)<<OLAA<<setw(30)<<A.overlapLCAO(A)<<setw(30)<<OLAA-A.overlapLCAO(A)<<endl;
	cout<<left<<setw(10)<<"<C|C> = "<<right<<setw(20)<<OLCC<<setw(30)<<C.overlapLCAO(C)<<setw(30)<<OLCC-C.overlapLCAO(C)<<endl;
	cout<<left<<setw(10)<<"<A|B> = "<<right<<setw(20)<<OLAB<<setw(30)<<A.overlapLCAO(B)<<setw(30)<<OLAB-A.overlapLCAO(B)<<endl;
	cout<<left<<setw(10)<<"<A|C> = "<<right<<setw(20)<<OLAC<<setw(30)<<A.overlapLCAO(C)<<setw(30)<<OLAC-A.overlapLCAO(C)<<endl;
	cout<<left<<setw(10)<<"<B|C> = "<<right<<setw(20)<<OLBC<<setw(30)<<B.overlapLCAO(C)<<setw(30)<<OLBC-B.overlapLCAO(C)<<endl;
*/
/*
	cout<<endl<<"***** TEST BECKE *****"<<endl;
	cout<<"<C|C> = "<<Becketest.overlapGTF(C,C)<<endl;
	cout<<"<A|B> = "<<Becketest.overlapGTF(A,B)<<endl;
	cout<<"<A|C> = "<<Becketest.overlapGTF(A,C)<<endl;
	cout<<"<B|C> = "<<Becketest.overlapGTF(B,C)<<endl;

	cout<<endl<<"***** TEST CHECK *****"<<endl;
	cout<<"<A|A> = "<<A.overlapGTF(A)<<endl;
	cout<<"<C|C> = "<<C.overlapGTF(C)<<endl;
	cout<<"<A|B> = "<<A.overlapGTF(B)<<endl;
	cout<<"<A|C> = "<<A.overlapGTF(C)<<endl;
	cout<<"<B|C> = "<<B.overlapGTF(C)<<endl;
*/

	//Orbitals Orb_test(wfx_test, Bin);

/*
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

	for(i=0; i<wfx_test.Number_of_Occupied_Molecular_Orbital(); i++)
		cout<<vlcao_test[i]<<endl;

	cout<<"Test Post Overlap ok"<<endl;

	for(i=0; i<wfx_test.Number_of_Occupied_Molecular_Orbital(); i++)
		for(j=0; j<wfx_test.Number_of_Occupied_Molecular_Orbital(); j++)
			cout<<"OverlapLCAO <"<<i<<"|"<<j<<"> = "<<vlcao_test[i].overlapLCAO(vlcao_test[j])<<endl;

	vector<vector<double>> vtest (wfx_test.Number_of_Occupied_Molecular_Orbital(), vector<double> (wfx_test.Number_of_Occupied_Molecular_Orbital()));

	for(i=0; i<wfx_test.Number_of_Occupied_Molecular_Orbital(); i++)
		for(j=0; j<wfx_test.Number_of_Occupied_Molecular_Orbital(); j++)
			vtest[i][j]=vlcao_test[i].overlapLCAO(vlcao_test[j]);

	cout<<"Test Overlap ok"<<endl;
*/
	//for(int i=0; i<Orb_test.NumberOfFunctions(); i++)
	//	cout<<Orb_test.lcao(i)<<endl;
/*
	for(int i=0; i<Orb_test.NumberOfFunctions(); i++)
		for(int j=0; j<Orb_test.NumberOfFunctions(); j++)
			Orb_test.PrintOverlapGTF(i,j);
*/

/*
	Orb_test.HOMO_LUMO(0,191);
	Orb_test.get_f();
	cout<<Orb_test<<endl;
*/
/*
	Orb_test.HOMO_LUMO();
	Orb_test.get_f();
	cout<<Orb_test<<endl;
*/
/*
	for(int i=0; i<10; i++)
	{
		for(int j=0; j<3; j++)
			cout<<Orb_test.get_f(i)[j]<<endl;
		cout<<endl;
	}
*/
	//cout<<"mu ="<<(Orb_test.eHOMO()+Orb_test.eLUMO())/2<<endl;

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
	cout<<"Test 4"<<endl;
*/
	cout<<"Temps d'execution : "<<time.get()<<" ms"<<endl;

	return 0;
}