#include<iostream>
#include<vector>
#include<Utils/Utils.h>
#include<Basis/GTF.h>
#include<Basis/CGTF.h>
#include<Orbitals/Orbitals.h>
#include<Utils/WFX.h>
#include<Utils/FCHK.h>
#include<Utils/MOLDENGAB.h>
#include<Utils/LOG.h>
#include<Becke/Becke.h>
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

	ifstream x;
	x.open("h2o.fchk");
	FCHK fchk_h2o (x);
	x.close();

	cout<<endl;

	ifstream y;
	y.open("h2ominus.fchk");
	FCHK fchk_h2ominus (y);
	y.close();

	cout<<endl;

	ifstream z;
	z.open("h2oplus.fchk");
	FCHK fchk_h2oplus (z);
	z.close();

	cout<<endl;

	ifstream u;
	u.open("h2o.gab");
	MOLDENGAB moldengab_h2o (u);
	u.close();
	//moldengab_h2o.PrintData();

	cout<<endl;

	ifstream v;
	v.open("h2ominus.gab");
	MOLDENGAB moldengab_h2ominus (v);
	v.close();
	//moldengab_h2ominus.PrintData();

	cout<<endl;

	ifstream w;
	w.open("h2oplus.gab");
	MOLDENGAB moldengab_h2oplus (w);
	w.close();
	//moldengab_h2oplus.PrintData();

	cout<<endl;
	ifstream l1;
	l1.open("h2o.log");
	LOG log_h2o (l1);
	l1.close();
	//log_h2o.PrintData();
	cout<<endl;

	ifstream l2;
	l2.open("h2ominus.log");
	LOG log_h2ominus (l2);
	l2.close();
	//log_h2ominus.PrintData();
	cout<<endl;

	ifstream l3;
	l3.open("h2oplus.log");
	LOG log_h2oplus (l3);
	l3.close();
	//log_h2oplus.PrintData();

													//WFX
/*	Becke wh2o (wfx_h2o, Bin, Table);
	Becke wh2ominus (wfx_h2ominus, Bin, Table);
	Becke wh2oplus (wfx_h2oplus, Bin, Table);
*/
	vector<vector<double>> h2oQE;
	vector<vector<double>> h2ominusQE;
	vector<vector<double>> h2oplusQE;
/*
	h2oQE=wh2o.PartialChargesAndEnergy(1, 21 ,3);
	h2ominusQE=wh2ominus.PartialChargesAndEnergy(1, 21 ,3);
	h2oplusQE=wh2oplus.PartialChargesAndEnergy(1, 21 ,3);

	double I, A;
	I=h2oplusQE[0][0]-h2oQE[0][0];
	A=h2ominusQE[0][0]-h2oQE[0][0];
	Descriptors test(wh2o.str(), h2oQE[1], h2ominusQE[1], h2oplusQE[1], I, A);
	cout<<test<<endl;

													//FCHK
	Becke fh2o (fchk_h2o, Bin, Table);
	Becke fh2ominus (fchk_h2ominus, Bin, Table);
	Becke fh2oplus (fchk_h2oplus, Bin, Table);

	h2oQE=fh2o.PartialChargesAndEnergy(1, 21 ,3);
	h2ominusQE=fh2ominus.PartialChargesAndEnergy(1, 21 ,3);
	h2oplusQE=fh2oplus.PartialChargesAndEnergy(1, 21 ,3);

	//double I, A;
	I=h2oplusQE[0][0]-h2oQE[0][0];
	A=h2ominusQE[0][0]-h2oQE[0][0];
	Descriptors test2(fh2o.str(), h2oQE[1], h2ominusQE[1], h2oplusQE[1], I, A);
	cout<<test2<<endl;
*/													//MOLDENGAB
	Becke mgh2o (moldengab_h2o, Bin, Table);
	Becke mgh2ominus (moldengab_h2ominus, Bin, Table);
	Becke mgh2oplus (moldengab_h2oplus, Bin, Table);

	h2oQE=mgh2o.PartialChargesAndEnergy(1, 21 ,3);
	h2ominusQE=mgh2ominus.PartialChargesAndEnergy(1, 21 ,3);
	h2oplusQE=mgh2oplus.PartialChargesAndEnergy(1, 21 ,3);

	double I, A;
	I= -mgh2o.eHOMO();
	A= -mgh2o.eLUMO();
	Descriptors test2(mgh2o.str(), h2oQE[1], h2ominusQE[1], h2oplusQE[1], I, A);
	cout<<test2<<endl;

/*													//LOG
	Becke lh2o (log_h2o, Bin, Table);
	Becke lh2ominus (log_h2ominus, Bin, Table);
	Becke lh2oplus (log_h2oplus, Bin, Table);

	h2oQE=lh2o.PartialChargesAndEnergy(1, 21 ,3);
	h2ominusQE=lh2ominus.PartialChargesAndEnergy(1, 21 ,3);
	h2oplusQE=lh2oplus.PartialChargesAndEnergy(1, 21 ,3);

	double I, A;
	I=h2oplusQE[0][0]-h2oQE[0][0];
	A=h2ominusQE[0][0]-h2oQE[0][0];
	Descriptors test2(lh2o.str(), h2oQE[1], h2ominusQE[1], h2oplusQE[1], I, A);
	cout<<test2<<endl;
*/
//	Becke Becketest(fchk_test, Bin, Table);
//	Orbitals Orbtest(fchk_test, Bin, Table);
//	cout<<Orbtest<<endl;
//	Orbtest.PrintDescriptors();
/*
	Orbitals Orb(wfx_h2o, Bin, Table);
	Orb.PrintOverlap(0,0);
	Orb.PrintOverlap(1,0);
	Orb.PrintOverlap(1,1);
	Orb.PrintDescriptors();
	cout<<Orb<<endl;
*/
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
	Orbitals h2owfxtest (wfx_h2o, Bin, Table);
	cout<<h2owfxtest<<endl;
	Orbitals h2ominuswfxtest (wfx_h2ominus, Bin, Table);
	Orbitals h2opluswfxtest (wfx_h2oplus, Bin, Table);

	Orbitals h2ofchktest (fchk_h2o, Bin, Table);
	Orbitals h2ominusfchktest (fchk_h2ominus, Bin, Table);
	Orbitals h2oplusfchktest (fchk_h2oplus, Bin, Table);

	Orbitals h2omoldengabtest (moldengab_h2o, Bin, Table);
	cout<<h2omoldengabtest<<endl;
	Orbitals h2ominusmoldengabtest (moldengab_h2ominus, Bin, Table);
	Orbitals h2oplusmoldengabtest (moldengab_h2oplus, Bin, Table);

	Orbitals h2ologtest (log_h2o, Bin, Table);
	Orbitals h2ominuslogtest (log_h2ominus, Bin, Table);
	Orbitals h2opluslogtest (log_h2oplus, Bin, Table);

	string aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk,ll;

	aa="h2o_wfx_test_save.molden";
	bb="h2ominus_wfx_test_save.molden";
	cc="h2oplus_wfx_test_save.molden";
	dd="h2o_fchk_test_save.molden";
	ee="h2ominus_fchk_test_save.molden";
	ff="h2oplus_fchk_test_save.molden";
	gg="h2o_moldengab_test_save.molden";
	hh="h2ominus_moldengab_test_save.molden";
	ii="h2oplus_moldengab_test_save.molden";
	jj="h2o_log_test_save.molden";
	kk="h2ominus_log_test_save.molden";
	ll="h2oplus_log_test_save.molden";

	h2owfxtest.Save(aa);
	h2ominuswfxtest.Save(bb);
	h2opluswfxtest.Save(cc);

	h2ofchktest.Save(dd);
	h2ominusfchktest.Save(ee);
	h2oplusfchktest.Save(ff);

	h2omoldengabtest.Save(gg);
	h2ominusmoldengabtest.Save(hh);
	h2oplusmoldengabtest.Save(ii);

	h2ologtest.Save(jj);
	h2ominuslogtest.Save(kk);
	h2opluslogtest.Save(ll);

	cout<<endl;

	ifstream t1;
	t1.open(aa);
	MOLDENGAB wfx_moldengab_h2o (t1);
	t1.close();

	cout<<endl;

	ifstream t2;
	t2.open(bb);
	MOLDENGAB wfx_moldengab_h2ominus (t2);
	t2.close();

	cout<<endl;

	ifstream t3;
	t3.open(cc);
	MOLDENGAB wfx_moldengab_h2oplus (t3);
	t3.close();

	cout<<endl;

	ifstream t4;
	t4.open(dd);
	MOLDENGAB fchk_moldengab_h2o (t4);
	t4.close();

	cout<<endl;

	ifstream t5;
	t5.open(ee);
	MOLDENGAB fchk_moldengab_h2ominus (t5);
	t5.close();

	cout<<endl;

	ifstream t6;
	t6.open(ff);
	MOLDENGAB fchk_moldengab_h2oplus (t6);
	t6.close();

	cout<<endl;

	ifstream t7;
	t7.open(gg);
	MOLDENGAB moldengab_moldengab_h2o (t7);
	t7.close();

	cout<<endl;

	ifstream t8;
	t8.open(hh);
	MOLDENGAB moldengab_moldengab_h2ominus (t8);
	t8.close();

	cout<<endl;

	ifstream t9;
	t9.open(ii);
	MOLDENGAB moldengab_moldengab_h2oplus (t9);
	t9.close();

	cout<<endl;

	ifstream t10;
	t10.open(jj);
	MOLDENGAB log_moldengab_h2o (t10);
	t10.close();

	cout<<endl;

	ifstream t11;
	t11.open(kk);
	MOLDENGAB log_moldengab_h2ominus (t11);
	t11.close();

	cout<<endl;

	ifstream t12;
	t12.open(ll);
	MOLDENGAB log_moldengab_h2oplus (t12);
	t12.close();


	Orbitals h2owfxmoldengab (wfx_moldengab_h2o, Bin, Table);
	cout<<h2owfxmoldengab<<endl;
/*	Orbitals h2ominuswfxmoldengab (wfx_moldengab_h2ominus, Bin, Table);
	Orbitals h2opluswfxmoldengab (wfx_moldengab_h2oplus, Bin, Table);

	Orbitals h2ofchkmoldengab (fchk_moldengab_h2o, Bin, Table);
	Orbitals h2ominusfchkmoldengab (fchk_moldengab_h2ominus, Bin, Table);
	Orbitals h2oplusfchkmoldengab (fchk_moldengab_h2oplus, Bin, Table);

	Orbitals h2omoldengabmoldengab (moldengab_moldengab_h2o, Bin, Table);
	Orbitals h2ominusmoldengabmoldengab (moldengab_moldengab_h2ominus, Bin, Table);
	Orbitals h2oplusmoldengabmoldengab (moldengab_moldengab_h2oplus, Bin, Table);

	Orbitals h2ologmoldengab (log_moldengab_h2o, Bin, Table);
	Orbitals h2ominuslogmoldengab (log_moldengab_h2ominus, Bin, Table);
	Orbitals h2opluslogmoldengab (log_moldengab_h2oplus, Bin, Table);
*/

/*	double I, A;
	vector<vector<double>> h2oQE;
	vector<vector<double>> h2ominusQE;
	vector<vector<double>> h2oplusQE;
*/
	Becke wfxmoldengabh2o (wfx_moldengab_h2o, Bin, Table);
	Becke wfxmoldengabh2ominus (wfx_moldengab_h2ominus, Bin, Table);
	Becke wfxmoldengabh2oplus (wfx_moldengab_h2oplus, Bin, Table);

	h2oQE=wfxmoldengabh2o.PartialChargesAndEnergy(1, 21 ,3);
	h2ominusQE=wfxmoldengabh2ominus.PartialChargesAndEnergy(1, 21 ,3);
	h2oplusQE=wfxmoldengabh2oplus.PartialChargesAndEnergy(1, 21 ,3);

	I=h2oplusQE[0][0]-h2oQE[0][0];
	A=h2ominusQE[0][0]-h2oQE[0][0];
	Descriptors testwfxmoldengab(wfxmoldengabh2o.str(), h2oQE[1], h2ominusQE[1], h2oplusQE[1], I, A);
	cout<<testwfxmoldengab<<endl;


	Becke fchkmoldengabh2o (fchk_moldengab_h2o, Bin, Table);
	Becke fchkmoldengabh2ominus (fchk_moldengab_h2ominus, Bin, Table);
	Becke fchkmoldengabh2oplus (fchk_moldengab_h2oplus, Bin, Table);

	h2oQE=fchkmoldengabh2o.PartialChargesAndEnergy(1, 21 ,3);
	h2ominusQE=fchkmoldengabh2ominus.PartialChargesAndEnergy(1, 21 ,3);
	h2oplusQE=fchkmoldengabh2oplus.PartialChargesAndEnergy(1, 21 ,3);

	I=h2oplusQE[0][0]-h2oQE[0][0];
	A=h2ominusQE[0][0]-h2oQE[0][0];
	Descriptors testfchkmoldengab(fchkmoldengabh2o.str(), h2oQE[1], h2ominusQE[1], h2oplusQE[1], I, A);
	cout<<testfchkmoldengab<<endl;


	Becke moldengabmoldengabh2o (moldengab_moldengab_h2o, Bin, Table);
	Becke moldengabmoldengabh2ominus (moldengab_moldengab_h2ominus, Bin, Table);
	Becke moldengabmoldengabh2oplus (moldengab_moldengab_h2oplus, Bin, Table);

	h2oQE=moldengabmoldengabh2o.PartialChargesAndEnergy(1, 21 ,3);
	h2ominusQE=moldengabmoldengabh2ominus.PartialChargesAndEnergy(1, 21 ,3);
	h2oplusQE=moldengabmoldengabh2oplus.PartialChargesAndEnergy(1, 21 ,3);

	I=h2oplusQE[0][0]-h2oQE[0][0];
	A=h2ominusQE[0][0]-h2oQE[0][0];
	Descriptors testmoldengabmoldengab(moldengabmoldengabh2o.str(), h2oQE[1], h2ominusQE[1], h2oplusQE[1], I, A);
	cout<<testmoldengabmoldengab<<endl;


	Becke logmoldengabh2o (log_moldengab_h2o, Bin, Table);
	Becke logmoldengabh2ominus (log_moldengab_h2ominus, Bin, Table);
	Becke logmoldengabh2oplus (log_moldengab_h2oplus, Bin, Table);

	h2oQE=logmoldengabh2o.PartialChargesAndEnergy(1, 21 ,3);
	h2ominusQE=logmoldengabh2ominus.PartialChargesAndEnergy(1, 21 ,3);
	h2oplusQE=logmoldengabh2oplus.PartialChargesAndEnergy(1, 21 ,3);

	I=h2oplusQE[0][0]-h2oQE[0][0];
	A=h2ominusQE[0][0]-h2oQE[0][0];
	Descriptors testlogmoldengab(logmoldengabh2o.str(), h2oQE[1], h2ominusQE[1], h2oplusQE[1], I, A);
	cout<<testlogmoldengab<<endl;

	cout<<"Temps d'execution : "<<time.get()<<" ms"<<endl;

	return 0;
}
