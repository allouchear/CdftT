#include<iostream>
#include<iomanip>
#include<analytic/WFX.h>

using namespace std;

//****************************************//
//************ MOPC class ****************//
//****************************************//

MOPC::MOPC()
{
	_MO_Number=0;
	_Coefficients.resize(0);
}

void MOPC::push_back(int a, vector<double> b)
{
	_MO_Number=a;
	_Coefficients=b;
}


//****************************************//
//************ NCEG class ****************//
//****************************************//

NCEG::NCEG()
{
	_symbol="None";
	_gradient.resize(0);
}

void NCEG::push_back(string a, vector<double> b)
{
	_symbol=a;
	_gradient=b;
}


//****************************************//
//*********** AEDF class *****************//
//****************************************//

AEDF::AEDF()
{
	_Number_of_EDF_Primitives=0;
	_EDF_Primitives_Centers.resize(0);
	_EDF_Primitives_Types.resize(0);
	_EDF_Primitives_Exponents.resize(0);
	_EDF_Primitives_Coefficients.resize(0);
}

void AEDF::push_back(int a, vector<int> b, vector<int> c, vector<double> d, vector<double> e)
{
	_Number_of_EDF_Primitives=a;
	_EDF_Primitives_Centers=b;
	_EDF_Primitives_Types=c;
	_EDF_Primitives_Exponents=d;
	_EDF_Primitives_Coefficients=e;
}


//****************************************//
//*********** WFX class ******************//
//****************************************//

WFX::WFX()
{
	_Title="None";														
	_Keywords="None";													
	_Number_of_Nuclei=0;												
	_Number_of_Primitives=0;	
	_Number_of_Occupied_Molecular_Orbital=0;	
	_Number_of_Perturbations=0;									
	_Nuclear_Names.resize(0);	
	_Atomic_Number.resize(0);
	_Nuclear_Charges.resize(0);	
	_Nuclear_Cartesian_Coordinates.resize(0);
	_Net_Charge=0.0;													
	_Number_of_Electrons=0;											
	_Number_of_Alpha_Electrons=0;										
	_Number_of_Beta_Electrons=0;										
	_Electronic_Spin_Multiplicity=0;
	_Model="None";							
	_Primitive_Centers.resize(0);										
	_Primitive_Types.resize(0);											
	_Primitive_Exponents.resize(0);									
	_Molecular_Orbital_Occupation_Numbers.resize(0);				
	_Molecular_Orbital_Energies.resize(0);
	_Molecular_Orbital_Spin_Types.resize(0);					
	_Molecular_Orbital_Primitive_Coefficients.resize(0);
	_Energy=0.0;
	_Virial_Ratio=0.0;
	_Nuclear_Cartesian_Energy_Gradients.resize(0);
	_Nuclear_Virial_of_Energy_Gradient_Based_Forces_on_Nuclei=0.0;
	_Full_Virial_Ratio=0.0;
	_Number_of_Core_Electrons=0;
	_Additionnal_Electron_Density_Function=AEDF();
}

WFX::WFX(ifstream& file)
{
	read_file_wfx(file);
}

vector<int> WFX::read_one_block_int(ifstream& f, string b, bool r, int nint)
{
	int data;
	vector<int> block_of_data(0);
	long int pos;
	int i,n;
	pos=LocaliseBlock(f, n, b);
	if(pos==-1)
	{
		cout<<b+" : data not found"<<endl;
		if(r)
		{
			cout<<"Data required, please check your file"<<endl;
			exit(1);
		}
		else
			return vector<int> (0);
	}
	f.seekg(pos);
	for(i=0; i<nint; i++)
	{
		f>>data;
		block_of_data.push_back(data);
	}
	return block_of_data;
}

vector<double> WFX::read_one_block_real(ifstream& f, string b, bool r)
{
	double data;
	vector<double> block_of_data(0);
	long int pos;
	int i,n;
	pos=LocaliseBlock(f, n, b);
	if(pos==-1)
	{
		cout<<b+" : data not found"<<endl;
		if(r)
		{
			cout<<"Data required, please check your file"<<endl;
			exit(1);
		}
		else
			return vector<double> (0);
	}
	f.seekg(pos);
	for(i=0; i<n; i++)
	{
		f>>data;
		block_of_data.push_back(data);
	}

	return block_of_data;
}

vector<double> WFX::read_one_block_real(ifstream& f, string b, bool r, int ndouble)
{
	double data;
	vector<double> block_of_data(0);
	long int pos;
	int i,n;
	pos=LocaliseBlock(f, n, b);
	if(pos==-1)
	{
		cout<<b+" : data not found"<<endl;
		if(r)
		{
			cout<<"Data required, please check your file"<<endl;
			exit(1);
		}
		else
			return vector<double> (0);
	}
	f.seekg(pos);
	for(i=0; i<ndouble; i++)
	{
		f>>data;
		block_of_data.push_back(data);
	}

	return block_of_data;
}

vector<string> WFX::read_one_block_string(ifstream& f, string b, bool r)
{
	string data;
	vector<string> block_of_data(0);
	long int pos;
	int i,n;
	pos=LocaliseBlock(f, n, b);
	if(pos==-1)
	{
		cout<<b+" : data not found"<<endl;
		if(r)
		{
			cout<<"Data required, please check your file"<<endl;
			exit(1);
		}
		else
			return vector<string> (0);
	}
	f.seekg(pos);
	for(i=0; i<n; i++)
	{
		getline(f,data);
		block_of_data.push_back(data);
	}

	return block_of_data;
}

int WFX::read_int(ifstream& f, string b, bool r)
{
	int data,i,n;
	long int pos;
	pos=LocaliseBlock(f, n, b);
	if(pos==-1)
	{
		cout<<b+" : data not found"<<endl;
		if(r)
		{
			cout<<"Data required, please check your file"<<endl;
			exit(1);
		}
		else
			return 0;
	}
	f.seekg(pos);
	for(i=0; i<n; i++)
		f>>data;

	return data;
}

double WFX::read_real(ifstream& f, string b, bool r)
{
	double data;
	int i,n;
	long int pos;
	pos=LocaliseBlock(f, n, b);
	if(pos==-1)
	{
		cout<<b+" : data not found"<<endl;
		if(r)
		{
			cout<<"Data required, please check your file"<<endl;
			exit(1);
		}
		else
			return 0.0;
	}
	f.seekg(pos);
	for(i=0; i<n; i++)
		f>>data;

	return data;
}

string WFX::read_string(ifstream& f, string b, bool r)
{
	string data;
	int i,n;
	long int pos;
	pos=LocaliseBlock(f, n, b);
	if(pos==-1)
	{
		cout<<b+" : data not found"<<endl;
		if(r)
		{
			cout<<"Data required, please check your file"<<endl;
			exit(1);
		}
		else
			return "None";
	}
	f.seekg(pos);
	for(i=0; i<n; i++)
		getline(f,data);

	return data;
}

vector<MOPC> WFX::read_MOPC_block(ifstream& f, string b, bool r)
{
	int i,j,n,nMO,ndat;
	double c;
	vector<double> cdat(0);
	MOPC data;
	vector<MOPC> block_of_data(0);
	string MO="MO Number";
	string p;
	string test;
	long int pos;
	pos=LocaliseMO(f,n,nMO,b);

	if(pos==-1)
	{
		cout<<b+" : data not found"<<endl;
		if(r)
		{
			cout<<"Data required, please check your file"<<endl;
			exit(1);
		}
		else
			return vector<MOPC> (0);
	}

	f.clear();
	f.seekg(pos);
	n=_Number_of_Primitives;

	for(i=0; i<nMO/2; i++)
	{
		getline(f,p);
		f>>ndat;
		getline(f,p);
		getline(f,p);

		for(j=0; j<n; j++)
		{
			f>>c;
			cdat.push_back(c);
		}
		data.push_back(ndat,cdat);
		block_of_data.push_back(data);
		cdat.resize(0);
		getline(f,p);
	}

	return block_of_data;
}

vector<NCEG> WFX::read_NCEG_block(ifstream& f, string b, bool r)
{
	int i,j,n;
	long int pos;

	string s;
	double g;
	vector<double> gdat(0);
	NCEG data;
	vector<NCEG> block_of_data(0);

	pos=LocaliseBlock(f,n,b);
	if(pos==-1)
	{
		cout<<b+" : data not found"<<endl;
		if(r)
		{
			cout<<"Data required, please check your file"<<endl;
			exit(1);
		}
		else
			return vector<NCEG> (0);
	}

	f.seekg(pos);
	for(i=0; i<n/2; i++)
	{
		f>>s;
		for(j=0; j<3; j++)
		{
			f>>g;
			gdat.push_back(g);
		}
		data.push_back(s,gdat);
		block_of_data.push_back(data);
	}

	return block_of_data;
}

AEDF WFX::read_AEDF_block(ifstream& f, string b, bool r)
{
	string snp="Number of EDF Primitives";
	string spc="EDF Primitive Centers";
	string spt="EDF Primitive Types";
	string spe="EDF Primitive Exponents";
	string spcoef="EDF Primitive Coefficients";
	int n=read_int(f, snp, r);
	vector<int> pc=read_one_block_int(f, spc, r, n);
	vector<int> pt=read_one_block_int(f, spt, r, n);
	vector<double> pe=read_one_block_real(f, spe, r);
	vector<double> pcoef=read_one_block_real(f, spcoef, r);

	AEDF block_of_data;
	block_of_data.push_back(n, pc, pt, pe, pcoef);

	return block_of_data;
}

long int LocaliseBlock(ifstream& f, int& n, string b)
{
	n=0;
	long int position;
	f.clear();
	f.seekg(0,f.beg);
	string test;
	bool ok=false;
	while(!f.eof())
	{	
		getline(f, test);
		if(test.find(b)!=string::npos)
		{
			ok=true;
			position=f.tellg();
			break;
		}
	}

	if(!ok) 
		return -1;

	while(!f.eof())
	{	
		getline(f, test);
		if(test.find(b)!=string::npos)
		{
			ok=true;
			break;
		}
		else n++;
	}
	if(!ok) 
		return -1;	

	return position;
}

long int LocaliseMO(ifstream& f, int& n, int& n_MO, string b)
{
	n=0;
	n_MO=0;
	string m="MO Number";
	long int position;
	position=LocaliseBlock(f,n,b);
	f.clear();
	f.seekg(position);
	string test;

	while(!f.eof())
	{
		getline(f,test);
		if(test.find(m)!=string::npos)
			n_MO++;
		else if(test.find(b)!=string::npos)
			break;
	}
	if(n==0 || n_MO==0)
		position=-1;

	return position;
}

void WFX::read_file_wfx(ifstream& file)
{
	string a="Title";
	string b="Keywords";
	string c="Number of Nuclei";
	string d="Number of Primitives";
	string e="Number of Occupied Molecular Orbitals";
	string f="Number of Perturbations";
	string g="Nuclear Names";
	string h="Atomic Number";
	string i="Nuclear Charges";
	string j="Nuclear Cartesian Coordinates";
	string k="Net Charge";
	string l="Number of Electrons";
	string m="Number of Alpha Electrons";
	string n="Number of Beta Electrons";
	string o="Electronic Spin Multiplicity";
	string p="Model";
	string q="Primitive Centers";
	string r="Primitive Types";
	string s="Primitive Exponents";
	string t="Molecular Orbital Occupation Numbers";
	string u="Molecular Orbital Energies";
	string v="Molecular Orbital Spin Types";
	string w="Molecular Orbital Primitive Coefficients";
	string x="Energy = T + Vne + Vee + Vnn";
	string y="Virial Ratio (-V/T)";
	string z="Nuclear Cartesian Energy Gradients";
	string aa="Nuclear Virial of Energy-Gradient-Based Forces on Nuclei, W";
	string bb="Full Virial Ratio, -(V - W)/T";
	string cc="Number of Core Electrons";
	string dd="Additional Electron Density Function (EDF)";

	_Title=read_string(file, a, true);
	_Keywords=read_string(file, b, true);
	_Number_of_Nuclei=read_int(file, c, true);
	_Number_of_Primitives=read_int(file, d, true);
	_Number_of_Occupied_Molecular_Orbital=read_int(file, e, true);
	_Number_of_Perturbations=read_int(file,f, true);
	_Nuclear_Names=read_one_block_string(file,g, true);
	_Atomic_Number=read_one_block_int(file,h, false,_Number_of_Nuclei);
	_Nuclear_Charges=read_one_block_real(file,i, true);
	_Nuclear_Cartesian_Coordinates=read_one_block_real(file,j, true, _Number_of_Nuclei*3);
	_Net_Charge=read_real(file,k, true); 
	_Number_of_Electrons=read_int(file,l, true);
	_Number_of_Alpha_Electrons=read_int(file,m, true);
	_Number_of_Beta_Electrons=read_int(file,n, true);
	_Electronic_Spin_Multiplicity=read_int(file,o, false);
	_Model=read_string(file,p, false);
	_Primitive_Centers=read_one_block_int(file,q, true, _Number_of_Primitives);
	_Primitive_Types=read_one_block_int(file,r, true, _Number_of_Primitives);
	_Primitive_Exponents=read_one_block_real(file,s, true, _Number_of_Primitives);
	_Molecular_Orbital_Occupation_Numbers=read_one_block_real(file,t, true);
	_Molecular_Orbital_Energies=read_one_block_real(file,u, true);
	_Molecular_Orbital_Spin_Types=read_one_block_string(file,v, true);
	_Molecular_Orbital_Primitive_Coefficients=read_MOPC_block(file,w, true);
	_Energy=read_real(file,x, true);
	_Virial_Ratio=read_real(file,y, true);
	_Nuclear_Cartesian_Energy_Gradients=read_NCEG_block(file,z, false);
	_Nuclear_Virial_of_Energy_Gradient_Based_Forces_on_Nuclei=read_real(file,aa, false);
	_Full_Virial_Ratio=read_real(file,bb, false);
	_Number_of_Core_Electrons=read_int(file,cc, false);
	_Additionnal_Electron_Density_Function=read_AEDF_block(file,dd, false);
}

void WFX::write_one_block_int(ofstream& f, vector<int> v, string b, bool r, int nint)
{
	if(!r && v.size()==0)
		return;
	f<<b<<endl;
	int i;
	for(i=0; i<int(v.size()); i++)
		f<<" "<<v[i];
	b.insert(1,"\\");
	f<<endl<<b<<endl;;
}

void WFX::write_one_block_real(ofstream& f, vector<double> v, string b, bool r)
{
	if(!r && v.size()==0)
		return;
	f<<b<<endl;
	int i;
	for(i=0; i<int(v.size()); i++)
		f<<" "<<v[i]<<"\t";
	b.insert(1,"\\");
	f<<endl<<b<<endl;
}

void WFX::write_one_block_real(ofstream& f, vector<double> v, string b, bool r, int ndouble)
{
	if(!r && v.size()==0)
		return;
	f<<b<<endl;
	int i;
	for(i=0; i<ndouble; i++)
		f<<" "<<v[i]<<"\t";
	b.insert(1,"\\");
	f<<endl<<b<<endl;
}

void WFX::write_one_block_string(ofstream& f, vector<string> v, string b, bool r)
{
	if(!r && v.size()==0)
		return;
	f<<b<<endl;
	int i;
	for(i=0; i<int(v.size()); i++)
		f<<v[i]<<endl;
	b.insert(1,"\\");
	f<<b<<endl;
}

void WFX::write_int(ofstream& f, int i, string b, bool r)
{
	if(!r && i==0)
		return;
	f<<b<<endl;
	f<<" "<<i<<endl;
	b.insert(1,"\\");
	f<<b<<endl;
}

void WFX::write_real(ofstream& f, double d, string b, bool r)
{
	if(!r && d==0.0)
		return;
	f<<b<<endl;
	f<<" "<<d<<endl;
	b.insert(1,"\\");
	f<<b<<endl;
}

void WFX::write_string(ofstream& f, string s, string b, bool r)
{
	if(!r && s=="None")
		return;
	f<<b<<endl;
	f<<" "<<s<<endl;
	b.insert(1,"\\");
	f<<b<<endl;
}

void WFX::write_MOPC_block(ofstream& f, vector<MOPC> v, bool r)
{
	if(!r && v.size()==0)
		return;
	string b="<Molecular Orbital Primitive Coefficients>";
	string m1="<MO Number>";
	string m2=m1;
	m2.insert(1,"\\");
	f<<b<<endl;
	int i,j;
	int n=_Number_of_Primitives;
	int nMO=v.back().MO_Number();
	for(i=0; i<nMO; i++)
	{
		f<<m1<<endl;
		f<<" "<<v[i].MO_Number()<<endl;
		f<<m2<<endl;
		for(j=0; j<n; j++)
		{
			f<<" "<<v[i].Coefficients()[j];
		}
		f<<endl;
	}
	b.insert(1,"\\");
	f<<b<<endl;
}

void WFX::write_NCEG_block(ofstream& f, vector<NCEG> v, bool r)
{
	if(!r && v.size()==0)
		return;
	string b="<Nuclear Cartesian Energy Gradients>";
	f<<b<<endl;
	int i,j;
	for(i=0; i<int(v.size()); i++)
	{
		f<<" "<<v[i].symbol();
		for(j=0; j<3; j++)
			f<<" "<<v[i].gradient()[i];
	}
	b.insert(1,"\\");
	f<<endl<<b<<endl;
}

void WFX::write_AEDF_block(ofstream& f, AEDF a, bool r)
{
	if(!r && a.EDF_Primitives_Coefficients().size()==0)
		return;
	string b="<Additional Electron Density Function (EDF)>";
	string ub1="<Number of EDF Primitives>";
	string ub2="<EDF Primitive Centers>";
	string ub3="<EDF Primitive Types>";
	string ub4="<EDF Primitive Exponents>";
	string ub5="<EDF Primitive Coefficients>";
	f<<b<<endl;
	write_int(f, a.Number_of_EDF_Primitives(), ub1, r);
	write_one_block_int(f, a.EDF_Primitives_Centers(), ub2, r, a.Number_of_EDF_Primitives());
	write_one_block_int(f, a.EDF_Primitives_Types(), ub3, r, a.Number_of_EDF_Primitives());
	write_one_block_real(f, a.EDF_Primitives_Exponents(), ub4, r);
	write_one_block_real(f, a.EDF_Primitives_Coefficients(), ub5, r);
	b.insert(1,"\\");
	f<<b<<endl;
}

void WFX::write_file_wfx(ofstream& file)
{
	cout<<std::scientific;
	cout<<std::setprecision(15);
	//cout<<std::left<<std::setw(20);


	string a="<Title>";
	string b="<Keywords>";
	string c="<Number of Nuclei>";
	string d="<Number of Primitives>";
	string e="<Number of Occupied Molecular Orbitals>";
	string f="<Number of Perturbations>";
	string g="<Nuclear Names>";
	string h="<Atomic Number>";
	string i="<Nuclear Charges>";
	string j="<Nuclear Cartesian Coordinates>";
	string k="<Net Charge>";
	string l="<Number of Electrons>";
	string m="<Number of Alpha Electrons>";
	string n="<Number of Beta Electrons>";
	string o="<Electronic Spin Multiplicity>";
	string p="<Model>";
	string q="<Primitive Centers>";
	string r="<Primitive Types>";
	string s="<Primitive Exponents>";
	string t="<Molecular Orbital Occupation Numbers>";
	string u="<Molecular Orbital Energies>";
	string v="<Molecular Orbital Spin Types>";
	string w="<Energy = T + Vne + Vee + Vnn>";
	string x="<Virial Ratio (-V/T)>";
	string y="<Nuclear Virial of Energy-Gradient-Based Forces on Nuclei, W>";
	string z="<Full Virial Ratio, -(V - W)/T>";
	string aa="<Number of Core Electrons>";

	write_string(file, _Title, a, true);
	write_string(file, _Keywords, b, true);
	write_int(file, _Number_of_Nuclei, c, true);
	write_int(file, _Number_of_Primitives, d, true);
	write_int(file, _Number_of_Occupied_Molecular_Orbital, e, true);
	write_int(file, _Number_of_Perturbations, f, true);
	write_one_block_string(file, _Nuclear_Names, g, true);
	write_one_block_int(file, _Atomic_Number, h, false, _Number_of_Nuclei);
	write_one_block_real(file, _Nuclear_Charges, i, true);
	write_one_block_real(file, _Nuclear_Cartesian_Coordinates, j, true, _Number_of_Nuclei*3);
	write_real(file, _Net_Charge, k, true); 
	write_int(file, _Number_of_Electrons, l, true);
	write_int(file, _Number_of_Alpha_Electrons, m, true);
	write_int(file, _Number_of_Beta_Electrons, n, true);
	write_int(file, _Electronic_Spin_Multiplicity, o, false);
	write_string(file, _Model, p, false);
	write_one_block_int(file, _Primitive_Centers, q, true, _Number_of_Primitives);
	write_one_block_int(file, _Primitive_Types, r, true, _Number_of_Primitives);
	write_one_block_real(file, _Primitive_Exponents, s, true);
	write_one_block_real(file, _Molecular_Orbital_Occupation_Numbers, t, true);
	write_one_block_real(file, _Molecular_Orbital_Energies, u, true);
	write_one_block_string(file, _Molecular_Orbital_Spin_Types, v, true);
	write_MOPC_block(file, _Molecular_Orbital_Primitive_Coefficients, true);
	write_real(file, _Energy, w, true);
	write_real(file, _Virial_Ratio, x, true);
	write_NCEG_block(file, _Nuclear_Cartesian_Energy_Gradients, false);
	write_real(file, _Nuclear_Virial_of_Energy_Gradient_Based_Forces_on_Nuclei, y, false);
	write_real(file, _Full_Virial_Ratio, z, false);
	write_int(file, _Number_of_Core_Electrons, aa, false);
	write_AEDF_block(file, _Additionnal_Electron_Density_Function, false);
}