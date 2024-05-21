#include<iostream>
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

vector<int> WFX::read_one_block_int(ifstream& f, string b)
{
	int data;
	vector<int> block_of_data(0);
	long int pos;
	int i,n;
	pos=LocaliseBlock(f, n, b);
	if(pos==-1)
	{
		cout<<b+" data not found"<<endl;
		return vector<int> (0);
	}
	f.seekg(pos);
	for(i=0; i<n; i++)
	{
		f>>data;
		block_of_data.push_back(data);
	}

	return block_of_data;
}

vector<double> WFX::read_one_block_real(ifstream& f, string b)
{
	double data;
	vector<double> block_of_data(0);
	long int pos;
	int i,n;
	pos=LocaliseBlock(f, n, b);
	if(pos==-1)
	{
		cout<<b+" data not found"<<endl;
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

vector<string> WFX::read_one_block_string(ifstream& f, string b)
{
	string data;
	vector<string> block_of_data(0);
	long int pos;
	int i,n;
	pos=LocaliseBlock(f, n, b);
	if(pos==-1)
	{
		cout<<b+" data not found"<<endl;
		return vector<string> (0);
	}
	f.seekg(pos);
	for(i=0; i<n; i++)
	{
		f>>data;
		block_of_data.push_back(data);
	}

	return block_of_data;
}

int WFX::read_int(ifstream& f, string b)
{
	int data,i,n;
	long int pos;
	pos=LocaliseBlock(f, n, b);
	if(pos==-1)
	{
		cout<<b+" data not found"<<endl;
		return 0;
	}
	f.seekg(pos);
	for(i=0; i<n; i++)
		f>>data;

	return data;
}

double WFX::read_real(ifstream& f, string b)
{
	double data;
	int i,n;
	long int pos;
	pos=LocaliseBlock(f, n, b);
	if(pos==-1)
	{
		cout<<b+" data not found"<<endl;
		return 0.0;
	}
	f.seekg(pos);
	for(i=0; i<n; i++)
		f>>data;

	return data;
}

string WFX::read_string(ifstream& f, string b)
{
	string data;
	int i,n;
	long int pos;
	pos=LocaliseBlock(f, n, b);
	if(pos==-1)
	{
		cout<<b+" data not found"<<endl;
		return "None";
	}
	f.seekg(pos);
	for(i=0; i<n; i++)
		f>>data;

	return data;
}

vector<MOPC> WFX::read_MOPC_block(ifstream& f, string b)
{
	int i,j,n,nMO,ndat;
	double c;
	vector<double> cdat(0);
	MOPC data;
	vector<MOPC> block_of_data(0);
	string MO="MO Number";
	long int pos;
	pos=LocaliseBlocks_MO(f,n,nMO,b);
	if(pos==-1)
	{
		cout<<b+" data not found"<<endl;
		return vector<MOPC> (0);
	}
	f.seekg(pos);
	for(i=0; i<nMO/2; i++)
	{
		f>>ndat;
		for(j=0; j<n/nMO; j++)
		{
			f>>c;
			cdat.push_back(c);
		}
		data.push_back(ndat,cdat);
		block_of_data.push_back(data);
	}

	return block_of_data;
}

vector<NCEG> WFX::read_NCEG_block(ifstream& f, string b)
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
		cout<<b+" data not found"<<endl;
		return vector<NCEG> (0);
	}

	f.seekg(pos);
	for(i=0; i<n/4; i++)
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

AEDF WFX::read_AEDF_block(ifstream& f, string b)
{
	string snp="Number of EDF Primitives";
	string spc="EDF Primitive Centers";
	string spt="EDF Primitive Types";
	string spe="EDF Primitive Exponents";
	string spcoef="EDF Primitive Coefficients";
	int n=read_int(f, snp);
	vector<int> pc=read_one_block_int(f, spc);
	vector<int> pt=read_one_block_int(f, spt);
	vector<double> pe=read_one_block_real(f, spe);
	vector<double> pcoef=read_one_block_real(f, spcoef);

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
			break;
		}
	}
	position=f.tellg();
	if(ok)
	{
		while(!f.eof())
		{
			getline(f,test);
			if(test.find(b)==string::npos)
				n++;
			else
				break;
		}
	}
	if(n==0)
		position=-1;

	return position;
}

long int LocaliseBlocks_MO(ifstream& f, int& n, int& n_MO, string b)
{
	n=0;
	n_MO=0;
	string m="MO Number";
	long int position;
	f.clear();
	f.seekg(0,f.beg);
	string test;
	bool ok=false;
	while(f.eof())
	{
		f>>test;
		if(test.find(b)!=string::npos)
		{
			ok=true;
			break;
		}
	}
	position=f.tellg();
	if(ok)
		while(f.eof())
		{
			f>>test;
			if(test.find(m)==string::npos)
				n_MO++;
			else if(test.find(b))
				break;
			else
				n++;
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
	string y="Virial ratio (-V/T)";
	string z="Nuclear Cartesian Energy Gradients";
	string aa="Nuclear Virial of Energy-Gradient-Based Forces on Nuclei, W";
	string bb="Full Virial Ratio, -(V - W)/T";
	string cc="Number of Core Electrons";
	string dd="Additional Electron Density Function (EDF)";

	_Title=read_string(file, a);
	_Keywords=read_string(file, b);
	_Number_of_Nuclei=read_int(file, c);
	_Number_of_Primitives=read_int(file, d);
	_Number_of_Occupied_Molecular_Orbital=read_int(file, e);
	_Number_of_Perturbations=read_int(file,f);
	_Nuclear_Names=read_one_block_string(file,g);
	_Atomic_Number=read_one_block_int(file,h);
	_Nuclear_Charges=read_one_block_real(file,i);
	_Nuclear_Cartesian_Coordinates=read_one_block_real(file,j);
	_Net_Charge=read_real(file,k); 
	_Number_of_Electrons=read_int(file,l);
	_Number_of_Alpha_Electrons=read_int(file,m);
	_Number_of_Beta_Electrons=read_int(file,n);
	_Electronic_Spin_Multiplicity=read_int(file,o);
	_Model=read_string(file,p);
	_Primitive_Centers=read_one_block_int(file,q);
	_Primitive_Types=read_one_block_int(file,r);
	_Primitive_Exponents=read_one_block_real(file,s);
	_Molecular_Orbital_Occupation_Numbers=read_one_block_real(file,t);
	_Molecular_Orbital_Energies=read_one_block_real(file,u);
	_Molecular_Orbital_Spin_Types=read_one_block_string(file,v);
	_Molecular_Orbital_Primitive_Coefficients=read_MOPC_block(file,w);
	_Energy=read_real(file,x);
	_Virial_Ratio=read_real(file,y);
	_Nuclear_Cartesian_Energy_Gradients=read_NCEG_block(file,z);
	_Nuclear_Virial_of_Energy_Gradient_Based_Forces_on_Nuclei=read_real(file,aa);
	_Full_Virial_Ratio=read_real(file,bb);
	_Number_of_Core_Electrons=read_int(file,cc);
	_Additionnal_Electron_Density_Function=read_AEDF_block(file,dd);
}