#include<iostream>
#include<sstream>
#include<common/Constants.h>
#include<analytic/Utils/LOG.h>

using namespace std;

LOG::LOG()
{
	_number_of_basis_functions=0;
	_number_of_cartesian_basis_functions=0;
	_number_of_primitive_gaussians=0;
	_number_of_alpha_electrons=0;
	_number_of_beta_electrons=0;
	_num_center=vector<int> ();
	_symbol=vector<string> ();
	_atomic_numbers=vector<int> ();
	_coordinates=vector<vector<double>> (); //Prendre la derniere standard, sinon prendre input  /!\ en angstrom, a convertir
	_mulliken_charges=vector<double> ();
	_MO_coefficients=vector<vector<double>> ();
	_MO_energy=vector<double> ();
	_shell_types=vector<string> ();
	_l_types=vector<int> ();
	_exposants=vector<double> ();
	_gtf_coefficients=vector<double> ();
	_cgtf_coefficients=vector<double> ();
}

LOG::LOG(ifstream& file)
{
	cout<<"test"<<endl;
}

void LOG::read_atoms_data(ifstream& f)
{
	long int pos, pos2;
	int n;
	double d;
	string p;
	_coordinates=vector<vector<double>> (3);

	pos=LocaliseDataLogBefore(f, "NAtoms=");
	f.seekg(pos);
	getline(f,p);
	f>>p;
	f>>_number_of_atoms;

	pos=LocaliseDataLog(f, "Standard orientation:");

	f.clear();
	f.seekg(0,f.beg);		//Cas où pos=-1;

	if(pos==-1)
		do{
			pos2=pos;
			pos=LocaliseNextDataLog(f,"Input orientation:");
		}while(pos!=-1);

	else
		do{
			pos2=pos;
			pos=LocaliseNextDataLog(f,"Standard orientation:");
		}while(pos!=-1);
	
	f.clear();
	f.seekg(pos2);
	getline(f,p);
	getline(f,p);
	getline(f,p);
	getline(f,p);

	for(int i=0; i<_number_of_atoms; i++)
	{
		getline(f,p);
		stringstream s(p);
		s>>n;
		_num_center.push_back(n);
		s>>n;
		_atomic_numbers.push_back(n);
		s>>p;

		for(int i=0; i<3; i++)
		{
			s>>d;
			d*=ANGTOBOHR;
			_coordinates[i].push_back(d);
		}
	}
}

void LOG::read_basis_data(ifstream& f)
{
	long int pos, pos2;
	int Nat, n;
	double d;
	string p;

	pos=LocaliseDataLog(f, "AO basis set in the form of general basis input (Overlap normalization):");
	
	do{
		pos2=pos;
		pos=LocaliseNextDataLog(f, "AO basis set in the form of general basis input (Overlap normalization):");
	}while(pos!=-1);

	f.clear();
	f.seekg(pos2);
	getline(f,p);
	stringstream s(p);
	s>>Nat;

	getline(f,p);

	int m=0;

	for(int i=0; i<_number_of_atoms; i++)
	{
		do
		{
			stringstream ss(p);
			ss>>p;
			_shell_types.push_back(p);
			ss>>n;
			ss>>d;
			_cgtf_coefficients.push_back(d);

			for(int j=0; j<n; i++)
			{
				getline(f,p);
				stringstream sss(p);
				sss>>d;
				_exposants.push_back(d);
				sss>>d;
				_gtf_coefficients.push_back(d);
				if(_shell_types[m]=="SP" || _shell_types[m]=="sp")					//A voir comment on organise les données !!!
				{
					sss>>d;
					_gtf_coefficients.push_back(d);
				}
			}
			getline(f,p);
			m++;
		}while(p.find("*")==string::npos);
	}

	pos=LocaliseNextDataLogBefore(f, "basis functions,");
	f.clear();
	f.seekg(pos);
	getline(f,p);

	

}

long int LocaliseDataLog(ifstream& f, string b)
{
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
			position=f.tellg();
			ok=true;
			break;
		}
	}

	if(!ok) 
		return -1;	

	return position;
}

long int LocaliseDataLogBefore(ifstream& f, string b)
{
	long int position;
	f.clear();
	f.seekg(0,f.beg);
	string test;
	bool ok=false;
	while(!f.eof())
	{	
		position=f.tellg();
		getline(f, test);
		if(test.find(b)!=string::npos)
		{
			ok=true;
			break;
		}
	}

	if(!ok) 
		return -1;	

	return position;
}

long int LocaliseNextDataLog(ifstream& f, string b)
{
	long int position;
	f.clear();
	string test;
	bool ok=false;
	while(!f.eof())
	{	
		getline(f, test);
		if(test.find(b)!=string::npos)
		{
			position=f.tellg();
			ok=true;
			break;
		}
	}

	if(!ok) 
		return -1;	

	return position;
}

long int LocaliseNextDataLogBefore(ifstream& f, string b)
{
	long int position;
	f.clear();
	string test;
	bool ok=false;
	while(!f.eof())
	{
		position=f.tellg();	
		getline(f, test);
		if(test.find(b)!=string::npos)
		{
			ok=true;
			break;
		}
	}

	if(!ok) 
		return -1;	

	return position;
}