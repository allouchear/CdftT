#include<iostream>
#include<sstream>
#include <Utils/FCHK.h>

using namespace std;

FCHK::FCHK()
{
	_number_of_electrons=0;
	_number_of_alpha_electrons=0;
	_number_of_beta_electrons=0;
	_number_of_basis_functions=0;
	_atomic_numbers=vector<int> ();
	_nuclear_charges=vector<double> ();
	_current_cartesian_coordinates=vector<double> (); // or vector<vector<double>>
	_primitive_exponents=vector<double> ();
	_contraction_coefficients=vector<double> (); // à voir
	_alpha_orbital_energies=vector<double> ();
	_alpha_mo_coefficients=vector<double> (); // or vector<vector<double>>
	_beta_orbital_energies=vector<double> ();
	_beta_mo_coefficients=vector<double> (); // or vector<vector<double>>
	_total_energy=0.0;
	_mulliken_charges=vector<double> ();
	_npa_charges=vector<double> ();
	_dipole_moment=vector<double> ();
	_shell_types=vector<int> ();
	_shell_to_atom_map=vector<int> ();
	_number_of_primitives_per_shell=vector<int> ();
	_sp_contraction_coefficients=vector<double> ();
	_coordinates_for_shells=vector<double> ();
	_alpha_occupation=vector<double>();
	_beta_occupation=vector<double>();
	_number_of_atoms=0;
	_number_of_contracted_shells=0;
	_number_of_primitive_shells=0;
	_highest_angular_momentum=0;

	_ok_alpha=0;
	_alpha_and_beta=false;
}

FCHK::FCHK(ifstream& file)
{
	_ok_alpha=0;
	_alpha_and_beta=false;
	read_file_fchk(file);
	if(_beta_orbital_energies==vector<double> ())
		_beta_orbital_energies=_alpha_orbital_energies;
	if(_beta_mo_coefficients==vector<double> ())
		_beta_mo_coefficients=_alpha_mo_coefficients;

	_alpha_occupation = vector<double> (_number_of_basis_functions, 0.0);
	_beta_occupation = vector<double> (_number_of_basis_functions, 0.0);

	for(int i=0; i<_number_of_basis_functions; i++)
	{
		if(i<_number_of_alpha_electrons)
			_alpha_occupation[i]=1.0;
		if(i<_number_of_beta_electrons)
			_beta_occupation[i]=1.0;
	}
}

int FCHK::read_one_int(ifstream& f, string b)
{
	string p;
	double data;

	long int pos=LocaliseData(f,b);

	if(pos==-1)
	{
		cout<<b+" : data not found"<<endl;
		cout<<"Data required, please check your file"<<endl;
		exit(1);
	}

	f.seekg(pos);
	getline(f,p);
	stringstream ss(p);

	do{
		ss>>p;
	}while(p.find("I")==string::npos);

	ss>>data;

	return data;
}

double FCHK::read_one_real(ifstream& f, string b)
{
	string p;
	double data;

	long int pos=LocaliseData(f,b);

	if(pos==-1)
	{
		cout<<b+" : data not found"<<endl;
		cout<<"Data required, please check your file"<<endl;
		exit(1);
	}

	f.seekg(pos);
	getline(f,p);
	stringstream ss(p);

	do{
		ss>>p;
	}while(p.find("R")==string::npos);

	ss>>data;

	return data;
}

vector<int> FCHK::read_one_block_int(ifstream& f, string b)
{
	string p;
	vector<int> data;
	int n,d;

	long int pos=LocaliseData(f,b);

	if(pos==-1)
	{
		cout<<b+" : data not found"<<endl;
		cout<<"Data required, please check your file"<<endl;
		exit(1);
	}

	f.seekg(pos);
	getline(f,p);
	stringstream ss(p);

	do{
		ss>>p;
	}while(p.find("N=")==string::npos);
		ss>>n;

	data=vector<int> (n);

	for(int i=0; i<n; i++)
	{
		f>>d;
		data[i]=d;
	}

	return data;
}

vector<double> FCHK::read_one_block_real(ifstream& f, string b)
{
	string p;
	vector<double> data;
	int n;
	double d;

	long int pos=LocaliseData(f,b);

	if(pos==-1)
	{
		cout<<b+" : data not found"<<endl;
		if((b=="Beta Orbital Energies" || b=="Beta MO coefficients") && _ok_alpha==2)
			return vector<double> ();
		else if(b=="NPA Charges")
			return vector<double> ();
		else
		{
			cout<<"Data required, please check your file"<<endl;
			exit(1);
		}
	}

	f.seekg(pos);
	getline(f,p);
	stringstream ss(p);

	do{
		ss>>p;
	}while(p.find("N=")==string::npos);
		ss>>n;

	data=vector<double> (n);

	for(int i=0; i<n; i++)
	{
		f>>d;
		data[i]=d;
	}

	if(b=="Alpha Orbital Energies" || b=="Alpha MO coefficients")
		_ok_alpha++;

	return data;
}

void FCHK::read_file_fchk(ifstream& file)
{
	string a = "Number of electrons";
	string b = "Number of alpha electrons";
	string c = "Number of beta electrons";
	string d = "Number of basis functions";
	string e = "Atomic numbers";
	string f = "Nuclear charges";
	string g = "Current cartesian coordinates";
	string h = "Primitive exponents";
	string i = "Contraction coefficients";
	string j = "Alpha Orbital Energies";
	string k = "Alpha MO coefficients";
	string l = "Beta Orbital Energies";
	string m = "Beta MO coefficients";
	string n = "Total Energy";
	string o = "Mulliken Charges";
	string p = "NPA Charges";
	string q = "Dipole Moment";
	string r = "Shell types";
	string s = "Shell to atom map";
	string t = "Number of primitives per shell";
	string u = "P(S=P) Contraction coefficients";
	string v = "Coordinates of each shell";
	string w = "Number of atoms";
	string x = "Number of contracted shells";
	string y = "Number of primitive shells";
	string z = "Highest angular momentum";

	_number_of_electrons=read_one_int(file, a);
	_number_of_alpha_electrons=read_one_int(file, b);
	_number_of_beta_electrons=read_one_int(file, c);
	_number_of_basis_functions=read_one_int(file, d);
	_atomic_numbers=read_one_block_int(file, e);
	_nuclear_charges=read_one_block_real(file, f);
	_current_cartesian_coordinates=read_one_block_real(file, g); // or vector<vector<double>>
	_primitive_exponents=read_one_block_real(file, h);
	_contraction_coefficients=read_one_block_real(file, i); // à voir
	_alpha_orbital_energies=read_one_block_real(file, j);
	_alpha_mo_coefficients=read_one_block_real(file, k); // or vector<vector<double>>
	_beta_orbital_energies=read_one_block_real(file, l);
	_beta_mo_coefficients=read_one_block_real(file, m); // or vector<vector<double>>
	_total_energy=read_one_real(file, n);
	_mulliken_charges=read_one_block_real(file, o);
	_npa_charges=read_one_block_real(file, p);
	_dipole_moment=read_one_block_real(file, q);
	_shell_types=read_one_block_int(file, r);
	_shell_to_atom_map=read_one_block_int(file, s);
	_number_of_primitives_per_shell=read_one_block_int(file, t);
	_sp_contraction_coefficients=read_one_block_real(file, u);
	_coordinates_for_shells=read_one_block_real(file, v);
	_number_of_atoms=read_one_int(file, w);
	_number_of_contracted_shells=read_one_int(file, x);
	_number_of_primitive_shells=read_one_int(file, y);
	_highest_angular_momentum=read_one_int(file, z);

	if((_beta_mo_coefficients.size()==0) && (_beta_orbital_energies.size()==0) && (_ok_alpha==2))
	{
		_beta_mo_coefficients=_alpha_mo_coefficients;
		_beta_orbital_energies=_alpha_orbital_energies;
		_alpha_and_beta=true;
	}	
}

long int LocaliseData(ifstream& f, string b)
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

void FCHK::PrintData()
{
	cout<<"Number of electrons : "<<_number_of_electrons<<endl;
	cout<<"Number of alpha electrons : "<<_number_of_alpha_electrons<<endl;
	cout<<"Number of beta electrons : "<<_number_of_beta_electrons<<endl;
	cout<<"Number of basis functions : "<<_number_of_basis_functions<<endl;
	for(size_t i=0; i<_atomic_numbers.size(); i++)
		cout<<"Atomic number "<<i<<" : "<<_atomic_numbers[i]<<endl;
	for(size_t i=0; i<_nuclear_charges.size(); i++)
		cout<<"Nuclear charges "<<i<<" : "<<_nuclear_charges[i]<<endl;
	for(int i=0; i<_number_of_atoms; i++)
	{
		cout<<"Current Cartesian Coordinates "<<i<<" : ";
		for(int j=0; j<3; j++)
			cout<<_current_cartesian_coordinates[j+3*i]<<"   ";
		cout<<endl;
	}
	for(size_t i=0; i<_primitive_exponents.size(); i++)
		cout<<"Primitive exponents "<<i<<" : "<<_primitive_exponents[i]<<endl;
	for(size_t i=0; i<_contraction_coefficients.size(); i++)
		cout<<"Contraction coefficients "<<i<<" : "<<_contraction_coefficients[i]<<endl;
	for(size_t i=0; i<_alpha_orbital_energies.size(); i++)
		cout<<"Alpha orbital energies "<<i<<" : "<<_alpha_orbital_energies[i]<<endl;
	for(size_t i=0; i<_alpha_mo_coefficients.size(); i++)
		cout<<"Alpha MO coefficients "<<i<<" : "<<_alpha_mo_coefficients[i]<<endl;
	for(size_t i=0; i<_beta_orbital_energies.size(); i++)
		cout<<"Beta orbital energies "<<i<<" : "<<_beta_orbital_energies[i]<<endl;
	for(size_t i=0; i<_beta_mo_coefficients.size(); i++)
		cout<<"Beta MO coefficients "<<i<<" : "<<_beta_mo_coefficients[i]<<endl;
	cout<<"Total energy : "<<_total_energy<<endl;
	for(size_t i=0; i<_mulliken_charges.size(); i++)
		cout<<"Mulliken charges "<<i<<" : "<<_mulliken_charges[i]<<endl;
	for(size_t i=0; i<_npa_charges.size(); i++)
		cout<<"NPA charges "<<i<<" : "<<_npa_charges[i]<<endl;
	for(size_t i=0; i<_dipole_moment.size(); i++)
		cout<<"Dipole moment "<<i<<" : "<<_dipole_moment[i]<<endl;
	for(size_t i=0; i<_shell_types.size(); i++)
		cout<<"Shell types "<<i<<" : "<<_shell_types[i]<<endl;
	for(size_t i=0; i<_shell_to_atom_map.size(); i++)
		cout<<"Shell to atom map "<<i<<" : "<<_shell_to_atom_map[i]<<endl;
	for(size_t i=0; i<_number_of_primitives_per_shell.size(); i++)
		cout<<"Number of primitives per shell "<<i<<" : "<<_number_of_primitives_per_shell[i]<<endl;
	for(size_t i=0; i<_sp_contraction_coefficients.size(); i++)
		cout<<"SP contraction coefficients "<<i<<" : "<<_sp_contraction_coefficients[i]<<endl;
	for(int i=0; i<_number_of_contracted_shells; i++)
	{
		cout<<"Coordinates for each shells "<<i<<" : ";
		for(int j=0; j<3; j++)
			cout<<_coordinates_for_shells[j+3*i]<<"   ";
		cout<<endl;
	}
	cout<<"Number of atoms : "<<_number_of_atoms<<endl;
	cout<<"Number of contracted shells : "<<_number_of_contracted_shells<<endl;
	cout<<"Number of primitive shells : "<<_number_of_primitive_shells<<endl;
	cout<<"Highest angular momentum : "<<_highest_angular_momentum<<endl;
}
