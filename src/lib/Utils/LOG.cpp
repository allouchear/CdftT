#include<iostream>
#include<sstream>
#include <Common/Constants.h>
#include <Utils/LOG.h>

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
	_coordinates=vector<vector<double>> ();
	_mulliken_charges=vector<double> ();
	_shell_types=vector<string> ();
	_l_types=vector<int> ();
	_exposants=vector<double> ();
	_cgtf_coefficients=vector<double> ();
	_factor_coefficients=vector<double> ();
	_alpha_occupation=vector<double> ();
	_beta_occupation=vector<double> ();
	_alpha_MO_coefs=vector<vector<double>> ();
	_beta_MO_coefs=vector<vector<double>> ();
	_alpha_energy=vector<double> ();
	_beta_energy=vector<double> ();
}

LOG::LOG(ifstream& file)
{
	_num_center=vector<int> ();
	_symbol=vector<string> ();
	_atomic_numbers=vector<int> ();
	_coordinates=vector<vector<double>> ();
	_mulliken_charges=vector<double> ();
	_shell_types=vector<string> ();
	_l_types=vector<int> ();
	_number_of_gtf=vector<int> ();
	_exposants=vector<double> ();
	_cgtf_coefficients=vector<double> ();
	_factor_coefficients=vector<double> ();
	_alpha_occupation=vector<double> ();
	_beta_occupation=vector<double> ();
	_alpha_MO_coefs=vector<vector<double>> ();
	_beta_MO_coefs=vector<vector<double>> ();
	_alpha_energy=vector<double> ();
	_beta_energy=vector<double> ();
	_n_at_basis=vector<int> ();
	_number_of_MO=0;

	read_atoms_data(file);
	read_basis_data(file);
	read_MO_data(file);
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
	stringstream sp(p);
	sp>>p;
	sp>>_number_of_atoms;
	
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
		s>>p;
		s>>n;
		_atomic_numbers.push_back(n);
		s>>p;

		for(int j=0; j<3; j++)
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
	int n,l;
	double d,c;
	string p;

	pos=LocaliseDataLog(f, "AO basis set in the form of general basis input (Overlap normalization):");
	
	do{
		pos2=pos;
		pos=LocaliseNextDataLog(f, "AO basis set in the form of general basis input (Overlap normalization):");
	}while(pos!=-1);

	f.clear();
	f.seekg(pos2);

	_n_at_basis=vector<int> (_number_of_atoms);
	int m=0;
	int v;

	for(int i=0; i<_number_of_atoms; i++)
	{
		getline(f,p);
		stringstream nc(p);
		nc>>l;
		getline(f,p);
		v=0;
		do
		{
			stringstream ss(p);
			ss>>p;
			_shell_types.push_back(p);
			ss>>n;
			_number_of_gtf.push_back(n);
			ss>>c;

			for(int j=0; j<n; j++)
			{
				v++;
				_num_center.push_back(l);
				_factor_coefficients.push_back(c);
				getline(f,p);

				if(p.find("D")!=string::npos)
					p.replace(p.find("D"),1,"E");
				else if(p.find("d")!=string::npos)
					p.replace(p.find("d"),1,"e");

				if(p.find("D")!=string::npos)
					p.replace(p.find("D"),1,"E");
				else if(p.find("d")!=string::npos)
					p.replace(p.find("d"),1,"e");

				if(_shell_types[m]=="SP" || _shell_types[m]=="sp")
				{
					if(p.find("D")!=string::npos)
						p.replace(p.find("D"),1,"E");
					else if(p.find("d")!=string::npos)
						p.replace(p.find("d"),1,"e");
				}

				stringstream sss(p);
				sss>>d;
				_exposants.push_back(d);
				sss>>d;
				_cgtf_coefficients.push_back(d);
				if(_shell_types[m]=="SP" || _shell_types[m]=="sp")					//A voir comment on organise les données !!!
				{
					sss>>d;
					_cgtf_sp_coefficients.push_back(d);
				}
				else
					_cgtf_sp_coefficients.push_back(0.0);
			}
			getline(f,p);
			m++;
		}while(p.find("*")==string::npos);
		_n_at_basis[i]=v;
	}

	pos=LocaliseNextDataLogBefore(f, "basis functions,");

	f.clear();
	f.seekg(pos);

	getline(f,p);

	stringstream ss(p);
	ss>>_number_of_basis_functions;
	ss>>p;
	ss>>p;
	ss>>_number_of_primitive_gaussians;
	ss>>p;
	ss>>p;
	ss>>_number_of_cartesian_basis_functions;

	getline(f,p);

	stringstream sss(p);
	sss>>_number_of_alpha_electrons;
	sss>>p;
	sss>>p;
	sss>>_number_of_beta_electrons;

	do{
		pos2=pos;
		pos=LocaliseNextDataLogBefore(f, "E(");
	}while(pos!=-1);

	f.clear();
	f.seekg(pos2);

	getline(f,p);
	stringstream ssss(p);
	ssss>>p;
	ssss>>p;
	ssss>>p;
	ssss>>p;
	ssss>>_energy;

	pos=LocaliseDataLogBefore(f, "Standard basis:", "General basis read from cards:");

	if(pos==-1)
	{
		_d_cart_sphe="sphe";
		_f_cart_sphe="sphe";
	}

	else
		do{
			pos2=pos;
			pos=LocaliseNextDataLogBefore(f, "Standard basis:", "General basis read from cards:");
		}while(pos!=-1);

	f.clear();
	f.seekg(pos2);

	getline(f,p);
	
	if(p.find("5D")!=string::npos || p.find("5d")!=string::npos)
		_d_cart_sphe="sphe";
	else if(p.find("6D")!=string::npos || p.find("6d")!=string::npos)
		_d_cart_sphe="cart";
	else
		_d_cart_sphe="sphe";

	if(p.find("7F")!=string::npos || p.find("7f")!=string::npos)
		_f_cart_sphe="sphe";
	else if(p.find("10F")!=string::npos || p.find("10f")!=string::npos)
		_f_cart_sphe="cart";
	else
		_f_cart_sphe="sphe";

	for(size_t i=0; i<_shell_types.size(); i++)
	{
		if(_shell_types[i]=="s" || _shell_types[i]=="S")
		{
			_l_types.push_back(0);
			_number_of_MO+=1;
		}

		else if(_shell_types[i]=="p" || _shell_types[i]=="P")
		{
			_l_types.push_back(1);
			_number_of_MO+=3;
		}

		else if(_shell_types[i]=="sp" || _shell_types[i]=="SP")
		{
			_l_types.push_back(-1);
			_number_of_MO+=4;
		}

		else if((_shell_types[i]=="d" || _shell_types[i]=="D") && _d_cart_sphe=="cart")
		{
			int N=int(toupper(_shell_types[i][0]))-int('D')+2;
			_l_types.push_back(N);
			_number_of_MO+= 2*N+1;
		}

		else if(_shell_types[i]=="d" || _shell_types[i]=="D")
		{
			int N=int(toupper(_shell_types[i][0]))-int('D')+2;
			_l_types.push_back(-N);
			_number_of_MO+= 2*N+1;
		}

		else if((_shell_types[i]=="f" || _shell_types[i]=="F") && _f_cart_sphe=="cart")
		{
			int N=int(toupper(_shell_types[i][0]))-int('D')+2;
			_l_types.push_back(N);
			_number_of_MO+= 2*N+1;
		}

		else if(_shell_types[i]=="f" || _shell_types[i]=="F")
		{
			int N=int(toupper(_shell_types[i][0]))-int('D')+2;
			_l_types.push_back(-N);
			_number_of_MO+= 2*N+1;
		}

		else
		{
			int N=int(toupper(_shell_types[i][0]))-int('D')+2;
			_l_types.push_back(-N);
			_number_of_MO+= 2*N+1;
		}
	}
	_number_of_MO_coefs=_number_of_MO;
	_beta_MO_coefs=_alpha_MO_coefs=vector<vector<double>> (_number_of_MO_coefs);
}

void LOG::read_MO_data(ifstream& f)
{
	long int pos, pos2;
	int n,m;
	double d;
	string p,p2,pp, name;

	pos=LocaliseDataLog(f, "Alpha Molecular Orbital Coefficients:");

	if(pos==-1)
	{
		_alpha_and_beta=true;
		pos=LocaliseDataLog(f, "Molecular Orbital Coefficients:");

		do{
			pos2=pos;
			pos=LocaliseNextDataLog(f, "Molecular Orbital Coefficients:");
		}while(pos!=-1);
	}

	else
	{
		_alpha_and_beta=false;
		do{
			pos2=pos;
			pos=LocaliseNextDataLog(f, "Alpha Molecular Orbital Coefficients:");
		}while(pos!=-1);
	}

	f.clear();
	f.seekg(pos2);

	getline(f,p);

	n=0;

	do{
		m=0;
		stringstream t(p);
		do{
			t>>p;
			m++;
		}while(!t.eof());

		getline(f,p);
		stringstream s(p);
		
		for(int i=0; i<m; i++)
		{
			s>>p;

			if(p.find("O")!=string::npos && _alpha_and_beta)
				_alpha_occupation.push_back(2.0);
			else if(p.find("O")!=string::npos)
				_alpha_occupation.push_back(1.0);
			else if(p.find("V")!=string::npos)
				_alpha_occupation.push_back(0.0);
		}

		getline(f,p);
		stringstream ss(p);
		ss>>p;
		ss>>p;

		for(int i=0; i<m; i++)
		{
			ss>>d;
			_alpha_energy.push_back(d);
		}

		for(int i=0; i<_number_of_MO; i++)
		{
			getline(f,p);
			stringstream sss(p);

			sss>>p2;

			if(p.find("1S")!=string::npos || p.find("1s")!=string::npos)
			{
				sss>>p;
				sss>>p;
				if(n==0)
					_symbol.push_back(p);
			}

			sss>>p;

			for(int j=0; j<m; j++)
			{
				sss>>d;
				_alpha_MO_coefs[n+j].push_back(d);
			}
		}
		n+=m;
		getline(f,p);
		stringstream nn;
		nn<<n+1;
		pp=nn.str();
	}while(p.find(pp)!=string::npos);

	if(!_alpha_and_beta)
	{
		f.clear();
		f.seekg(pos2);
		pos=LocaliseNextDataLog(f, "Beta Molecular Orbital Coefficients:");

		f.clear();
		f.seekg(pos);

		getline(f,p);

		n=0;

		do{
			m=0;
			stringstream t(p);
			do{
				t>>p;
				m++;
			}while(!t.eof());

			getline(f,p);
			stringstream s(p);
		
			for(int i=0; i<m; i++)
			{
				s>>p;

				if(p.find("O")!=string::npos)
					_beta_occupation.push_back(1.0);
				else if(p.find("V")!=string::npos)
					_beta_occupation.push_back(0.0);
			}

			getline(f,p);
			stringstream ss(p);
			ss>>p;
			ss>>p;

			for(int i=0; i<m; i++)
			{
				ss>>d;
				_beta_energy.push_back(d);
			}

			for(int i=0; i<_number_of_MO; i++)
			{
				getline(f,p);
				stringstream sss(p);

				if(p.find("1S")!=string::npos || p.find("1s")!=string::npos)
				{
					sss>>p;
					sss>>p;
				}

				sss>>p;
				sss>>p;
		
				for(int j=0; j<m; j++)
				{
					sss>>d;
					_beta_MO_coefs[n+j].push_back(d);
				}
			}
			n+=m;
			getline(f,p);
			stringstream nn;
			nn<<n+1;
			pp=nn.str();
		}while(p.find(pp)!=string::npos);
	}

	else
	{
		_beta_occupation=_alpha_occupation;
		_beta_MO_coefs=_alpha_MO_coefs;
		_beta_energy=_alpha_energy;
	}

	pos=LocaliseNextDataLog(f, "Mulliken charges:");

	f.clear();
	f.seekg(pos);

	getline(f,p);

	for(int i=0; i<_number_of_atoms; i++)
	{
		getline(f,p);
		stringstream mc(p);
		mc>>p;
		mc>>p;
		mc>>d;
		_mulliken_charges.push_back(d);
	}
}

void LOG::PrintData()
{
	cout<<"Number of atoms = "<<_number_of_atoms<<endl;
	cout<<"Number of basis functions = "<<_number_of_basis_functions<<endl;
	cout<<"Number of cartesian basis functions = "<<_number_of_cartesian_basis_functions<<endl;
	cout<<"Number of primitive gaussians = "<<_number_of_primitive_gaussians<<endl;
	cout<<"Number of alpha electrons = "<<_number_of_alpha_electrons<<endl;
	cout<<"Number of beta electrons = "<<_number_of_beta_electrons<<endl;
	cout<<"Energy = "<<_energy<<endl;

	for(size_t i=0; i<_num_center.size(); i++)
		cout<<"Number center ["<<i<<"] = "<<_num_center[i]<<endl;
	for(size_t i=0; i<_symbol.size(); i++)
		cout<<"Symbol ["<<i<<"] = "<<_symbol[i]<<endl;
	for(size_t i=0; i<_atomic_numbers.size(); i++)
		cout<<"Atomic number ["<<i<<"] = "<<_atomic_numbers[i]<<endl;
	for(size_t i=0; i<_coordinates.size(); i++)
	{
		cout<<"Coordinates atom ["<<i<<"] = ";
		for(size_t j=0; j<_coordinates[i].size(); j++)
			cout<<_coordinates[i][j]<<"    ";
		cout<<endl;
	}

	for(size_t i=0; i<_mulliken_charges.size(); i++)
		cout<<"Mulliken charges "<<i<<" = "<<_mulliken_charges[i]<<endl;
	for(size_t i=0; i<_shell_types.size(); i++)
		cout<<"Shell types "<<i<<" = "<<_shell_types[i]<<endl;
	for(size_t i=0; i<_l_types.size(); i++)
		cout<<"Ltypes ["<<i<<"] = "<<_l_types[i]<<endl;
	for(size_t i=0; i<_number_of_gtf.size(); i++)
		cout<<"Number of GTF in CGTF "<<i<<" = "<<_number_of_gtf[i]<<endl;
	for(size_t i=0; i<_exposants.size(); i++)
		cout<<"Exposant "<<i<<" = "<<_exposants[i]<<endl;
	for(size_t i=0; i<_cgtf_coefficients.size(); i++)
		cout<<"CGTF coefficient "<<i<<" = "<<_cgtf_coefficients[i]<<endl;
	for(size_t i=0; i<_cgtf_sp_coefficients.size(); i++)
		cout<<"CGTF SP coefficient "<<i<<" = "<<_cgtf_sp_coefficients[i]<<endl;
	for(size_t i=0; i<_factor_coefficients.size(); i++)
		cout<<"Factor coefficient "<<i<<" = "<<_factor_coefficients[i]<<endl;

	for(size_t i=0; i<_alpha_occupation.size(); i++)
		cout<<"Alpha occupation "<<i<<" = "<<_alpha_occupation[i]<<endl;
	for(size_t i=0; i<_beta_occupation.size(); i++)
		cout<<"Beta occupation "<<i<<" = "<<_beta_occupation[i]<<endl;
	for(size_t i=0; i<_alpha_MO_coefs.size(); i++)
		for(size_t j=0; j<_alpha_MO_coefs[i].size(); j++)
			cout<<"Alpha MO coefficients ["<<i<<"]["<<j<<"] = "<<_alpha_MO_coefs[i][j]<<endl;
	for(size_t i=0; i<_beta_MO_coefs.size(); i++)
		for(size_t j=0; j<_beta_MO_coefs[i].size(); j++)
			cout<<"Beta MO coefficients ["<<i<<"]["<<j<<"] = "<<_beta_MO_coefs[i][j]<<endl;
	for(size_t i=0; i<_alpha_energy.size(); i++)
		cout<<"Alpha energy "<<i<<" = "<<_alpha_energy[i]<<endl;
	for(size_t i=0; i<_beta_energy.size(); i++)
		cout<<"Beta energy "<<i<<" = "<<_beta_energy[i]<<endl;

	cout<<"D cart/sphe = "<<_d_cart_sphe<<endl;
	cout<<"F cart/sphe = "<<_f_cart_sphe<<endl;

	cout<<"Number of MO = "<<_number_of_MO<<endl;
	cout<<"Number of MO coefficients = "<<_number_of_MO_coefs<<endl;

	for(size_t i=0; i<_n_at_basis.size(); i++)
		cout<<"n in basis "<<i<<" = "<<_n_at_basis[i]<<endl;

	if(_alpha_and_beta)
		cout<<"Alpha == Beta"<<endl;
	else
		cout<<"Alpha != Beta"<<endl;
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

long int LocaliseDataLogBefore(ifstream& f, string b1, string b2)
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
		if(test.find(b1)!=string::npos || test.find(b2)!=string::npos)
		{	
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

long int LocaliseNextDataLogBefore(ifstream& f, string b1, string b2)
{
	long int position;
	f.clear();
	string test;
	bool ok=false;
	while(!f.eof())
	{
		position=f.tellg();	
		getline(f, test);
		if(test.find(b1)!=string::npos || test.find(b2)!=string::npos)
		{
			ok=true;
			break;
		}
	}

	if(!ok) 
		return -1;	

	return position;
}