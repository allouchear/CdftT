#include<iostream>
#include<sstream>
#include <Utils/MOLDENGAB.h>
#include <Common/Constants.h>

using namespace std;

MOLDENGAB::MOLDENGAB()
{
	_symbol=vector<string> ();
	_atomic_number=vector<int> ();
	_coord=vector<vector<double>> ();
	_shell_types=vector<string> ();
	_L_types=vector<int> ();
	_exposants=vector<double> ();
	_number_of_gtf=vector<int> ();
	_cgtf_coefs=vector<double> ();
	_factor_coefs=vector<double> ();
	_MO_energy=vector<double> ();
	_MO_coefs=vector<vector<double>> ();
	_occupation=vector<double> ();
	_spin_types=vector<string> ();
	_coord_type="None";
	_number_of_atoms=0;
	_number_of_MO_coefs=0;
	_number_of_MO=0;
	_alpha_and_beta=true;
	_basis_or_gto="None";
	_alpha_occupation=vector<double> ();
	_beta_occupation=vector<double> ();
	_alpha_MO_coefs=vector<vector<double>> ();
	_beta_MO_coefs=vector<vector<double>> ();
	_alpha_energies=vector<double> ();
	_beta_energies=vector<double> ();
	_format="None";
	_cart_sphe="none";
	_mixte=false;
}

MOLDENGAB::MOLDENGAB(ifstream& file)
{
	_symbol=vector<string> ();
	_atomic_number=vector<int> ();
	_coord=vector<vector<double>> ();
	_shell_types=vector<string> ();
	_L_types=vector<int> ();
	_exposants=vector<double> ();
	_number_of_gtf=vector<int> ();
	_num_center=vector<int> ();
	_cgtf_coefs=vector<double> ();
	_factor_coefs=vector<double> ();
	_MO_energy=vector<double> ();
	_alpha_energies=vector<double> ();
	_beta_energies=vector<double> ();
	_MO_coefs=vector<vector<double>> ();
	_alpha_MO_coefs=vector<vector<double>> ();
	_beta_MO_coefs=vector<vector<double>> ();
	_occupation=vector<double> ();
	_alpha_occupation=vector<double> ();
	_beta_occupation=vector<double> ();
	_spin_types=vector<string> ();
	_coord_type="None";
	_number_of_MO_coefs=0;
	_number_of_MO=0;
	_alpha_and_beta=true;
	_mixte=false;

	file.clear();
	file.seekg(0,file.beg);
	string p;
	getline(file,p);

	if(p.find("[Molden Format]")!=string::npos)
	{
		_format="molden";
		_basis_or_gto="[GTO]";
	}

	else if(p.find("[Gabedit Format]")!=string::npos)
	{
		_format="gabedit";
		_basis_or_gto="[Basis]";
		if(p.find("Sphe")!=string::npos)
			_cart_sphe="sphe";
		else if(p.find("Cart")!=string::npos)
			_cart_sphe="cart";
		else
		{
			cout<<"Error, can't recognize data format (sphe/cart)."<<endl;
			cout<<"Please check your file."<<endl;
			exit(1);
		}
	}

	else
	{
		cout<<"Error, can't recognize file format."<<endl;
		cout<<"Please check your file."<<endl;
		exit(1);
	}

	read_atom_data(file);

	_number_of_atoms=_atomic_number.size();
	_n_at_basis=vector<int> (_number_of_atoms);
	
	read_basis_data(file);
	read_MO_data(file);

	if(_coord_type=="Angs")
		for(int i=0; i<_number_of_atoms; i++)
			for(int j=0; j<3; j++)
				_coord[i][j]*=ANGTOBOHR;

	if(_alpha_and_beta)
	{
		_alpha_energies=_beta_energies=_MO_energy;
		_alpha_occupation=_beta_occupation=_occupation;
		_alpha_MO_coefs=_beta_MO_coefs=_MO_coefs;
	}

	else
		for(int i=0; i<_number_of_MO; i++)
		{
			if(_spin_types[i]=="Alpha")
			{
				_alpha_energies.push_back(_MO_energy[i]);
				_alpha_occupation.push_back(_occupation[i]);
				_alpha_MO_coefs.push_back(_MO_coefs[i]);
			}
			
			else if(_spin_types[i]=="Beta")
			{
				_beta_energies.push_back(_MO_energy[i]);
				_beta_occupation.push_back(_occupation[i]);
				_beta_MO_coefs.push_back(_MO_coefs[i]);
			}
		}
}

void MOLDENGAB::read_atom_data(ifstream& f)
{
	string p;
	int an;
	vector<double> c (3);
	long int pos=LocaliseDataMolGabBefore(f,"[Atoms]");

	if(pos==-1)
	{
		cout<<"Atoms data not found"<<endl;
		cout<<"Data required, please check your file"<<endl;
		exit(1);
	}

	f.seekg(pos);
	getline(f,p);
	stringstream t(p);
	t>>p;
	t>>_coord_type;
	
	getline(f,p);

	do{
		stringstream s(p);
		s>>p;	
		_symbol.push_back(p);
		s>>p;
		s>>an;
		_atomic_number.push_back(an);
		s>>c[0];
		s>>c[1];
		s>>c[2];
		_coord.push_back(c);
		getline(f,p);
	}while(p.find("[")==string::npos);
}

void MOLDENGAB::read_basis_data(ifstream& f)
{
	string p;

	long int pos=LocaliseDataMolGab(f, _basis_or_gto);

	if(pos==-1)
	{
		cout<<"Basis (GTO) data not found"<<endl;
		cout<<"Data required, please check your file"<<endl;
		exit(1);
	}

	f.seekg(pos);
	int n= _atomic_number.size();

	for(int i=0; i<n; i++)
		read_one_basis_data(f);

	int m=_shell_types.size();
	
	long int posd, posf, posg;

	if(_format=="molden")
	{
		posd = LocaliseDataMolGab(f, "[5D]");
		posf = LocaliseDataMolGab(f, "[7F]");
		posg = LocaliseDataMolGab(f, "[9G]");

		if((posd==-1 && (posf!=-1 || posg!=1)) || (posf==-1 && (posd!=-1 || posg!=1)) || (posg==-1 && (posf!=-1 || posd!=1)))
			_mixte=true;
	}

	else if(_format=="gabedit" && _cart_sphe=="sphe")
		posd = posf = posg = 0;

	else
		posd = posf = posg =-1;

	for(int i=0; i<m; i++)
	{
		if(_shell_types[i]=="s" || _shell_types[i]=="S")
		{
			_L_types.push_back(0);
			_number_of_MO_coefs+=1;
		}

		else if(_shell_types[i]=="p" || _shell_types[i]=="P")
		{
			_L_types.push_back(1);
			_number_of_MO_coefs+=3;
		}

		else if((_shell_types[i]=="d" || _shell_types[i]=="D") && posd==-1)
		{
			_L_types.push_back(2);
			_number_of_MO_coefs+=6;
		}
		
		else if((_shell_types[i]=="d" || _shell_types[i]=="D") && posd!=-1)
		{
			_L_types.push_back(-2);
			_number_of_MO_coefs+=5;
		}
		
		else if((_shell_types[i]=="f" || _shell_types[i]=="F") && posf==-1)
		{
			_L_types.push_back(3);
			_number_of_MO_coefs+=10;
		}

		else if((_shell_types[i]=="f" || _shell_types[i]=="F") && posf!=-1)
		{
			_L_types.push_back(-3);
			_number_of_MO_coefs+=7;
		}
		
		else if((_shell_types[i]=="g" || _shell_types[i]=="G") && posg==-1)
		{
			_L_types.push_back(4);
			_number_of_MO_coefs+=15;
		}

		else if((_shell_types[i]=="g" || _shell_types[i]=="G") && posg!=-1)
		{
			_L_types.push_back(-4);
			_number_of_MO_coefs+=9;
		}

		else if(_format=="molden")
		{
			int N=int(toupper(_shell_types[i][0]))-int('F')+3;
			_L_types.push_back(-N);
			_number_of_MO_coefs+= 2*abs(N)+1;
		}

		else if(_cart_sphe=="sphe")
		{
			int N=int(toupper(_shell_types[i][0]))-int('F')+3;
			_L_types.push_back(-N);
			_number_of_MO_coefs+= 2*abs(N)+1;
		}

		else if(_cart_sphe=="cart")
		{
			int N=int(toupper(_shell_types[i][0]))-int('F')+3;
			_L_types.push_back(N);
			_number_of_MO_coefs+= 2*N+1;
		}

		else
		{
			cout<<"Error, shell type no recognize."<<endl;
			cout<<"Please check your file"<<endl;
			exit(1);
		}
	}

	_number_of_MO=_number_of_MO_coefs;
}

void MOLDENGAB::read_one_basis_data(istream& f)
{
	string p, t;
	int n,m;
	int v=0;
	double pc,c;

	getline(f,p);
	stringstream k(p);
	k>>m;
	getline(f,p);
	
	do{
		_num_center.push_back(m);
		stringstream s(p);
		s>>p;
		_shell_types.push_back(p);
		s>>n;
		_number_of_gtf.push_back(n);
		s>>c;

		for(int i=0; i<n; i++)
		{
			v++;
			_factor_coefs.push_back(c);
			getline(f,p);
			if(p.find("D")!=string::npos)
				p.replace(p.find("D"),1,"E");
			else if(p.find("d")!=string::npos)
				p.replace(p.find("d"),1,"e");

			if(p.find("D")!=string::npos)
				p.replace(p.find("D"),1,"E");
			else if(p.find("d")!=string::npos)
				p.replace(p.find("d"),1,"e");

			stringstream ss(p);
			ss>>pc;
			_exposants.push_back(pc);
			ss>>pc;
			_cgtf_coefs.push_back(pc);
		}

		getline(f,p);
		t=p;
		while(t.find(" ")!=string::npos)
			t.erase(t.find(" "),1);

	}while(!t.empty());
	m--;
	_n_at_basis[m]=v;
}

void MOLDENGAB::read_MO_data(ifstream& f)
{
	string p,t;
	double a;
	vector<double> aa;

	if(LocaliseDataMolGab(f,"Spin= Beta")!=-1)
	{
		_alpha_and_beta=false;
	}

	long int pos=LocaliseDataMolGab(f, "[MO]");

	if(pos==-1)
	{
		cout<<"Basis (GTO) data not found"<<endl;
		cout<<"Data required, please check your file"<<endl;
		exit(1);
	}

	f.seekg(pos);
	getline(f,p);

	do{
		for(int j=0; j<4; j++)
		{
			if(j!=0)
				getline(f,p);
			stringstream s(p);
			s>>p;

			if(p=="Ene=")
			{
				s>>a;
				_MO_energy.push_back(a);
			}
			else if(p=="Spin=")
			{
				s>>p;
				_spin_types.push_back(p);
			}
			else if(p=="Occup=")
			{
				s>>a;
				_occupation.push_back(a);
			}
		}

		for(int k=0; k<_number_of_MO_coefs; k++)
		{
			getline(f,p);
			stringstream ss(p);
			ss>>p;
			ss>>a;
			aa.push_back(a);
		}
		_MO_coefs.push_back(aa);
		aa=vector<double> ();
		getline(f,p);
		t=p;
		while(t.find(" ")!=string::npos)
			t.erase(t.find(" "),1);
	}while(!t.empty());

	_number_of_MO=_MO_energy.size();
}

void MOLDENGAB::PrintData()
{
	cout<<"Number of atoms = "<<_number_of_atoms<<endl;
	for(size_t i=0; i<_symbol.size(); i++)
		cout<<"Symbol "<<i<<" = "<<_symbol[i]<<endl;
	for(size_t i=0; i<_atomic_number.size(); i++)
		cout<<"Atomic number "<<i<<" = "<<_atomic_number[i]<<endl;
	for(size_t i=0; i<_coord.size(); i++)
		for(size_t j=0; j<_coord[i].size(); j++)
			cout<<"Coordinates "<<j<<" for atom "<<i<<" = "<<_coord[i][j]<<endl;
	for(size_t i=0; i<_num_center.size(); i++)
		cout<<"Num center "<<i<<" = "<<_num_center[i]<<endl;
	for(size_t i=0; i<_n_at_basis.size(); i++)
		cout<<"N at basis "<<i<<" = "<<_n_at_basis[i]<<endl;
	for(size_t i=0; i<_shell_types.size(); i++)
		cout<<"Shell type "<<i<<" = "<<_shell_types[i]<<endl;
	for(size_t i=0; i<_L_types.size(); i++)
		cout<<"L type "<<i<<" = "<<_L_types[i]<<endl;
	for(size_t i=0; i<_exposants.size(); i++)
		cout<<"Exposant "<<i<<" = "<<_exposants[i]<<endl;
	for(size_t i=0; i<_number_of_gtf.size(); i++)
		cout<<"Number of GTF "<<i<<" = "<<_number_of_gtf[i]<<endl;
	for(size_t i=0; i<_cgtf_coefs.size(); i++)
		cout<<"GTF coefficient "<<i<<" = "<<_cgtf_coefs[i]<<endl;
	for(size_t i=0; i<_factor_coefs.size(); i++)
		cout<<"CGTF coefficient "<<i<<" = "<<_factor_coefs[i]<<endl;
	for(size_t i=0; i<_alpha_energies.size(); i++)
		cout<<"Alpha MO energy "<<i<<" = "<<_MO_energy[i]<<endl;
	for(size_t i=0; i<_beta_energies.size(); i++)
		cout<<"Beta MO energy "<<i<<" = "<<_MO_energy[i]<<endl;
	for(size_t i=0; i<_alpha_MO_coefs.size(); i++)
		for(size_t j=0; j<_alpha_MO_coefs[i].size(); j++)
			cout<<"Alpha MO coefficient ["<<i<<"]["<<j<<"] = "<<_alpha_MO_coefs[i][j]<<endl;
	for(size_t i=0; i<_beta_MO_coefs.size(); i++)
		for(size_t j=0; j<_beta_MO_coefs[i].size(); j++)
			cout<<"Beta MO coefficient ["<<i<<"]["<<j<<"] = "<<_beta_MO_coefs[i][j]<<endl;
	for(size_t i=0; i<_alpha_occupation.size(); i++)
		cout<<"Alpha Occupation "<<i<<" = "<<_alpha_occupation[i]<<endl;
	for(size_t i=0; i<_beta_occupation.size(); i++)
		cout<<"Beta Occupation "<<i<<" = "<<_beta_occupation[i]<<endl;
	for(size_t i=0; i<_spin_types.size(); i++)
		cout<<"Spin types "<<i<<" = "<<_spin_types[i]<<endl;
	cout<<"Coordinates type = "<<_coord_type<<endl;
	cout<<"Number of MO = "<<_number_of_MO<<endl;
	cout<<"Number of MO coefficients = "<<_number_of_MO_coefs<<endl;

	if(_alpha_and_beta)
		cout<<"Alpha == Beta"<<endl;
	else
		cout<<"Alpha != Beta"<<endl;
}

long int LocaliseDataMolGab(ifstream& f, string b)
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
			ok=true;
			position=f.tellg();
			break;
		}
	}

	if(!ok) 
		return -1;	

	return position;
}

long int LocaliseDataMolGabBefore(ifstream& f, string b)
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
