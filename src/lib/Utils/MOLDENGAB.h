#ifndef CDFTT_MOLDENGAB_H_INCLUDED
#define CDFTT_MOLDENGAB_H_INCLUDED

#include<iostream>
#include<fstream>
#include<string>
#include<vector>

using namespace std;

class MOLDENGAB{

private:
	vector<string> _symbol;
	vector<int> _atomic_number;
	vector<vector<double>> _coord;
	vector<string> _shell_types;
	vector<int> _L_types;
	vector<double> _exposants;
	vector<int> _number_of_gtf;
	vector<int> _num_center;
	vector<double> _cgtf_coefs;
	vector<double> _factor_coefs;
	vector<double> _MO_energy;
	vector<vector<double>> _MO_coefs;
	vector<double> _occupation;
	vector<string> _spin_types;
	string _coord_type;
	string _basis_or_gto;
	int _number_of_atoms;
	int _number_of_MO_coefs;
	int _number_of_MO;
	bool _alpha_and_beta; //true de base
	vector<int> _n_at_basis;

//Pour séparer alpha et beta
	vector<double> _alpha_occupation;
	vector<double> _beta_occupation;
	vector<vector<double>> _alpha_MO_coefs;
	vector<vector<double>> _beta_MO_coefs;
	vector<double> _alpha_energies;
	vector<double> _beta_energies;
//Pour un fichier gab
	string _format;
	string _cart_sphe;

public:
	MOLDENGAB();
	MOLDENGAB(ifstream&);
	~MOLDENGAB() {}
	void read_atom_data(ifstream& f);
	void read_basis_data(ifstream& f);
	void read_MO_data(ifstream& f);
	void read_one_basis_data(istream& f);
	void PrintData();

	vector<string> Symbol() {return _symbol;}
	vector<int> AtomicNumbers() {return _atomic_number;}
	vector<vector<double>> Coordinates() {return _coord;}
	vector<string> ShellTypes() {return _shell_types;}		//Pas utile
	vector<int> Ltypes() {return _L_types;}
	vector<double> Exposants() {return _exposants;}
	vector<int> NumberOfGtf() {return _number_of_gtf;}
	vector<int> Numcenter() {return _num_center;}
	vector<double> CgtfCoefficients() {return _cgtf_coefs;}
	vector<double> FactorCoefficients() {return _factor_coefs;}
	vector<double> MO_Energy() {return _MO_energy;}				 //Pas utile
	vector<vector<double>> MO_coefficients() {return _MO_coefs;} //Pas utile
	vector<double> OccupationNumber() {return _occupation;}		//Pas utile
	vector<string> SpinTypes() {return _spin_types;}	//Pas utile
	string CoordinatesType() {return _coord_type;}		//Pas utile
	string BasisOrGTO() {return _basis_or_gto;}			//Pas utile
	int NumberOfAtoms() {return _number_of_atoms;}
	int NumberOfMOCoefs() {return _number_of_MO_coefs;}	
	int NumberOfMO() {return _number_of_MO;}			
	bool AlphaAndBeta() {return _alpha_and_beta;}
	vector<double> AlphaOccupation() {return _alpha_occupation;}
	vector<double> BetaOccupation() {return _beta_occupation;}
	vector<vector<double>> AlphaMOCoefs() {return _alpha_MO_coefs;}
	vector<vector<double>> BetaMOCoefs() {return _beta_MO_coefs;}
	vector<double> AlphaEnergies() {return _alpha_energies;}
	vector<double> BetaEnergies() {return _beta_energies;}
	string Format() {return _format;}			//Pas utile
	string CartOrSphe() {return _cart_sphe;}	//Pas utile
	vector<int> NatBasis() {return _n_at_basis;}
};

long int LocaliseDataMolGab(ifstream& f, string b);
long int LocaliseDataMolGabBefore(ifstream& f, string b);

#endif