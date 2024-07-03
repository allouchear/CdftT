#ifndef CDFTT_LOG_H_INCLUDED
#define CDFTT_LOG_H_INCLUDED

#include<iostream>
#include<vector>
#include<fstream>
#include<string>

using namespace std;

class LOG{

private:
	int _number_of_atoms;
	int _number_of_basis_functions;
	int _number_of_cartesian_basis_functions;
	int _number_of_primitive_gaussians;
	int _number_of_alpha_electrons;
	int _number_of_beta_electrons;
	double _energy;

	vector<int> _num_center;
	vector<string> _symbol;
	vector<int> _atomic_numbers;
	vector<vector<double>> _coordinates;

	vector<double> _mulliken_charges;
	
	vector<string> _shell_types;
	vector<int> _l_types;
	vector<int> _number_of_gtf;
	vector<double> _exposants;
	vector<double> _cgtf_coefficients;
	vector<double> _cgtf_sp_coefficients;
	vector<double> _factor_coefficients;

	vector<double> _alpha_occupation;
	vector<double> _beta_occupation;
	vector<vector<double>> _alpha_MO_coefs;
	vector<vector<double>> _beta_MO_coefs;
	vector<double> _alpha_energy;
	vector<double> _beta_energy;

	string _d_cart_sphe;
	string _f_cart_sphe;

	int _number_of_MO;
	int _number_of_MO_coefs;
	vector<int> _n_at_basis;

	bool _alpha_and_beta;

public:
	LOG();
	LOG(ifstream&);
	~LOG() {}

	int NumberOfAtoms() {return _number_of_atoms;}
	int NumberOfBasisFunctions() {return _number_of_basis_functions;}
	int NumberOfCartesianBasisFunctions() {return _number_of_cartesian_basis_functions;}
	int NumberOfPrimitiveGaussians() {return _number_of_primitive_gaussians;}
	int NumberOfAlphaElectrons() {return _number_of_alpha_electrons;}
	int NumberOfBetaElectrons() {return _number_of_beta_electrons;}
	double Energy() {return _energy;}
	vector<int> NumCenter() {return _num_center;}
	vector<string> Symbol() {return _symbol;}
	vector<int> AtomicNumbers() {return _atomic_numbers;}
	vector<vector<double>> Coordinates() {return _coordinates;}
	vector<double> MullikenCharges() {return _mulliken_charges;}
	vector<string> ShellTypes() {return _shell_types;}
	vector<int> Ltypes() {return _l_types;}
	vector<int> NumberOfGtf() {return _number_of_gtf;}
	vector<double> Exposants() {return _exposants;}
	vector<double> CgtfCoefficients() {return _cgtf_coefficients;}
	vector<double> CgtfSpCoefficients() {return _cgtf_sp_coefficients;}
	vector<double> FactorCoefficients() {return _factor_coefficients;}
	vector<double> AlphaOccupation() {return _alpha_occupation;}
	vector<double> BetaOccupation() {return _beta_occupation;}
	vector<vector<double>> AlphaMOcoefs() {return _alpha_MO_coefs;}
	vector<vector<double>> BetaMOcoefs() {return _beta_MO_coefs;}
	vector<double> AlphaEnergy() {return _alpha_energy;}
	vector<double> BetaEnergy() {return _beta_energy;}
	string D_cart_sphe() {return _d_cart_sphe;}
	string F_cart_sphe() {return _f_cart_sphe;}
	int NumberOfMO() {return _number_of_MO;}
	int NumberOfMOcoefs() {return _number_of_MO_coefs;}
	vector<int> NatBasis() {return _n_at_basis;}
	bool AlphaAndBeta() {return _alpha_and_beta;}

	void read_atoms_data(ifstream&);
	void read_basis_data(ifstream&);
	void read_MO_data(ifstream&);
	void PrintData();
};

long int LocaliseDataLog(ifstream&, string);
long int LocaliseDataLogBefore(ifstream&, string);
long int LocaliseDataLogBefore(ifstream&, string, string);
long int LocaliseNextDataLog(ifstream&, string);
long int LocaliseNextDataLogBefore(ifstream&, string);
long int LocaliseNextDataLogBefore(ifstream&, string, string);

#endif