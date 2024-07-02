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

	vector<int> _num_center;
	vector<string> _symbol;
	vector<int> _atomic_numbers;
	vector<vector<double>> _coordinates; //Prendre la derniere standard, sinon prendre input  /!\ en angstrom, a convertir

	vector<double> _mulliken_charges;
	vector<vector<double>> _MO_coefficients;
	vector<double> _MO_energy;
	
	vector<string> _shell_types;
	vector<int> _l_types;
	vector<double> _exposants;
	vector<double> _gtf_coefficients;
	vector<double> _cgtf_coefficients;

public:
	LOG();
	LOG(ifstream&);
	~LOG() {}

	void read_atoms_data(ifstream&);
	void read_basis_data(ifstream&);
};

long int LocaliseDataLog(ifstream&, string);
long int LocaliseDataLogBefore(ifstream&, string);
long int LocaliseNextDataLog(ifstream&, string);
long int LocaliseNextDataLogBefore(ifstream&, string);

#endif