#ifndef CDFTT_FCHK_H_INCLUDED
#define CDFTT_FCHK_H_INCLUDED

#include<iostream>
#include<string>
#include<vector>
#include <Utils/Utils.h>

using namespace std;

class FCHK{

private:
	int _number_of_electrons;
	int _number_of_alpha_electrons;
	int _number_of_beta_electrons;
	int _number_of_basis_functions;
	vector<int> _atomic_numbers;
	vector<double> _nuclear_charges;
	vector<double> _current_cartesian_coordinates; // or vector<vector<double>>
	vector<double> _primitive_exponents;
	vector<double> _contraction_coefficients; // à voir
	vector<double> _alpha_orbital_energies;
	vector<double> _alpha_mo_coefficients; // or vector<vector<double>>
	vector<double> _beta_orbital_energies;
	vector<double> _beta_mo_coefficients; // or vector<vector<double>>
	double _total_energy;
	vector<double> _mulliken_charges;
	vector<double> _npa_charges;
	vector<double> _dipole_moment;

	vector<int> _shell_types;
	vector<int> _shell_to_atom_map;
	vector<int> _number_of_primitives_per_shell;
	vector<double> _sp_contraction_coefficients;
	vector<double> _coordinates_for_shells;
	vector<double> _alpha_occupation;
	vector<double> _beta_occupation;
	int _number_of_atoms;
	int _number_of_contracted_shells;
	int _number_of_primitive_shells;
	int _highest_angular_momentum;

	int _ok_alpha;
	bool _alpha_and_beta;

public:
	FCHK();
	FCHK(ifstream&);
	~FCHK() {}
	int NumberOfElectrons() {return _number_of_electrons;}
	int NumberOfAlphaElectrons() {return _number_of_alpha_electrons;}
	int NumberOfBetaElectrons() {return _number_of_beta_electrons;}
	int NumberOfBasisFunctions() {return _number_of_basis_functions;}
	vector<int> AtomicNumbers() {return _atomic_numbers;}
	vector<double> NuclearCharges() {return _nuclear_charges;}
	vector<double> CurrentCartesianCoordinates() {return _current_cartesian_coordinates;}
	vector<double> PrimitiveExponents() {return _primitive_exponents;}
	vector<double> ContractionCoefficients() {return _contraction_coefficients;}
	vector<double> AlphaOrbitalEnergies() {return _alpha_orbital_energies;}
	vector<double> AlphaMOCoefficients() {return _alpha_mo_coefficients;}
	vector<double> BetaOrbitalEnergies() {return _beta_orbital_energies;}
	vector<double> BetaMOCoefficients() {return _beta_mo_coefficients;}
	double TotalEnergy() {return _total_energy;}
	vector<double> MullikenCharges() {return _mulliken_charges;}
	vector<double> NpaCharges() {return _npa_charges;}
	vector<double> DipoleMoment() {return _dipole_moment;}
	vector<int> ShellTypes() {return _shell_types;}
	vector<int> ShellToAtomMap() {return _shell_to_atom_map;}
	vector<int> NumberOfPrimitivesPerShell() {return _number_of_primitives_per_shell;}
	vector<double> spContractionCoefficients() {return _sp_contraction_coefficients;}
	vector<double> CoordinatesForShells() {return _coordinates_for_shells;}
	vector<double> AlphaOccupation() {return _alpha_occupation;}
	vector<double> BetaOccupation() {return _beta_occupation;}
	int NumberOfPrimitivesShells() {return _number_of_primitive_shells;}
	int NumberOfAtoms() {return _number_of_atoms;}
	int NumberOfContractedShells() {return _number_of_contracted_shells;}
	int HighestAngularMomentum() {return _highest_angular_momentum;}
	bool AlphaAndBeta() {return _alpha_and_beta;}
	int read_one_int(ifstream&, string);
	double read_one_real(ifstream&, string);
	vector<int> read_one_block_int(ifstream&, string);
	vector<double> read_one_block_real(ifstream&, string);
	void read_file_fchk(ifstream&);
	void PrintData();
};

long int LocaliseData(ifstream&, string);

#endif