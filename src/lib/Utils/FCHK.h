#ifndef CDFTT_FCHK_H_INCLUDED
#define CDFTT_FCHK_H_INCLUDED

#include<iostream>
#include<string>
#include<vector>
#include <Utils/Utils.h>

using namespace std;

	//! A LOG class.
	/*! This class will be used to read in the log format. */

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
		//! A default constructor.
		/*! This constructor is used to set all of the parameters on 0 or "None" value. */
	FCHK();
	
		//! A constructor taking one argument.
		/*! This constructor is used to set all of the parameters with the data in the file. */
	FCHK(ifstream&);

		//! A default desctructor.
		/*! We don't use it. */
	~FCHK() {}

		//! A normal member taking no arguments and returning an int value.
		/*! \return The number of electrons. */
	int NumberOfElectrons() {return _number_of_electrons;}

		//! A normal member taking no arguments and returning an int value.
		/*! \return The number of alpha electrons. */
	int NumberOfAlphaElectrons() {return _number_of_alpha_electrons;}

		//! A normal member taking no arguments and returning an int value.
		/*! \return The number of beta electrons. */
	int NumberOfBetaElectrons() {return _number_of_beta_electrons;}

		//! A normal member taking no arguments and returning an int value.
		/*! \return The number of basis functions. */
	int NumberOfBasisFunctions() {return _number_of_basis_functions;}

		//! A normal member taking no arguments and returning a vector<int> value.
		/*! \return The table of atomic number of each atoms. */
	vector<int> AtomicNumbers() {return _atomic_numbers;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of nuclear charges of each atoms. */
	vector<double> NuclearCharges() {return _nuclear_charges;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of coordinates of each atoms. */
	vector<double> CurrentCartesianCoordinates() {return _current_cartesian_coordinates;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of exponents of each CGTF. */
	vector<double> PrimitiveExponents() {return _primitive_exponents;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of CGTF coefficients. */
	vector<double> ContractionCoefficients() {return _contraction_coefficients;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of alpha molecular orbitals energy for each molecular orbital. */
	vector<double> AlphaOrbitalEnergies() {return _alpha_orbital_energies;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of alpha molecular orbitals coefficients for each molecular orbital. */
	vector<double> AlphaMOCoefficients() {return _alpha_mo_coefficients;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of beta molecular orbitals energy for each molecular orbital. */
	vector<double> BetaOrbitalEnergies() {return _beta_orbital_energies;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of beta molecular orbitals coefficients for each molecular orbital. */
	vector<double> BetaMOCoefficients() {return _beta_mo_coefficients;}

		//! A normal member taking no arguments and returning a double value.
		/*! \return The total energy. */
	double TotalEnergy() {return _total_energy;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of Mulliken charges of each atoms. */
	vector<double> MullikenCharges() {return _mulliken_charges;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of NPA charges of each atoms. */
	vector<double> NpaCharges() {return _npa_charges;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of dipole moment of each atoms. */
	vector<double> DipoleMoment() {return _dipole_moment;}

		//! A normal member taking no arguments and returning a vector<int> value.
		/*! \return The table of shell type of each shell. */
	vector<int> ShellTypes() {return _shell_types;}

		//! A normal member taking no arguments and returning a vector<int> value.
		/*! \return The table of shell to atom map. */
	vector<int> ShellToAtomMap() {return _shell_to_atom_map;}

		//! A normal member taking no arguments and returning a vector<int> value.
		/*! \return The table of number of primitives per shell. */
	vector<int> NumberOfPrimitivesPerShell() {return _number_of_primitives_per_shell;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of CGTF coefficients (only for sp contraction). */
	vector<double> spContractionCoefficients() {return _sp_contraction_coefficients;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of coordinates for each shells. */
	vector<double> CoordinatesForShells() {return _coordinates_for_shells;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of alpha occupation number for each molecular orbital. */
	vector<double> AlphaOccupation() {return _alpha_occupation;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of beta occupation number for each molecular orbital. */
	vector<double> BetaOccupation() {return _beta_occupation;}

		//! A normal member taking no arguments and returning an int value.
		/*! \return The number of primitives shells. */
	int NumberOfPrimitivesShells() {return _number_of_primitive_shells;}

		//! A normal member taking no arguments and returning an int value.
		/*! \return The number of atoms. */
	int NumberOfAtoms() {return _number_of_atoms;}

		//! A normal member taking no arguments and returning an int value.
		/*! \return The number of contracted shells. */
	int NumberOfContractedShells() {return _number_of_contracted_shells;}

		//! A normal member taking no arguments and returning an int value.
		/*! \return The highest angular momentum. */
	int HighestAngularMomentum() {return _highest_angular_momentum;}

		//! A normal member taking no arguments and returning a boolean value.
		/*! \return The boolean of alpha and beta (true if alpha and beta have the same molecular orbitals coefficients). */
	bool AlphaAndBeta() {return _alpha_and_beta;}

		//! A normal member taking two arguments and returning an int value.
		/*! \return One int read. */
	int read_one_int(ifstream&, string);

		//! A normal member taking two arguments and returning a double value.
		/*! \return One real read. */
	double read_one_real(ifstream&, string);

		//! A normal member taking two arguments and returning a vector<int> value.
		/*! \return One block of int read. */
	vector<int> read_one_block_int(ifstream&, string);

		//! A normal member taking two arguments and returning a vector<double> value.
		/*! \return One block of real read. */
	vector<double> read_one_block_real(ifstream&, string);

		//! A normal member taking three arguments and returning a void value.
		/*! Read the file and set parameters on it. */
	void read_file_fchk(ifstream&);

		//! A normal member taking one argument and returning a void value.
		/*! Print all the data save in the object. */
	void PrintData();
};

	//! A function taking two arguments and returning a long int value.
	/*! \return The position after the string search in the file. */
long int LocaliseData(ifstream&, string);

#endif
