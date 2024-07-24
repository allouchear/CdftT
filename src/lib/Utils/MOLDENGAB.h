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
		//! A default constructor.
		/*! This constructor is used to set all of the parameters on 0 or "None" value. */
	MOLDENGAB();

		//! A constructor taking one argument.
		/*! This constructor is used to set all of the parameters with the data in the file. */
	MOLDENGAB(ifstream&);

		//! A default desctructor.
		/*! We don't use it. */
	~MOLDENGAB() {}

		//! A normal member taking one argument and returning a void value.
		/*! Set all the atomic data on the value in the file. */
	void read_atom_data(ifstream& f);

		//! A normal member taking one argument and returning a void value.
		/*! Set all the basis data on the value in the file. */
	void read_basis_data(ifstream& f);

		//! A normal member taking one argument and returning a void value.
		/*! Set all the molecular orbitals data on the value in the file. */
	void read_MO_data(ifstream& f);

		//! A normal member taking one argument and returning a void value.
		/*! Read one basis data in the file. */
	void read_one_basis_data(istream& f);

		//! A normal member taking one argument and returning a void value.
		/*! Print all the data save in the object. */
	void PrintData();

		//! A normal member taking no arguments and returning a vector<string> value.
		/*! \return The table of symbol of each atoms. */
	vector<string> Symbol() {return _symbol;}

		//! A normal member taking no arguments and returning a vector<int> value.
		/*! \return The table of atomic number of each atoms. */
	vector<int> AtomicNumbers() {return _atomic_number;}

		//! A normal member taking no arguments and returning a vector<vector<double>> value.
		/*! \return The table of coordinates of each atoms. */
	vector<vector<double>> Coordinates() {return _coord;}

		//! A normal member taking no arguments and returning a vector<string> value.
		/*! \return The table of shell type of each shell. */
	vector<string> ShellTypes() {return _shell_types;}

		//! A normal member taking no arguments and returning a vector<int> value.
		/*! \return The table of L type of each shell. */
	vector<int> Ltypes() {return _L_types;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of exponents of each CGTF. */
	vector<double> Exposants() {return _exposants;}

		//! A normal member taking no arguments and returning a vector<int> value.
		/*! \return The table of number of GTF for each centers (atoms). */
	vector<int> NumberOfGtf() {return _number_of_gtf;}

		//! A normal member taking no arguments and returning a vector<int> value.
		/*! \return The table of num center (the size of this table is equal to the number of atoms). */
	vector<int> Numcenter() {return _num_center;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of CGTF coefficients. */
	vector<double> CgtfCoefficients() {return _cgtf_coefs;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of CGTF factor coefficients. */
	vector<double> FactorCoefficients() {return _factor_coefs;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of molecular orbitals energy. */
	vector<double> MO_Energy() {return _MO_energy;}

		//! A normal member taking no arguments and returning a vector<vector<double>> value.
		/*! \return The table of CGTF factor coefficients. */
	vector<vector<double>> MO_coefficients() {return _MO_coefs;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of occupation number for each molecular orbital. */
	vector<double> OccupationNumber() {return _occupation;}

		//! A normal member taking no arguments and returning a vector<string> value.
		/*! \return The table of spin type (Alpha or Beta) for each molecular orbital. */
	vector<string> SpinTypes() {return _spin_types;}

		//! A normal member taking no arguments and returning a string value.
		/*! \return The unis of coordinates (Angs or AU). */
	string CoordinatesType() {return _coord_type;}

		//! A normal member taking no arguments and returning a string value.
		/*! \return If their is Basis or GTO (because their is different in a .gab or a .molden). */
	string BasisOrGTO() {return _basis_or_gto;}

		//! A normal member taking no arguments and returning an int value.
		/*! \return The number of atoms. */
	int NumberOfAtoms() {return _number_of_atoms;}

		//! A normal member taking no arguments and returning an int value.
		/*! \return The number of molecular orbital coefficients. */
	int NumberOfMOCoefs() {return _number_of_MO_coefs;}

		//! A normal member taking no arguments and returning an int value.
		/*! \return The number of molecular orbitals. */
	int NumberOfMO() {return _number_of_MO;}

		//! A normal member taking no arguments and returning a boolean value.
		/*! \return The boolean of alpha and beta (true if alpha and beta have the same molecular orbitals coefficients). */
	bool AlphaAndBeta() {return _alpha_and_beta;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of alpha occupation number for each molecular orbital. */
	vector<double> AlphaOccupation() {return _alpha_occupation;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of beta occupation number for each molecular orbital. */
	vector<double> BetaOccupation() {return _beta_occupation;}

		//! A normal member taking no arguments and returning a vector<vector<double>> value.
		/*! \return The table of alpha molecular orbitals coefficients for each molecular orbital. */
	vector<vector<double>> AlphaMOCoefs() {return _alpha_MO_coefs;}

		//! A normal member taking no arguments and returning a vector<vector<double>> value.
		/*! \return The table of beta molecular orbitals coefficients for each molecular orbital. */
	vector<vector<double>> BetaMOCoefs() {return _beta_MO_coefs;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of alpha molecular orbitals energy for each molecular orbital. */
	vector<double> AlphaEnergies() {return _alpha_energies;}

		//! A normal member taking no arguments and returning a vector<double> value.
		/*! \return The table of alpha molecular orbitals energy for each molecular orbital. */
	vector<double> BetaEnergies() {return _beta_energies;}

		//! A normal member taking no arguments and returning a string value.
		/*! \return The format of the file (gabedit or molden). */
	string Format() {return _format;}

		//! A normal member taking no arguments and returning a string value.
		/*! \return The type of coordinates for each shell */
	string CartOrSphe() {return _cart_sphe;}

		//! A normal member taking no arguments and returning a vector<int> value.
		/*! \return The table of number of gtf in each centers. */
	vector<int> NatBasis() {return _n_at_basis;}
};

	//! A function taking two arguments and returning a long int value.
	/*! \return The position after the string search in the file. */
long int LocaliseDataMolGab(ifstream& f, string b);

	//! A function taking three arguments and returning a long int value.
	/*! \return The position before the string search in the file. */
long int LocaliseDataMolGabBefore(ifstream& f, string b);

#endif
