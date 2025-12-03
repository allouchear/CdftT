#ifndef CDFTT_LOG_H_INCLUDED
#define CDFTT_LOG_H_INCLUDED

#include<iostream>
#include<vector>
#include<fstream>
#include<string>

using namespace std;

    //! A LOG class.
    /*! This class will be used to read in the log format. */
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
    bool _mixte;

public:

        //! A default constructor.
        /*! This constructor is used to set all of the parameters on 0 or "None" value. */
    LOG();

        //! A constructor taking one argument.
        /*! This constructor is used to set all of the parameters with the data in the file. */
    LOG(ifstream&);

        //! A default desctructor.
        /*! We don't use it. */
    ~LOG() {}

        //! A normal member taking no arguments and returning an int value.
        /*! \return The number of atoms. */
    int NumberOfAtoms() {return _number_of_atoms;}

        //! A normal member taking no arguments and returning an int value.
        /*! \return The number of basis functions in the file. */
    int NumberOfBasisFunctions() {return _number_of_basis_functions;}

        //! A normal member taking no arguments and returning an int value.
        /*! \return The number of cartesian basis functions in the file. */
    int NumberOfCartesianBasisFunctions() {return _number_of_cartesian_basis_functions;}

        //! A normal member taking no arguments and returning an int value.
        /*! \return The number of primitive gaussians in the file. */
    int NumberOfPrimitiveGaussians() {return _number_of_primitive_gaussians;}

        //! A normal member taking no arguments and returning an int value.
        /*! \return The number of alpha electrons. */
    int NumberOfAlphaElectrons() {return _number_of_alpha_electrons;}

        //! A normal member taking no arguments and returning an int value.
        /*! \return The number of beta electrons. */
    int NumberOfBetaElectrons() {return _number_of_beta_electrons;}

        //! A normal member taking no arguments and returning a double value.
        /*! \return The total energy. */
    double Energy() {return _energy;}

        //! A normal member taking no arguments and returning a vector<int> value.
        /*! \return The table of num center (the size of this table is equal to the number of atoms). */
    vector<int> NumCenter() {return _num_center;}

        //! A normal member taking no arguments and returning a vector<string> value.
        /*! \return The table of symbol of each atoms. */
    vector<string> Symbol() {return _symbol;}

        //! A normal member taking no arguments and returning a vector<int> value.
        /*! \return The table of atomic number of each atoms. */
    vector<int> AtomicNumbers() {return _atomic_numbers;}

        //! A normal member taking no arguments and returning a vector<vector<double>> value.
        /*! \return The table of coordinates of each atoms. */
    vector<vector<double>> Coordinates() {return _coordinates;}

        //! A normal member taking no arguments and returning a vector<double> value.
        /*! \return The table of Mulliken charges of each atoms. */
    vector<double> MullikenCharges() {return _mulliken_charges;}

        //! A normal member taking no arguments and returning a vector<string> value.
        /*! \return The table of shell type of each shell. */
    vector<string> ShellTypes() {return _shell_types;}

        //! A normal member taking no arguments and returning a vector<int> value.
        /*! \return The table of L type of each shell. */
    vector<int> Ltypes() {return _l_types;}

        //! A normal member taking no arguments and returning a vector<int> value.
        /*! \return The table of number of GTF for each centers (atoms). */
    vector<int> NumberOfGtf() {return _number_of_gtf;}

        //! A normal member taking no arguments and returning a vector<double> value.
        /*! \return The table of exponents of each CGTF. */
    vector<double> Exposants() {return _exposants;}

        //! A normal member taking no arguments and returning a vector<double> value.
        /*! \return The table of CGTF coefficients. */
    vector<double> CgtfCoefficients() {return _cgtf_coefficients;}

        //! A normal member taking no arguments and returning a vector<double> value.
        /*! \return The table of CGTF coefficients (only for sp contraction). */
    vector<double> CgtfSpCoefficients() {return _cgtf_sp_coefficients;}

        //! A normal member taking no arguments and returning a vector<double> value.
        /*! \return The table of CGTF factor coefficients. */
    vector<double> FactorCoefficients() {return _factor_coefficients;}

        //! A normal member taking no arguments and returning a vector<double> value.
        /*! \return The table of alpha occupation number for each molecular orbital. */
    vector<double> AlphaOccupation() {return _alpha_occupation;}

        //! A normal member taking no arguments and returning a vector<double> value.
        /*! \return The table of beta occupation number for each molecular orbital. */
    vector<double> BetaOccupation() {return _beta_occupation;}

        //! A normal member taking no arguments and returning a vector<vector<double>> value.
        /*! \return The table of alpha molecular orbitals coefficients for each molecular orbital. */
    vector<vector<double>> AlphaMOcoefs() {return _alpha_MO_coefs;}

        //! A normal member taking no arguments and returning a vector<vector<double>> value.
        /*! \return The table of beta molecular orbitals coefficients for each molecular orbital. */
    vector<vector<double>> BetaMOcoefs() {return _beta_MO_coefs;}

        //! A normal member taking no arguments and returning a vector<double> value.
        /*! \return The table of alpha molecular orbitals energy for each molecular orbital. */
    vector<double> AlphaEnergy() {return _alpha_energy;}

        //! A normal member taking no arguments and returning a vector<double> value.
        /*! \return The table of alpha molecular orbitals energy for each molecular orbital. */
    vector<double> BetaEnergy() {return _beta_energy;}

        //! A normal member taking no arguments and returning a string value.
        /*! \return The type of coordinates for shell D. */
    string D_cart_sphe() {return _d_cart_sphe;}

        //! A normal member taking no arguments and returning a string value.
        /*! \return The type of coordinates for shell F. */
    string F_cart_sphe() {return _f_cart_sphe;}

        //! A normal member taking no arguments and returning an int value.
        /*! \return The number of molecular orbitals. */
    int NumberOfMO() {return _number_of_MO;}

        //! A normal member taking no arguments and returning an int value.
        /*! \return The number of molecular orbital coefficients. */
    int NumberOfMOcoefs() {return _number_of_MO_coefs;}

        //! A normal member taking no arguments and returning a vector<int> value.
        /*! \return The table of number of gtf in each centers. */
    vector<int> NatBasis() {return _n_at_basis;}

        //! A normal member taking no arguments and returning a boolean value.
        /*! \return The boolean of alpha and beta (true if alpha and beta have the same molecular orbitals coefficients). */
    bool AlphaAndBeta() {return _alpha_and_beta;}

        //! A normal member taking one argument and returning a void value.
        /*! Set all the atomic data on the value in the file. */
    void read_atoms_data(ifstream&);

        //! A normal member taking one argument and returning a void value.
        /*! Set all the basis data on the value in the file. */
    void read_basis_data(ifstream&);

        //! A normal member taking one argument and returning a void value.
        /*! Set all the molecular orbitals data on the value in the file. */
    void read_MO_data(ifstream&);

        //! A normal member taking no arguments and returning a void value.
        /*! Print all the data save in the object. */
    void PrintData();

        //! A normal member taking no arguments and returning a bool value.
        /*! \return If their is a mixte basis. (Exemple : D is cartisian and F is spherical). */
    bool Mixte() {return _mixte;}
};

    //! A function taking two arguments and returning a long int value.
    /*! \return The position after the string search in the file. */
long int LocaliseDataLog(ifstream&, string);

    //! A function taking two arguments and returning a long int value.
    /*! \return The position before the string search in the file. */
long int LocaliseDataLogBefore(ifstream&, string);

    //! A function taking two arguments and returning a long int value.
    /*! \return The position before the first string find in the file. */
long int LocaliseDataLogBefore(ifstream&, string, string);

    //! A function taking two arguments and returning a long int value.
    /*! \return The next position after the string search in the file. */
long int LocaliseNextDataLog(ifstream&, string);

    //! A function taking two arguments and returning a long int value.
    /*! \return The next position before the string search in the file. */
long int LocaliseNextDataLogBefore(ifstream&, string);

    //! A function taking two arguments and returning a long int value.
    /*! \return The position before the first string find in the file. */
long int LocaliseNextDataLogBefore(ifstream&, string, string);

#endif
