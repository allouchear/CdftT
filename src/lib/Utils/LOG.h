#ifndef CDFTT_LOG_H_INCLUDED
#define CDFTT_LOG_H_INCLUDED

#include <fstream>
#include <string>
#include <vector>

    //! A LOG class.
    /*! This class will be used to read in the log format. */
class LOG
{
    private:
        int _number_of_atoms;
        int _number_of_basis_functions;
        int _number_of_cartesian_basis_functions;
        int _number_of_primitive_gaussians;
        int _number_of_alpha_electrons;
        int _number_of_beta_electrons;
        double _energy;
        std::vector<int> _num_center;
        std::vector<std::string> _symbol;
        std::vector<int> _atomic_numbers;
        std::vector<std::vector<double>> _coordinates;
        std::vector<double> _mulliken_charges;
        std::vector<std::string> _shell_types;
        std::vector<int> _l_types;
        std::vector<int> _number_of_gtf;
        std::vector<double> _exposants;
        std::vector<double> _cgtf_coefficients;
        std::vector<double> _cgtf_sp_coefficients;
        std::vector<double> _factor_coefficients;
        std::vector<double> _alpha_occupation;
        std::vector<double> _beta_occupation;
        std::vector<std::vector<double>> _alpha_MO_coefs;
        std::vector<std::vector<double>> _beta_MO_coefs;
        std::vector<double> _alpha_energy;
        std::vector<double> _beta_energy;
        std::string _d_cart_sphe;
        std::string _f_cart_sphe;
        int _number_of_MO;
        int _number_of_MO_coefs;
        std::vector<int> _n_at_basis;
        bool _alpha_and_beta;
        bool _mixte;

        
    public:
            //! A default constructor.
            /*! This constructor is used to set all of the parameters on 0 or "None" value. */
        LOG();

            //! A constructor taking one argument.
            /*! This constructor is used to set all of the parameters with the data in the file. */
        LOG(std::ifstream&);

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

            //! A normal member taking no arguments and returning a std::vector<int> value.
            /*! \return The table of num center (the size of this table is equal to the number of atoms). */
        std::vector<int> NumCenter() {return _num_center;}

            //! A normal member taking no arguments and returning a std::vector<std::string> value.
            /*! \return The table of symbol of each atoms. */
        std::vector<std::string> Symbol() {return _symbol;}

            //! A normal member taking no arguments and returning a std::vector<int> value.
            /*! \return The table of atomic number of each atoms. */
        std::vector<int> AtomicNumbers() {return _atomic_numbers;}

            //! A normal member taking no arguments and returning a std::vector<std::vector<double>> value.
            /*! \return The table of coordinates of each atoms. */
        std::vector<std::vector<double>> Coordinates() {return _coordinates;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of Mulliken charges of each atoms. */
        std::vector<double> MullikenCharges() {return _mulliken_charges;}

            //! A normal member taking no arguments and returning a std::vector<std::string> value.
            /*! \return The table of shell type of each shell. */
        std::vector<std::string> ShellTypes() {return _shell_types;}

            //! A normal member taking no arguments and returning a std::vector<int> value.
            /*! \return The table of L type of each shell. */
        std::vector<int> Ltypes() {return _l_types;}

            //! A normal member taking no arguments and returning a std::vector<int> value.
            /*! \return The table of number of GTF for each centers (atoms). */
        std::vector<int> NumberOfGtf() {return _number_of_gtf;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of exponents of each CGTF. */
        std::vector<double> Exposants() {return _exposants;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of CGTF coefficients. */
        std::vector<double> CgtfCoefficients() {return _cgtf_coefficients;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of CGTF coefficients (only for sp contraction). */
        std::vector<double> CgtfSpCoefficients() {return _cgtf_sp_coefficients;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of CGTF factor coefficients. */
        std::vector<double> FactorCoefficients() {return _factor_coefficients;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of alpha occupation number for each molecular orbital. */
        std::vector<double> AlphaOccupation() {return _alpha_occupation;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of beta occupation number for each molecular orbital. */
        std::vector<double> BetaOccupation() {return _beta_occupation;}

            //! A normal member taking no arguments and returning a std::vector<std::vector<double>> value.
            /*! \return The table of alpha molecular orbitals coefficients for each molecular orbital. */
        std::vector<std::vector<double>> AlphaMOcoefs() {return _alpha_MO_coefs;}

            //! A normal member taking no arguments and returning a std::vector<std::vector<double>> value.
            /*! \return The table of beta molecular orbitals coefficients for each molecular orbital. */
        std::vector<std::vector<double>> BetaMOcoefs() {return _beta_MO_coefs;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of alpha molecular orbitals energy for each molecular orbital. */
        std::vector<double> AlphaEnergy() {return _alpha_energy;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of alpha molecular orbitals energy for each molecular orbital. */
        std::vector<double> BetaEnergy() {return _beta_energy;}

            //! A normal member taking no arguments and returning a std::string value.
            /*! \return The type of coordinates for shell D. */
        std::string D_cart_sphe() {return _d_cart_sphe;}

            //! A normal member taking no arguments and returning a std::string value.
            /*! \return The type of coordinates for shell F. */
        std::string F_cart_sphe() {return _f_cart_sphe;}

            //! A normal member taking no arguments and returning an int value.
            /*! \return The number of molecular orbitals. */
        int NumberOfMO() {return _number_of_MO;}

            //! A normal member taking no arguments and returning an int value.
            /*! \return The number of molecular orbital coefficients. */
        int NumberOfMOcoefs() {return _number_of_MO_coefs;}

            //! A normal member taking no arguments and returning a std::vector<int> value.
            /*! \return The table of number of gtf in each centers. */
        std::vector<int> NatBasis() {return _n_at_basis;}

            //! A normal member taking no arguments and returning a boolean value.
            /*! \return The boolean of alpha and beta (true if alpha and beta have the same molecular orbitals coefficients). */
        bool AlphaAndBeta() {return _alpha_and_beta;}

            //! A normal member taking one argument and returning a void value.
            /*! Set all the atomic data on the value in the file. */
        void read_atoms_data(std::ifstream&);

            //! A normal member taking one argument and returning a void value.
            /*! Set all the basis data on the value in the file. */
        void read_basis_data(std::ifstream&);

            //! A normal member taking one argument and returning a void value.
            /*! Set all the molecular orbitals data on the value in the file. */
        void read_MO_data(std::ifstream&);

            //! A normal member taking no arguments and returning a void value.
            /*! Print all the data save in the object. */
        void PrintData();

            //! A normal member taking no arguments and returning a bool value.
            /*! \return If their is a mixte basis. (Exemple : D is cartisian and F is spherical). */
        bool Mixte() {return _mixte;}
};

    //! A function taking two arguments and returning a long int value.
    /*! \return The position after the std::string search in the file. */
long int LocaliseDataLog(std::ifstream&, std::string);

    //! A function taking two arguments and returning a long int value.
    /*! \return The position before the std::string search in the file. */
long int LocaliseDataLogBefore(std::ifstream&, std::string);

    //! A function taking two arguments and returning a long int value.
    /*! \return The position before the first std::string find in the file. */
long int LocaliseDataLogBefore(std::ifstream&, std::string, std::string);

    //! A function taking two arguments and returning a long int value.
    /*! \return The next position after the std::string search in the file. */
long int LocaliseNextDataLog(std::ifstream&, std::string);

    //! A function taking two arguments and returning a long int value.
    /*! \return The next position before the std::string search in the file. */
long int LocaliseNextDataLogBefore(std::ifstream&, std::string);

    //! A function taking two arguments and returning a long int value.
    /*! \return The position before the first std::string find in the file. */
long int LocaliseNextDataLogBefore(std::ifstream&, std::string, std::string);

#endif
