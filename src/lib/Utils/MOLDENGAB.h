#ifndef CDFTT_MOLDENGAB_H_INCLUDED
#define CDFTT_MOLDENGAB_H_INCLUDED

#include <fstream>
#include <istream>
#include <string>
#include <vector>


class MOLDENGAB
{
    private:
        std::vector<std::string> _symbol;
        std::vector<int> _atomic_number;
        std::vector<std::vector<double>> _coord;
        std::vector<std::string> _shell_types;
        std::vector<int> _L_types;
        std::vector<double> _exposants;
        std::vector<int> _number_of_gtf;
        std::vector<int> _num_center;
        std::vector<double> _cgtf_coefs;
        std::vector<double> _factor_coefs;
        std::vector<double> _MO_energy;
        std::vector<std::vector<double>> _MO_coefs;
        std::vector<double> _occupation;
        std::vector<std::string> _spin_types;
        std::string _coord_type;
        std::string _basis_or_gto;
        int _number_of_atoms;
        int _number_of_MO_coefs;
        int _number_of_MO;
        bool _alpha_and_beta; //true de base
        std::vector<int> _n_at_basis;
    //Pour séparer alpha et beta
        std::vector<double> _alpha_occupation;
        std::vector<double> _beta_occupation;
        std::vector<std::vector<double>> _alpha_MO_coefs;
        std::vector<std::vector<double>> _beta_MO_coefs;
        std::vector<double> _alpha_energies;
        std::vector<double> _beta_energies;
    //Pour un fichier gab
        std::string _format;
        std::string _cart_sphe;

        bool _mixte;

        
    public:
            //! A default constructor.
            /*! This constructor is used to set all of the parameters on 0 or "None" value. */
        MOLDENGAB();

            //! A constructor taking one argument.
            /*! This constructor is used to set all of the parameters with the data in the file. */
        MOLDENGAB(std::ifstream& f);

            //! A default desctructor.
            /*! We don't use it. */
        ~MOLDENGAB() {}

            //! A normal member taking one argument and returning a void value.
            /*! Set all the atomic data on the value in the file. */
        void read_atom_data(std::ifstream& f);

            //! A normal member taking one argument and returning a void value.
            /*! Set all the basis data on the value in the file. */
        void read_basis_data(std::ifstream& f);

            //! A normal member taking one argument and returning a void value.
            /*! Set all the molecular orbitals data on the value in the file. */
        void read_MO_data(std::ifstream& f);

            //! A normal member taking one argument and returning a void value.
            /*! Read one basis data in the file. */
        void read_one_basis_data(std::istream& f);

            //! A normal member taking one argument and returning a void value.
            /*! Print all the data save in the object. */
        void PrintData();

            //! A normal member taking no arguments and returning a std::vector<std::string> value.
            /*! \return The table of symbol of each atoms. */
        std::vector<std::string> Symbol() {return _symbol;}

            //! A normal member taking no arguments and returning a std::vector<int> value.
            /*! \return The table of atomic number of each atoms. */
        std::vector<int> AtomicNumbers() {return _atomic_number;}

            //! A normal member taking no arguments and returning a std::vector<std::vector<double>> value.
            /*! \return The table of coordinates of each atoms. */
        std::vector<std::vector<double>> Coordinates() {return _coord;}

            //! A normal member taking no arguments and returning a std::vector<std::string> value.
            /*! \return The table of shell type of each shell. */
        std::vector<std::string> ShellTypes() {return _shell_types;}

            //! A normal member taking no arguments and returning a std::vector<int> value.
            /*! \return The table of L type of each shell. */
        std::vector<int> Ltypes() {return _L_types;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of exponents of each CGTF. */
        std::vector<double> Exposants() {return _exposants;}

            //! A normal member taking no arguments and returning a std::vector<int> value.
            /*! \return The table of number of GTF for each centers (atoms). */
        std::vector<int> NumberOfGtf() {return _number_of_gtf;}

            //! A normal member taking no arguments and returning a std::vector<int> value.
            /*! \return The table of num center (the size of this table is equal to the number of atoms). */
        std::vector<int> Numcenter() {return _num_center;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of CGTF coefficients. */
        std::vector<double> CgtfCoefficients() {return _cgtf_coefs;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of CGTF factor coefficients. */
        std::vector<double> FactorCoefficients() {return _factor_coefs;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of molecular orbitals energy. */
        std::vector<double> MO_Energy() {return _MO_energy;}

            //! A normal member taking no arguments and returning a std::vector<std::vector<double>> value.
            /*! \return The table of CGTF factor coefficients. */
        std::vector<std::vector<double>> MO_coefficients() {return _MO_coefs;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of occupation number for each molecular orbital. */
        std::vector<double> OccupationNumber() {return _occupation;}

            //! A normal member taking no arguments and returning a std::vector<std::string> value.
            /*! \return The table of spin type (Alpha or Beta) for each molecular orbital. */
        std::vector<std::string> SpinTypes() {return _spin_types;}

            //! A normal member taking no arguments and returning a std::string value.
            /*! \return The unis of coordinates (Angs or AU). */
        std::string CoordinatesType() {return _coord_type;}

            //! A normal member taking no arguments and returning a std::string value.
            /*! \return If their is Basis or GTO (because their is different in a .gab or a .molden). */
        std::string BasisOrGTO() {return _basis_or_gto;}

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

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of alpha occupation number for each molecular orbital. */
        std::vector<double> AlphaOccupation() {return _alpha_occupation;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of beta occupation number for each molecular orbital. */
        std::vector<double> BetaOccupation() {return _beta_occupation;}

            //! A normal member taking no arguments and returning a std::vector<std::vector<double>> value.
            /*! \return The table of alpha molecular orbitals coefficients for each molecular orbital. */
        std::vector<std::vector<double>> AlphaMOCoefs() {return _alpha_MO_coefs;}

            //! A normal member taking no arguments and returning a std::vector<std::vector<double>> value.
            /*! \return The table of beta molecular orbitals coefficients for each molecular orbital. */
        std::vector<std::vector<double>> BetaMOCoefs() {return _beta_MO_coefs;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of alpha molecular orbitals energy for each molecular orbital. */
        std::vector<double> AlphaEnergies() {return _alpha_energies;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of alpha molecular orbitals energy for each molecular orbital. */
        std::vector<double> BetaEnergies() {return _beta_energies;}

            //! A normal member taking no arguments and returning a std::string value.
            /*! \return The format of the file (gabedit or molden). */
        std::string Format() {return _format;}

            //! A normal member taking no arguments and returning a std::string value.
            /*! \return The type of coordinates for each shell */
        std::string CartOrSphe() {return _cart_sphe;}

            //! A normal member taking no arguments and returning a std::vector<int> value.
            /*! \return The table of number of gtf in each centers. */
        std::vector<int> NatBasis() {return _n_at_basis;}

            //! A normal member taking no arguments and returning a bool value.
            /*! \return If their is a mixte basis. (Exemple : D is cartisian and F is spherical). */
        bool Mixte() {return _mixte;}
};

    //! A function taking two arguments and returning a long int value.
    /*! \return The position after the std::string search in the file. */
long int LocaliseDataMolGab(std::ifstream& f, std::string b);

    //! A function taking three arguments and returning a long int value.
    /*! \return The position before the std::string search in the file. */
long int LocaliseDataMolGabBefore(std::ifstream& f, std::string b);

#endif
