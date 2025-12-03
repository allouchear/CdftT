#ifndef CDFTT_FCHK_H_INCLUDED
#define CDFTT_FCHK_H_INCLUDED

#include<iostream>
#include<string>
#include<vector>
#include <Utils/Utils.h>

using std::ifstream;
using std::string;
using std::vector;

/**
 * @brief FCHK class.
 *
 * This class is used to read and store data from a formatted checkpoint (FCHK) file.
 */
class FCHK
{
    private:
        /** @brief Total number of electrons in the system. */
        int _number_of_electrons;

        /** @brief Number of alpha (spin-up) electrons. */
        int _number_of_alpha_electrons;

        /** @brief Number of beta (spin-down) electrons. */
        int _number_of_beta_electrons;

        /** @brief Total number of basis functions. */
        int _number_of_basis_functions;

        /** @brief Atomic numbers of all atoms in the system. */
        vector<int> _atomic_numbers;

        /** @brief Nuclear charges of all atoms in the system. */
        vector<double> _nuclear_charges;

        /** @brief Cartesian coordinates of all atoms. */
        vector<double> _current_cartesian_coordinates; // or vector<vector<double>> --> ??

        /** @brief Primitive exponents for contracted Gaussian-type functions (CGTFs). */
        vector<double> _primitive_exponents;

        /** @brief Contraction coefficients for CGTFs. */
        vector<double> _contraction_coefficients; // à voir --> ??

        /** @brief Alpha molecular orbital energies. */
        vector<double> _alpha_orbital_energies;

        /** @brief Alpha molecular orbital coefficients. */
        vector<double> _alpha_mo_coefficients; // or vector<vector<double>> --> ??

        /** @brief Beta molecular orbital energies. */
        vector<double> _beta_orbital_energies;

        /** @brief Beta molecular orbital coefficients. */
        vector<double> _beta_mo_coefficients; // or vector<vector<double>> --> ??

        /** @brief Total energy of the system. */
        double _total_energy;

        /** @brief Mulliken charges for all atoms. */
        vector<double> _mulliken_charges;

        /** @brief Natural Population Analysis (NPA) charges for all atoms. */
        vector<double> _npa_charges;

        /** @brief Total dipole moment components. */
        vector<double> _dipole_moment;

        /** @brief Shell types for all basis functions. */
        vector<int> _shell_types; // int ?? alors que l_type est string dans CGTF

        /** @brief Mapping of shells to atoms. */
        vector<int> _shell_to_atom_map;

        /** @brief Number of primitives per shell. */
        vector<int> _number_of_primitives_per_shell;

        /** @brief Contraction coefficients for sp shells. */
        vector<double> _sp_contraction_coefficients;

        /** @brief Coordinates for all shells. */
        vector<double> _coordinates_for_shells;

        /** @brief Alpha occupation numbers for molecular orbitals. */
        vector<double> _alpha_occupation;

        /** @brief Beta occupation numbers for molecular orbitals. */
        vector<double> _beta_occupation;

        /** @brief Total number of atoms in the system. */
        int _number_of_atoms;

        /** @brief Total number of contracted shells. */
        int _number_of_contracted_shells;

        /** @brief Total number of primitive shells. */
        int _number_of_primitive_shells;

        /** @brief Highest angular momentum among all shells. */
        int _highest_angular_momentum;

        /** @brief Internal flag for alpha orbital processing. */ // ??
        int _ok_alpha;

        /** @brief True if alpha and beta orbitals share the same coefficients. */ // ??
        bool _alpha_and_beta;

        /** @brief True if sp shells are present. */ // ??
        bool _sp;

        /** @brief True if mixed basis sets (Cartesian and spherical) are used. */
        bool _mixte;

    public:
        /**
         * @brief Default constructor.
         *
         * Initializes all attributes to default values (0 for numeric, empty for vectors, false for booleans).
         */
        FCHK();

        /**
         * @brief Constructor from file.
         *
         * Reads data from a formatted checkpoint (FCHK) file and initializes attributes.
         *
         * @param file Input file stream opened on an FCHK file.
         */
        FCHK(ifstream& file);

        /**
         * @brief Default destructor.
         *
         * Not used explicitly.
         */
        ~FCHK() {}

        /**
         * @brief Returns the total number of electrons in the system.
         */
        int NumberOfElectrons() { return _number_of_electrons; }

        /**
         * @brief Returns the number of alpha (spin-up) electrons.
         */
        int NumberOfAlphaElectrons() { return _number_of_alpha_electrons; }

        /**
         * @brief Returns the number of beta (spin-down) electrons.
         */
        int NumberOfBetaElectrons() { return _number_of_beta_electrons; }

        /**
         * @brief Returns the total number of basis functions.
         */
        int NumberOfBasisFunctions() { return _number_of_basis_functions; }

        /**
         * @brief Returns the atomic numbers of all atoms in the system.
         */
        vector<int> AtomicNumbers() { return _atomic_numbers; }

        /**
         * @brief Returns the nuclear charges of all atoms in the system.
         */
        vector<double> NuclearCharges() { return _nuclear_charges; }

        /**
         * @brief Returns the Cartesian coordinates of all atoms in the system.
         */
        vector<double> CurrentCartesianCoordinates() { return _current_cartesian_coordinates; }

        /**
         * @brief Returns the primitive exponents for contracted Gaussian-type functions (CGTFs).
         */
        vector<double> PrimitiveExponents() { return _primitive_exponents; }

        /**
         * @brief Returns the contraction coefficients for CGTFs.
         */
        vector<double> ContractionCoefficients() { return _contraction_coefficients; }

        /**
         * @brief Returns the alpha molecular orbital energies.
         */
        vector<double> AlphaOrbitalEnergies() { return _alpha_orbital_energies; }

        /**
         * @brief Returns the alpha molecular orbital coefficients.
         */
        vector<double> AlphaMOCoefficients() { return _alpha_mo_coefficients; }

        /**
         * @brief Returns the beta molecular orbital energies.
         */
        vector<double> BetaOrbitalEnergies() { return _beta_orbital_energies; }

        /**
         * @brief Returns the beta molecular orbital coefficients.
         */
        vector<double> BetaMOCoefficients() { return _beta_mo_coefficients; }

        /**
         * @brief Returns the total energy of the system.
         */
        double TotalEnergy() { return _total_energy; }

        /**
         * @brief Returns the Mulliken charges for all atoms.
         */
        vector<double> MullikenCharges() { return _mulliken_charges; }

        /**
         * @brief Returns the Natural Population Analysis (NPA) charges for all atoms.
         */
        vector<double> NpaCharges() { return _npa_charges; }

        /**
         * @brief Returns the dipole moment components.
         */
        vector<double> DipoleMoment() { return _dipole_moment; }

        /**
         * @brief Returns the shell types for all basis functions (S, P, D, F, etc.).
         */
        vector<int> ShellTypes() { return _shell_types; }

        /**
         * @brief Returns the mapping of shells to atoms.
         */
        vector<int> ShellToAtomMap() { return _shell_to_atom_map; }

        /**
         * @brief Returns the number of primitives per shell.
         */
        vector<int> NumberOfPrimitivesPerShell() { return _number_of_primitives_per_shell; }

        /**
         * @brief Returns the contraction coefficients for sp shells.
         */
        vector<double> spContractionCoefficients() { return _sp_contraction_coefficients; }

        /**
         * @brief Returns the coordinates for all shells.
         */
        vector<double> CoordinatesForShells() { return _coordinates_for_shells; }

        /**
         * @brief Returns the alpha occupation number for molecular orbitals.
         */
        vector<double> AlphaOccupation() { return _alpha_occupation; }

        /**
         * @brief Returns the beta occupation number for molecular orbitals.
         */
        vector<double> BetaOccupation() { return _beta_occupation; }

        /**
         * @brief Returns the number of primitive shells.
         */
        int NumberOfPrimitivesShells() { return _number_of_primitive_shells; }

        /**
         * @brief Returns the total number of atoms in the system.
         */
        int NumberOfAtoms() { return _number_of_atoms; }

        /**
         * @brief Returns the total number of contracted shells.
         */
        int NumberOfContractedShells() { return _number_of_contracted_shells; }

        /**
         * @brief Returns the highest angular momentum among all shells.
         */
        int HighestAngularMomentum() { return _highest_angular_momentum; }

        /**
         * @brief Returns true if alpha and beta orbitals share the same coefficients.
         */
        bool AlphaAndBeta() { return _alpha_and_beta; }

        /**
         * @brief Reads an integer value from the FCHK file.
         *
         * @param file Input file stream.
         * @param label Label to search for in the file.
         * @return The integer value read from the file.
         */
        int read_one_int(ifstream& file, string label);

        /**
         * @brief Reads a double value from the FCHK file.
         *
         * @param file Input file stream.
         * @param label Label to search for in the file.
         * @return The double value read from the file.
         */
        double read_one_real(ifstream& file, string label);

        /**
         * @brief Reads a block of integers from the FCHK file.
         *
         * @param file Input file stream.
         * @param label Label to search for in the file.
         * @return A vector of integers read from the file.
         */
        vector<int> read_one_block_int(ifstream& file, string label);

        /**
         * @brief Reads a block of doubles from the FCHK file.
         *
         * @param file Input file stream.
         * @param label Label to search for in the file.
         * @return A vector of doubles read from the file.
         */
        vector<double> read_one_block_real(ifstream& file, string label);

        /**
         * @brief Reads all data from the FCHK file.
         *
         * Initializes all attributes based on the data in the file.
         *
         * @param file Input file stream opened on an FCHK file.
         */
        void read_file_fchk(ifstream& file);

        /**
         * @brief Prints all stored data to the standard output.
         */
        void PrintData();

        /**
         * @brief Returns true if mixed basis sets are used.
         *
         * @return True if mixed basis sets (Cartesian and spherical) are used, false otherwise.
         */
        bool Mixte() { return _mixte; }
};

/**
 * @brief Searches for a label in the FCHK file and returns its position.
 *
 * @param file Input file stream.
 * @param label Label to search for in the file.
 * @return The position after the label in the file, or -1 if not found.
 */
long int LocaliseData(ifstream& file, string label);

#endif
