#ifndef CDFTT_EXCITEDSTATE_HPP_INCLUDED
#define CDFTT_EXCITEDSTATE_HPP_INCLUDED

#include <array>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <Cube/Grid.h>
#include <Orbitals/Orbitals.h>
#include <Orbitals/SlaterDeterminant.hpp>
#include <Utils/Enums.hpp>


/** @brief ExcitedState class.
 * 
 *  Manages excited states and their associated electronic transitions.
 */
class ExcitedState
{
    private:
        typedef std::pair<int, SpinType> SpinOrbital;
        
        /** @brief Electronic transitions associated with the excited state. */
        std::vector<std::tuple<SpinOrbital, SpinOrbital, double>> _electronicTransitions;

        /** @brief Energy of the excited state. */
        double _energy;

        /** @brief Number of the excited state. */
        int _number;

        /** @brief Slater determinants associated with the excited state (one for each transition). */
        std::vector<std::pair<SlaterDeterminant, double>> _slaterDeterminants;

        /** @brief argsort of the coefficients of the Slater determinants */
        std::vector<size_t> _argsortCoefs;

        /** @brief Slater Determinants indices associated to each degree of excitation */
        std::vector<std::vector<size_t>> _excitationDegree;

        //----------------------------------------------------------------------------------------------------//
        // STATIC FIELDS
        //----------------------------------------------------------------------------------------------------//

        /** @brief Static ground state Slater determinant instance shared among all SlaterDeterminant objects. */
        static SlaterDeterminant _s_GS_SD;

    
        //----------------------------------------------------------------------------------------------------//
        // PRIVATE STATIC METHODS
        //----------------------------------------------------------------------------------------------------//
        
        /**
         * @brief computes the RDM-1 matrix using the Gamma method, between to states i,j (i=j for eletronc density, i!=j for transtition density)
         * 
         * @param[in] matrix to compute
         * @param[in] psi_i state i
         * @param[in] psi_j state_j
         * @param[in] orbitals Orbitals object containing the molecular orbitals information
         * @param[in] ignoredMos list of molecular orbitals to not take into account in calculation
         * @param[in] showProgress Whether to display a progress bar during matrix creation.
         */
        static void computeGammaMatrix(std::vector<std::vector<std::vector<double>>>& gammaMatrix, const ExcitedState& psi_i, const ExcitedState& psi_j, const Orbitals& orbitals, const std::vector<int>& ignoredMos, bool showProgress);

        /**
         * @brief computes the RDM-1 matrix using the X method, between to states i,j (i=j for eletronc density, i!=j for transtition density, for transition must have i=0)
         * 
         * @param[in] matrix to compute
         * @param[in] psi_i state i
         * @param[in] psi_j state_j
         * @param[in] orbitals Orbitals object containing the molecular orbitals information
         * @param[in] ignoredMos list of molecular orbitals to not take into account in calculation
         * @param[in] showProgress Whether to display a progress bar during matrix creation.
         */
        static void computeXMatrix(std::vector<std::vector<std::vector<double>>>& xMatrix, const ExcitedState& psi1, const ExcitedState& psi2, const Orbitals& orbitals, const std::vector<int>& ignoredMos, bool showProgress);
    

    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTORS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Constructor.
         *
         * @param[in] number Number of the excited state.
         * @param[in] energy Energy of the excited state, in Hartree.
         */
        ExcitedState(const int number, const double energy);

        /**
         * @brief Constructor for the ground state.
         *
         * @param[in] energy Energy of the ground state, in Hartree.
         * @param[in] slaterDeterminant Slater determinant associated with the ground state.
         */
        ExcitedState(const double energy, const SlaterDeterminant& slaterDeterminant);

        /**
         * @brief Constructor for the ground state.
         *
         * @param[in] energy Energy of the ground state, in Hartree.
         * @param[in] slaterDeterminant Slater determinant associated with the ground state.
         * @param[in] GS if true sets the static _s_GS_SD attribute to slaterDeterminant
         */
        ExcitedState(const double energy, const SlaterDeterminant& slaterDeterminant, bool GS);

        //----------------------------------------------------------------------------------------------------//
        // GETTERS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Returns the excited state's energy, in Hartree.
         */
        double get_energy() const;

        /**
         * @brief Returns the number of the excited state.
         */
        int get_number() const;

        /**
         * @brief Returns the Slater determinants associated with the excited state along with their respective coefficient.
         */
        const std::vector<std::pair<SlaterDeterminant, double>>& get_slaterDeterminants() const;

        /**
         * @brief Returns the argsort of the coefficients of the Slater determinants
         */
        std::vector<size_t> get_argsortCoefs() const;

        /**
         * @brief Returns the list of indices of the Slater determinants corresponding to each degree of excitation
         *        get_excitationDegree()[i] -> indices for every SD having i excitations
         */
        const std::vector<std::vector<size_t>>& get_excitationDegree() const;
        //----------------------------------------------------------------------------------------------------//
        // SETTERS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief sets _argsortCoefs, as a vector of indices giving the _slaterDeterminants of coefficients above a given treshold in descending order
         * 
         * @param[in] treshold limit on the coefficient of the Slater determinant, keep only Slater determinants above it
         */        
        void set_argsortCoefs(double treshold);

        /**
         * @brief set _excitationDegree, _excitationDegree[i] contains the list of all SD (their index in _slaterDeterminants) having i excitations (with respect to ground state SD)
         */
        void set_excitationDegree();

        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Adds a Slater determinant to the excited state along with its coefficient.
         * 
         * @param[in] slaterDeterminant Slater determinant to add.
         * @param[in] coefficient Coefficient associated with the Slater determinant.
         */
        void addSlaterDeterminant(const SlaterDeterminant& slaterDeterminant, const double coefficient);

        /**
         * @brief gives the degree of excitation of the SD, with respect to the static Ground State Determinant
         * 
         * @param[in] SD reference to the Slater determinant of which we want to get the degree of excitation
         * 
         * @return the degree of excitation of SD
         */
        size_t getExcitation(const SlaterDeterminant& SD) const;

        /**
         * @brief Adds an electronic transition to the excited state.
         *
         * @param[in] initialOrbital Initial orbital state (number and spin).
         * @param[in] finalOrbital Final orbital state (number and spin).
         * @param[in] coefficient Coefficient of the transition.
         */
        void addTransition(const SpinOrbital& initialOrbital, const SpinOrbital& finalOrbital, const double coefficient);

        /**
         * @brief Computes the Slater determinant of the excited state from the electronic transitions.
         *
         * @param[in]  SlaterDeterminant Reference to the ground state Slater determinant.
         */
        void computeSlaterDeterminants(const SlaterDeterminant& groundStateSlaterDeterminant);

        /**
         * @brief Returns the number of electronic transitions associated with the excited state.
         */
        int getNumberOfTransitions() const;

        /**
         * @brief Returns whether the state is a ground state or an excited state.
         */
        bool isGroundState() const;

        /**
         * @brief Computes the density grid for the excited state.
         * 
         * @param[in] orbitals Orbitals object containing the molecular orbitals information.
         * @param[out] grid Grid object on which to compute the density.
         * @param[in] showProgress Whether to display a progress bar during grid creation.
         */
        void makeDensityGrid(Orbitals& orbitals, Grid& grid, bool showProgress = false);

        /** @brief Prints lambda diagnostic for the excited state.
         *
         * Ref: Excitation energies in density functional theory: An evaluation and a diagnostic test
         * M. J. G. Peach, P. Benfield, T. Helgaker, D. J. Tozer
         * J. Chem. Phys. 128, 044118 (2008); DOI 10.1063/1.2831900
         *
         * @param grid Grid used to compute the diagnostic.
         */
        void printLambdaDiagnostic(const Grid& grid) const;

        //----------------------------------------------------------------------------------------------------//
        // STATIC METHODS
        //----------------------------------------------------------------------------------------------------//

        ////////////////////////////////////////
        // GROUND STATE ENERGY READING METHODS

        /**
         * @brief Reads the energy of the ground state from a file.
         * 
         * @param[in] fileName Name of the file to read.
         * @param[out] groundStateEnergy Energy of the ground state, in Hartree.
         * 
         * @return True if reading was successful, false otherwise.
         */
        static bool readGroundStateEnergy(const std::string& fileName, double& groundStateEnergy);

        /**
         * @brief Reads the energy of the ground state from an Orca .out file.
         * 
         * @param[in] orcaOutFileName Name of the Orca output file to read.
         * @param[out] energy Energy of the ground state, in Hartree.
         * 
         * @return True if reading was successful, false otherwise.
         */
        static bool readGroundStateEnergyFromOutFile(const std::string& orcaOutFileName, double& energy);

        /**
         * @brief Reads the energy of the ground state from a transitions file.
         * 
         * @param[in] transitionsFileName Name of the transitions file to read.
         * @param[out] groundStateEnergy Energy of the ground state, in Hartree.
         * 
         * @return True if reading was successful, false otherwise.
         */
        static bool readGroundStateEnergyFromTransitionsFile(const std::string& transitionsFileName, double& groundStateEnergy);

        ////////////////////////////////
        // TRANSITIONS READING METHODS

        /**
         * @brief Reads transitions from the provided file and populates a vector of ExcitedState objects.
         *
         * @param[in] fileName Name of the file to read.
         * @param[out] excitedStates Vector of ExcitedState objects populated from the file.
         * @param[in] groundStateEnergy Energy of the ground state, in Hartree.
         * @param[in] maxNumberOfExcitedStates Maximum number of excited states to read (if -1, all are read).
         * 
         * @return True if reading was successful, false otherwise.
         */
        static bool readTransitions(const std::string& fileName, std::vector<ExcitedState>& excitedStates, const double groundStateEnergy, const double maxNumberOfExcitedStates = -1, const std::vector<int>& statesNumbersToKeep = std::vector<int>());

        /**
         * @brief Reads a transitions file and populates a vector of ExcitedState objects.
         *
         * @param[in] transitionsFileName Name of the transitions file to read.
         * @param[out] excitedStates Vector of ExcitedState objects populated from the file.
         * @param[in] groundStateEnergy Energy of the ground state, in Hartree.
         * @param[in] maxNumberOfExcitedStates Maximum number of excited states to read (if -1, all are read).
         * 
         * @return True if reading was successful, false otherwise.
         */
        static bool readTransitionsFile(const std::string& transitionsFileName, std::vector<ExcitedState>& excitedStates, const double groundStateEnergy, const double maxNumberOfExcitedStates = -1, const std::vector<int>& statesNumbersToKeep = std::vector<int>());

        /**
         * @brief Reads transitions from an Orca .out file and populates a vector of ExcitedState objects. *
         * @param[in] orcaOutFileName Name of the Orca output file to read.
         * @param[out] excitedStates Vector of ExcitedState objects populated from the file.
         * @param[in] groundStateEnergy Energy of the ground state, in Hartree. 
         * @param[in] maxNumberOfExcitedStates Maximum number of excited states to read (if -1, all are read).
         * 
         * @return True if reading was successful, false otherwise.
         */
        static bool readTransitionsFromOutFile(const std::string& orcaOutFileName, std::vector<ExcitedState>& excitedStates, const double groundStateEnergy, const double maxNumberOfExcitedStates = -1, const std::vector<int>& statesNumbersToKeep = std::vector<int>());

        ///////////////////////////////
        // LOADING AND SAVING METHODS

        /**
         * @brief Loads excited states from a file.
         * 
         * @param[in] excitedStatesFileName Name of the transitions file to read.
         * @param[out] states Vector of ExcitedState objects populated from the file.
         * @param[out] slaterDeterminants Vector of SlaterDeterminant objects populated from the file.
         * @param[in] maxNumberOfExcitedStates Maximum number of excited states to read (if -1, all are read).
         * 
         * @return True if reading was successful, false otherwise.
         */
        static bool loadExcitedStatesFromFile(const std::string& statesFileName, std::vector<ExcitedState>& states, std::vector<SlaterDeterminant>& slaterDeterminants, int maxNumberOfExcitedStates = -1);

        /**
         * @brief Saves excited states to a file.
         * 
         * @param[in] statesFileName Name of the file to save the excited states to.
         * @param[in] states Vector of ExcitedState objects to save.
         */
        static bool saveExcitedStatesToFile(const std::string& statesFileName, const std::vector<ExcitedState>& states);

        
        /////////////////////////
        // OTHER STATIC METHODS

        /**
         * @brief Builds a vector of states after a perturbation, based on the unperturbed states and the eigenvalues and eigenvectors.
         *
         * The method assumes that the first unperturbed state is the ground state.
         *
         * @param[in] unperturbedStates Vector of unperturbed excited states.
         * @param[in] energies Vector of energy values (eigenvalues).
         * @param[in] eigenvectors Matrix of eigenvectors.
         * 
         * @return Vector of perturbed excited states.
         */
        static std::vector<ExcitedState> buildPerturbedStates(const std::vector<ExcitedState>& unperturbedStates, const std::vector<double>& energies, const std::vector<std::vector<double>>& eigenvectors);

        /**
         * @brief Calculates the ionic potential between two excited states.
         *
         * @param[in] psi_i First excited state.
         * @param[in] psi_j Second excited state.
         * @param[in] ionicMatrixes Ionic matrix < phi_i | V_ion/electrons | phi_j > (the first index corresponds to alpha spin, the second to beta spin).
         * 
         * @return The ionic potential matrix element < psi_i | V_ion/electrons | psi_j >.
         */
        static double ionicPotential(const ExcitedState& psi_i, const ExcitedState& psi_j, const std::vector<std::vector<std::vector<double>>>& ionicMatrixes);


        /** 
         * @brief calls the right matrix construction function based on the RDMethod chosen
         * 
         * @param[in] matrix to compute
         * @param[in] rdmMethod method to use to compute the matrix
         * @param[in] psi_i state i
         * @param[in] psi_j state_j
         * @param[in] orbitals Orbitals object containing the molecular orbitals information
         * @param[in] ignoredMos list of molecular orbitals to not take into account in calculation
         * @param[in] showProgress Whether to display a progress bar during matrix creation.
        */
        static void reducedDensityMatrix(std::vector<std::vector<std::vector<double>>>& rdmMatrix, RDMMethod rdmMethod, const ExcitedState& psi_i, const ExcitedState& psi_j, const Orbitals& orbitals, const std::vector<int>& ignoredMos, bool showProgress);


        //----------------------------------------------------------------------------------------------------//
        // OPERATOR OVERLOADS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Overloads the output stream redirection operator for an ExcitedState.
         *
         * Prints energy of the excited state and transitions associated with it.
         *
         * @param[in,out] stream Output stream.
         * @param[in] excitedState ExcitedState to print.
         * 
         * @return Reference to the output stream.
         */
        friend std::ostream& operator<<(std::ostream& stream, const ExcitedState& excitedState);
};

#endif // CDFTT_EXCITEDSTATE_HPP_INCLUDED