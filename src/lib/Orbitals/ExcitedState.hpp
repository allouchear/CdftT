#ifndef CDFTT_EXCITEDSTATE_HPP_INCLUDED
#define CDFTT_EXCITEDSTATE_HPP_INCLUDED

#include <array>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

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
        typedef std::pair<int, SpinType> OrbitalState;

        /** @brief Energy of the excited state. */
        double _energy;

        /** @brief Indicates if the state is the ground state. */
        bool _isGroundState;
        
        /** @brief Electronic transitions associated with the excited state. */
        std::vector<std::tuple<OrbitalState, OrbitalState, double>> _electronicTransitions;

        /** @brief Slater determinants associated with the excited state (one for each transition). */
        std::vector<SlaterDeterminant> _slaterDeterminants;
    

    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTORS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Constructor.
         *
         * @param energy Energy of the excited state, in Hartree.
         */
        ExcitedState(const double energy);

        /**
         * @brief Constructor for the ground state.
         *
         * @param[in] energy Energy of the ground state, in Hartree.
         * @param[in] slaterDeterminant Slater determinant associated with the ground state.
         */
        ExcitedState(const double energy, const SlaterDeterminant& slaterDeterminant);

        //----------------------------------------------------------------------------------------------------//
        // GETTERS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Returns the excited state's energy, in Hartree.
         */
        double get_energy() const;

        /**
         * @brief Returns whether the state is the ground state.
         */
        bool isGroundState() const;

        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Adds an electronic transition to the excited state.
         *
         * @param[in] initialOrbital Initial orbital state (number and spin).
         * @param[in] finalOrbital Final orbital state (number and spin).
         * @param[in] coefficient Coefficient of the transition.
         */
        void addTransition(const OrbitalState& initialOrbital, const OrbitalState& finalOrbital, const double coefficient);

        /**
         * @brief Computes the Slater determinant of the excited state.
         *
         * @param[in]  SlaterDeterminant Reference to the ground state Slater determinant.
         */
        void computeSlaterDeterminants(const SlaterDeterminant& groundStateSlaterDeterminant);

        /**
         * @brief Returns the number of electronic transitions associated with the excited state.
         */
        int getNumberOfTransitions() const;

        /**
         * @brief Returns the Slater determinants associated with the excited state along with their respective coefficient.
         */
        std::vector<std::pair<SlaterDeterminant, double>> getSlaterDeterminantsAndCoefficients() const;

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

        /**
         * @brief Reads a transitions file and populates a vector of ExcitedState objects.
         *
         * @param[in] transitionsFileName Name of the transitions file to read.
         * @param[out] excitedStates Vector of ExcitedState objects populated from the file.
         * @param[in] groundStateEnergy Energy of the ground state, in Hartree.
         */
        static bool readTransitionsFile(const std::string& transitionsFileName, std::vector<ExcitedState>& excitedStates, const double groundStateEnergy);

        /**
         * @brief Reads transitions from a .log file and populates a vector of ExcitedState objects.
         *
         * @param[in] logFileName Name of the log file to read.
         * @param[out] excitedStates Vector of ExcitedState objects populated from the file.
         * @param[in] groundStateEnergy Energy of the ground state, in Hartree.
         * @return True if reading was successful, false otherwise.
         */
        static bool readTransitionsFromLogFile(const std::string& logFileName, std::vector<ExcitedState>& excitedStates, const double groundStateEnergy);

        /**
         * @brief Reads transitions from an Orca .out file and populates a vector of ExcitedState objects. *
         * @param[in] orcaOutFileName Name of the Orca output file to read.
         * @param[out] excitedStates Vector of ExcitedState objects populated from the file.
         * @param[in] groundStateEnergy Energy of the ground state, in Hartree. 
         * @param[in] alphaAndBeta Whether beta transitions are the same as alpha ones (true) or not. 
         * @return True if reading was successful, false otherwise.
         */
        static bool readTransitionsFromOrcaOutFile(const std::string& orcaOutFileName, std::vector<ExcitedState>& excitedStates, const double groundStateEnergy, const bool alphaAndBeta);

        /**
         * @brief Reads transitions from the provided file and populates a vector of ExcitedState objects.
         *
         * @param[in] fileName Name of the file to read.
         * @param[out] excitedStates Vector of ExcitedState objects populated from the file.
         * @param[in] groundStateEnergy Energy of the ground state, in Hartree.
         * @return True if reading was successful, false otherwise.
         */
        static bool readTransitions(const std::string& fileName, std::vector<ExcitedState>& excitedStates, const double groundStateEnergy, const bool alphaAndBeta = true);

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
         * @return Reference to the output stream.
         */
        friend std::ostream& operator<<(std::ostream& stream, const ExcitedState& excitedState);
};

#endif // CDFTT_EXCITEDSTATE_HPP_INCLUDED