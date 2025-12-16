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
 *  Manages excited states and their associated transitions.
 */
class ExcitedState
{
    private:
        typedef std::pair<int, SpinType> OrbitalState;

        /** @brief Energy of the excited state. */
        double _energy;
        
        /** @brief Electronic transitions associated with the excited state. */
        std::vector<std::tuple<OrbitalState, OrbitalState, double>> _transitions;

        /** @brief Slater determinants associated with the excited state (one for each transition). */
        std::vector<SlaterDeterminant> _slaterDeterminants;
    

    public:
        /**
         * @brief Constructor.
         * 
         * @param energy Energy of the excited state, in Hartree.
         */
        ExcitedState(const double energy);

        /**
         * @brief Returns the excited state's energy, in Hartree.
         */
        double get_energy() const;

        /**
         * @brief Returns the Slater determinants associated with the excited state.
         */
        const std::vector<SlaterDeterminant>& get_slaterDeterminants() const;

        /**
         * @brief Returns the number of electronic transitions associated with the excited state.
         */
        int getNumberOfTransitions() const;

        /**
         * @brief Adds an electronic transition to the excited state.
         * 
         * @param initialOrbital Initial orbital state (number and spin).
         * @param finalOrbital Final orbital state (number and spin).
         * @param coefficient Coefficient of the transition.
         */
        void addTransition(const OrbitalState& initialOrbital, const OrbitalState& finalOrbital, const double coefficient);

        /**
         * @brief Computes the Slater determinant of the excited state.
         * 
         * @param[in]  SlaterDeterminant Reference to the ground state Slater determinant.
         */
        void computeSlaterDeterminants(const SlaterDeterminant& groundStateSlaterDeterminant);

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
        // STATIC FUNCTIONS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Reads a transitions file and populates a vector of ExcitedState objects.
         * 
         * @param[in] transitionsFileName Name of the transitions file to read.
         * @param[out] excitedStates Vector of ExcitedState objects populated from the file.
         */
        static void readTransitionsFile(std::string& transitionsFileName, std::vector<ExcitedState>& excitedStates);


        //----------------------------------------------------------------------------------------------------//
        // FRIEND FUNCTIONS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Overloads the output stream redirection operator for an ExcitedState.
         *
         * Prints energy of the excited state and transitions associated with it.
         *
         * @param stream Output stream.
         * @param excitedState ExcitedState to print.
         * @return Reference to the output stream.
         */
        friend std::ostream& operator<<(std::ostream& stream, const ExcitedState& excitedState);
};

#endif // CDFTT_EXCITEDSTATE_HPP_INCLUDED