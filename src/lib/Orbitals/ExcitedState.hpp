#ifndef CDFTT_EXCITEDSTATE_HPP_INCLUDED
#define CDFTT_EXCITEDSTATE_HPP_INCLUDED

#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <Orbitals/Orbitals.h>
#include <Utils/Enums.hpp>


/** @brief ExcitedState class.
 * 
 *  Manages excited states and their associated transitions.
 */
class ExcitedState
{
    private:
        typedef std::pair<int, SpinType> OrbitalState;
        
        double _energy;
        std::vector<std::tuple<OrbitalState, OrbitalState, double>> _transitions;
    

    public:
        ExcitedState(const double energy);

        double get_energy() const;

        int numberOfTransitions() const;

        void addTransition(const OrbitalState& initialOrbital, const OrbitalState& finalOrbital, const double coefficient);

        void printLambdaDiagnostic(const Grid& grid) const;


        //----------------------------------------------------------------------------------------------------//
        // STATIC FUNCTIONS
        //----------------------------------------------------------------------------------------------------//

        /** @brief Computes lambda diagnostic.
         *
         * Ref: Excitation energies in density functional theory: An evaluation and a diagnostic test
         * M. J. G. Peach, P. Benfield, T. Helgaker, D. J. Tozer
         * J. Chem. Phys. 128, 044118 (2008); DOI 10.1063/1.2831900
         *
         * @param excitedStates Vector of excited states.
         * @param orbitals Orbitals object containing molecular orbitals information.
         * @return Lambda diagnostic value.
         */
        //static double computeLambdaDiagnostic(const std::vector<ExcitedState>& excitedStates, Orbitals& orbitals);

        static void readTransitionsFile(std::string& transitionsFileName, std::vector<ExcitedState>& excitedStates);


        //----------------------------------------------------------------------------------------------------//
        // FRIEND FUNCTIONS
        //----------------------------------------------------------------------------------------------------//

        friend std::ostream& operator<<(std::ostream& stream, const ExcitedState& excitedState);
};

#endif // CDFTT_EXCITEDSTATE_HPP_INCLUDED