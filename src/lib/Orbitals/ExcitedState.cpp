#include <cmath>
#include <iostream>
#include <fstream>
#include <regex>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <Common/Constants.h>
#include <Orbitals/ExcitedState.hpp>
#include <Utils/Enums.hpp>


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//----------------------------------------------------------------------------------------------------//

ExcitedState::ExcitedState(const double energy) :
    _energy(energy)
{ }


//----------------------------------------------------------------------------------------------------//
// GETTERS
//----------------------------------------------------------------------------------------------------//

double ExcitedState::get_energy() const
{
    return _energy;
}


//----------------------------------------------------------------------------------------------------//
// PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

int ExcitedState::numberOfTransitions() const
{
    return static_cast<int>(_transitions.size());
}

void ExcitedState::addTransition(const OrbitalState& initialOrbital, const OrbitalState& finalOrbital, const double coefficient)
{
    _transitions.push_back(std::make_tuple(initialOrbital, finalOrbital, coefficient));
}

void ExcitedState::printLambdaDiagnostic(const Grid& grid) const
{
    // Get grid's infinitesimal volume element
    double dv = grid.get_domain().get_dv();


    // Compute lambda
    double sum_lambdaNumerator = 0.0;
    double sum_lambdaDenominator = 0.0;

    for (const auto& transition : _transitions)
    {
        // Unpack transition
        int initialOrbitalNumber = std::get<0>(transition).first;
        int finalOrbitalNumber = std::get<1>(transition).first;
        double kappa = std::get<2>(transition);
        double kappaSquared = kappa * kappa;

        // Integrate over the grid
        double sum_phiInitialTimesPhiInitial = 0.0; // Phi_initial squared norm < phi_initial | phi_initial >
        double sum_phiFinalTimesPhiFinal = 0.0; // Phi_final squared norm < phi_final | phi_final >

        double sum_phiInitialTimesPhiFinal = 0.0; // Overlap integral < phi_initial | phi_final >
        double sum_phiInitialAbsTimesPhiFinalAbs = 0.0; // Integral < |phi_initial| | |phi_final| >

        #ifdef ENABLE_OMP
        #pragma omp parallel for reduction(+:sum_phiInitialTimesPhiInitial, sum_phiFinalTimesPhiFinal, sum_phiInitialTimesPhiFinal, sum_phiInitialAbsTimesPhiFinalAbs)
        #endif
        for (int i = 0; i < grid.get_domain().get_N1(); ++i)
        {
            for (int j = 0; j < grid.get_domain().get_N2(); ++j)
            {
                for (int k = 0; k < grid.get_domain().get_N3(); ++k)
                {
                    double phiInitial = grid.value(i, j, k, initialOrbitalNumber - 1);
                    double phiFinal = grid.value(i, j, k, finalOrbitalNumber - 1);

                    sum_phiInitialTimesPhiInitial += phiInitial * phiInitial;
                    sum_phiFinalTimesPhiFinal += phiFinal * phiFinal;

                    sum_phiInitialTimesPhiFinal += phiInitial * phiFinal;
                    sum_phiInitialAbsTimesPhiFinalAbs += std::abs(phiInitial * phiFinal);
                }
            }
        }

        sum_phiInitialTimesPhiInitial *= dv;
        sum_phiFinalTimesPhiFinal *= dv;
        sum_phiInitialTimesPhiFinal *= dv;
        sum_phiInitialAbsTimesPhiFinalAbs *= dv;


        // Accumulate lambda sums
        sum_lambdaNumerator += kappaSquared * sum_phiInitialAbsTimesPhiFinalAbs;
        sum_lambdaDenominator += kappaSquared;


        // Print integrals
        std::cout << "< " << initialOrbitalNumber << " | " << initialOrbitalNumber << " > = " << sum_phiInitialTimesPhiInitial << std::endl;
        std::cout << "< " << finalOrbitalNumber << " | " << finalOrbitalNumber << " > = " << sum_phiFinalTimesPhiFinal << std::endl;
        std::cout << "< " << initialOrbitalNumber << " | " << finalOrbitalNumber << " > = " << sum_phiInitialTimesPhiFinal << std::endl;
        std::cout << "< |" << initialOrbitalNumber << "| | |" << finalOrbitalNumber << "| > = " << sum_phiInitialAbsTimesPhiFinalAbs << std::endl << std::endl;
    }


    // Print lambda diagnostic
    std::cout << "Lambda = " << sum_lambdaNumerator / sum_lambdaDenominator << std::endl;
}


//----------------------------------------------------------------------------------------------------//
// STATIC FUNCTIONS
//----------------------------------------------------------------------------------------------------//

/*
double ExcitedState::computeLambdaDiagnostic(const std::vector<ExcitedState>& excitedStates, Orbitals& orbitals)
{
    double sum_kappaSquaredTimesOverlapIntegral = 0.0;
    double sum_kappaSquared = 0.0;

    for (const ExcitedState& excitedState : excitedStates)
    {
        for (const auto& transition : excitedState._transitions)
        {
            // Unpack transition
            const auto& initialOrbital = std::get<0>(transition);
            const auto& finalOrbital = std::get<1>(transition);
            double coefficient = std::get<2>(transition);

            double kappaSquared = coefficient * coefficient;
            sum_kappaSquared += kappaSquared;
            
            double O_ia;

            sum_kappaSquaredTimesOverlapIntegral += kappaSquared * O_ia;
        }
    }

    return sum_kappaSquaredTimesOverlapIntegral / sum_kappaSquared;
}
*/

void ExcitedState::readTransitionsFile(std::string& transitionsFileName, std::vector<ExcitedState>& excitedStates)
{
    std::ifstream transitionsFile(transitionsFileName);
    if (transitionsFile)
    {
        std::string line;
        while (!transitionsFile.eof())
        {
            // Read line
            std::getline(transitionsFile, line);
            line = trim_whitespaces(line, true, true);

            if (line.empty())
            {
                // Skip empty lines in the beginning of the file
                continue;
            }
            else if (line[0] == '#')
            {
                // Comment line: skip
                continue;
            }
            else
            {
                // New excited state: read energy
                std::regex energyRegex("(?:energy)\\s+(\\d*\\.?\\d+)\\s+(eV|H)", std::regex_constants::icase);
                std::smatch energyRegexMatch;
                if (std::regex_search(line, energyRegexMatch, energyRegex))
                {
                    double energy = std::stod(energyRegexMatch[1]);
                    std::string energyUnit = energyRegexMatch[2];

                    // For the unit, we only analyze the first letter (eV or H)
                    if (std::toupper(energyUnit[0]) == 'E')
                    {
                        energy *= Constants::EV_TO_HARTREE;
                    }
                    else if (std::toupper(energyUnit[0]) != 'H')
                    {
                        std::stringstream errorMessage;
                        errorMessage << "Error: unknown energy unit \"" << energyUnit << "\" in transitions file " << transitionsFileName << '.' << std::endl;
                        errorMessage << "Please use eV or H as energy unit.";

                        print_error(errorMessage.str());

                        std::exit(1);
                    }

                    ExcitedState excitedState(energy);

                    // Read transitions
                    do
                    {
                        std::getline(transitionsFile, line);
                        line = trim_whitespaces(line, true, true);

                        if (line[0] == '#')
                        {
                            // Skip comment lines
                            continue;
                        }
                        else if (!line.empty())
                        {
                            std::regex transitionRegex("(\\d+)\\s+([aAbB])\\s+(\\d+)\\s+([aAbB])\\s+(-?\\d*\\.?\\d+)");
                            std::smatch transitionRegexMatch;
                            if (std::regex_search(line, transitionRegexMatch, transitionRegex))
                            {
                                std::pair<int, SpinType> initialOrbital;
                                std::pair<int, SpinType> finalOrbital;

                                initialOrbital.first = std::stoi(transitionRegexMatch[1]);
                                initialOrbital.second = (transitionRegexMatch[2] == "a" || transitionRegexMatch[2] == "A") ? SpinType::ALPHA : SpinType::BETA;

                                finalOrbital.first = std::stoi(transitionRegexMatch[3]);
                                finalOrbital.second = (transitionRegexMatch[4] == "a" || transitionRegexMatch[4] == "A") ? SpinType::ALPHA : SpinType::BETA;

                                double coefficient = std::stod(transitionRegexMatch[5]);

                                excitedState.addTransition(initialOrbital, finalOrbital, coefficient);
                            }
                            else
                            {
                                std::stringstream errorMessage;
                                errorMessage << "Error: could not read transition in transitions file " << transitionsFileName << '.' << std::endl;
                                errorMessage << "Please check the documentation for the format of the file.";

                                print_error(errorMessage.str());

                                std::exit(1);
                            }
                        }
                    } while (!transitionsFile.eof() && !line.empty());

                    // Check that at least one transition was read
                    if (excitedState.numberOfTransitions() > 0)
                    {
                        // Add excited state to the list
                        excitedStates.push_back(excitedState);
                    }
                    else
                    {
                        std::stringstream errorMessage;
                        errorMessage << "Error: no transition found for excited state with energy " << excitedState.get_energy() << " in transitions file " << transitionsFileName << '.' << std::endl;
                        errorMessage << "Please check the documentation for the format of the file.";

                        print_error(errorMessage.str());

                        std::exit(1);
                    }
                }
                else
                {
                    std::stringstream errorMessage;
                    errorMessage << "Error: could not read excited state energy in transitions file " << transitionsFileName << '.' << std::endl;
                    errorMessage << "Please check the documentation for the format of the file.";

                    print_error(errorMessage.str());

                    std::exit(1);
                }
            }
        }
    }
    else
    {
        std::stringstream errorMessage;
        errorMessage << "Error: could not open transitions file " << transitionsFileName << '.' << std::endl;
        errorMessage << "Please check that the file exists and is readable.";

        print_error(errorMessage.str());

        std::exit(1);
    }
}

//----------------------------------------------------------------------------------------------------//
// FRIEND FUNCTIONS
//----------------------------------------------------------------------------------------------------//

std::ostream& operator<<(std::ostream& stream, const ExcitedState& excitedState)
{
    stream << "Excited State Energy: " << excitedState._energy << " Hartree." << std::endl;
    stream << "Transitions:" << std::endl;

    for (const auto& transition : excitedState._transitions)
    {
        const ExcitedState::OrbitalState& initialOrbital = std::get<0>(transition);
        const ExcitedState::OrbitalState& finalOrbital = std::get<1>(transition);
        const double& coefficient =  std::get<2>(transition);

        stream << "  From Orbital " << initialOrbital.first << to_char(initialOrbital.second)
               << " to Orbital " << finalOrbital.first << to_char(finalOrbital.second)
               << " with Coefficient: " << coefficient << std::endl;
    }

    return stream;
}