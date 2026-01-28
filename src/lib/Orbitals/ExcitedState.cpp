#include <cmath>
#include <iostream>
#include <fstream>
#include <regex>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#ifdef ENABLE_OMP
#include <omp.h>
#endif

#include <Common/Constants.h>
#include <Orbitals/ExcitedState.hpp>
#include <Utils/Enums.hpp>


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//----------------------------------------------------------------------------------------------------//

ExcitedState::ExcitedState(const double energy) :
    _energy(energy),
    _electronicTransitions(),
    _slaterDeterminants()
{ }

ExcitedState::ExcitedState(const double energy, const SlaterDeterminant& slaterDeterminant) :
    _energy(energy),
    _electronicTransitions(),
    _slaterDeterminants(1, slaterDeterminant)
{ }


//----------------------------------------------------------------------------------------------------//
// GETTERS
//----------------------------------------------------------------------------------------------------//

double ExcitedState::get_energy() const
{
    return _energy;
}


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

void ExcitedState::addTransition(const OrbitalState& initialOrbital, const OrbitalState& finalOrbital, const double coefficient)
{
    _electronicTransitions.push_back(std::make_tuple(initialOrbital, finalOrbital, coefficient));
}

void ExcitedState::computeSlaterDeterminants(const SlaterDeterminant& groundStateSlaterDeterminant)
{
    // Apply the transitions to the Slater determinant
    for (const auto& transition : _electronicTransitions)
    {
        // Copy ground state Slater determinant
        SlaterDeterminant slaterDeterminantTransition(groundStateSlaterDeterminant);

        // Unpack transition
        int initialOrbitalNumber = std::get<0>(transition).first;
        SpinType initialOrbitalSpin = std::get<0>(transition).second;
        int finalOrbitalNumber = std::get<1>(transition).first;
        SpinType finalOrbitalSpin = std::get<1>(transition).second;

        // Update Slater determinant based on the transition
        slaterDeterminantTransition.updateFromTransition(initialOrbitalNumber, initialOrbitalSpin, finalOrbitalNumber, finalOrbitalSpin);

        // Store Slater determinant
        _slaterDeterminants.push_back(slaterDeterminantTransition);
    }
}

int ExcitedState::getNumberOfTransitions() const
{
    return static_cast<int>(_electronicTransitions.size());
}

std::vector<std::pair<SlaterDeterminant, double>> ExcitedState::getSlaterDeterminantsAndCoefficients() const
{
    std::vector<std::pair<SlaterDeterminant, double>> slaterDeterminantsAndCoefficients;

    // Handle ground state case
    if (_electronicTransitions.empty())
    {
        slaterDeterminantsAndCoefficients.emplace_back(_slaterDeterminants[0], 1.0);
    }
    else // Excited state case
    {
        for (size_t i = 0; i < _slaterDeterminants.size(); ++i)
        {
            slaterDeterminantsAndCoefficients.emplace_back(_slaterDeterminants[i], std::get<2>(_electronicTransitions[i]));
        }
    }

    return slaterDeterminantsAndCoefficients;
}

void ExcitedState::printLambdaDiagnostic(const Grid& grid) const
{
    // Get grid's infinitesimal volume element
    double dv = grid.get_domain().get_dv();


    // Compute lambda
    double sum_lambdaNumerator = 0.0;
    double sum_lambdaDenominator = 0.0;

    for (const auto& transition : _electronicTransitions)
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
// STATIC METHODS
//----------------------------------------------------------------------------------------------------//

bool ExcitedState::readTransitionsFile(const std::string& transitionsFileName, std::vector<ExcitedState>& excitedStates, const double groundStateEnergy)
{
    bool ok = true;

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
                        ok = false;

                        std::stringstream errorMessage;
                        errorMessage << "Error: unknown energy unit \"" << energyUnit << "\" in transitions file " << transitionsFileName << '.' << std::endl;
                        errorMessage << "Please use eV or H as energy unit.";

                        print_error(errorMessage.str());

                        std::exit(1);
                    }

                    ExcitedState excitedState(energy + groundStateEnergy);

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
                                ok = false;

                                std::stringstream errorMessage;
                                errorMessage << "Error: could not read transition in transitions file " << transitionsFileName << '.' << std::endl;
                                errorMessage << "Please check the documentation for the format of the file.";

                                print_error(errorMessage.str());

                                std::exit(1);
                            }
                        }
                    } while (!transitionsFile.eof() && !line.empty());

                    // Check that at least one transition was read
                    if (excitedState.getNumberOfTransitions() > 0)
                    {
                        // Add excited state to the list
                        excitedStates.push_back(excitedState);
                    }
                    else
                    {
                        ok = false;

                        std::stringstream errorMessage;
                        errorMessage << "Error: no transition found for excited state with energy " << excitedState.get_energy() << " in transitions file " << transitionsFileName << '.' << std::endl;
                        errorMessage << "Please check the documentation for the format of the file.";

                        print_error(errorMessage.str());

                        std::exit(1);
                    }
                }
                else
                {
                    ok = false;

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
        ok = false;

        std::stringstream errorMessage;
        errorMessage << "Error: could not open transitions file " << transitionsFileName << '.' << std::endl;
        errorMessage << "Please check that the file exists and is readable.";

        print_error(errorMessage.str());

        std::exit(1);
    }

    return ok;
}

bool ExcitedState::readTransitionsFromLogFile(const std::string& logFileName, std::vector<ExcitedState>& excitedStates, const double groundStateEnergy)
{
    bool ok = true;

    std::ifstream logFile(logFileName);
    if (logFile)
    {
        std::string line;
        while (!logFile.eof())
        {
            // Read line
            std::getline(logFile, line);
            line = trim_whitespaces(line, true, true);

            if (line.empty())
            {
                continue;
            }
            else
            {
                // New excited state: read energy
                std::regex energyRegex("Excited State\\s+\\d+:.*\\s+(?:(-?\\d*\\.?\\d+) eV).*");
                std::smatch energyRegexMatch;
                if (std::regex_search(line, energyRegexMatch, energyRegex))
                {
                    double energy = std::stod(energyRegexMatch[1]) * Constants::EV_TO_HARTREE;

                    ExcitedState excitedState(energy + groundStateEnergy);

                    do
                    {
                        std::getline(logFile, line);
                        line = trim_whitespaces(line, true, true);

                        if (!line.empty())
                        {
                            // First, consider the case where the spins are specified
                            std::regex transitionRegexAlphaBeta("(\\d+)(A|B)\\s+->\\s+(\\d+)(A|B)\\s+(-?\\d*\\.?\\d+)");
                            std::smatch transitionRegexAlphaBetaMatch;
                            if (std::regex_search(line, transitionRegexAlphaBetaMatch, transitionRegexAlphaBeta))
                            {
                                std::pair<int, SpinType> initialOrbital;
                                std::pair<int, SpinType> finalOrbital;

                                initialOrbital.first = std::stoi(transitionRegexAlphaBetaMatch[1]);
                                initialOrbital.second = (transitionRegexAlphaBetaMatch[2] == "A" ? SpinType::ALPHA : SpinType::BETA);

                                finalOrbital.first = std::stoi(transitionRegexAlphaBetaMatch[3]);
                                finalOrbital.second = (transitionRegexAlphaBetaMatch[4] == "A" ? SpinType::ALPHA : SpinType::BETA);

                                double coefficient = std::stod(transitionRegexAlphaBetaMatch[5]);

                                // debug
                                std::cout << "Found alpha and beta transition: " << initialOrbital.first << to_char(initialOrbital.second)
                                          << " -> " << finalOrbital.first << to_char(finalOrbital.second)
                                          << " with coefficient " << coefficient << std::endl;

                                excitedState.addTransition(initialOrbital, finalOrbital, coefficient);
                            }
                            else
                            {
                                // Then, consider the case where spins are not specified: both alpha and beta transitions are assumed
                                std::regex transitionRegex("(\\d+)\\s+->\\s+(\\d+)\\s+(-?\\d*\\.?\\d+)");
                                std::smatch transitionRegexMatch;
                                if (std::regex_search(line, transitionRegexMatch, transitionRegex))
                                {
                                    std::pair<int, SpinType> initialOrbital_alpha;
                                    std::pair<int, SpinType> finalOrbital_alpha;
                                    std::pair<int, SpinType> initialOrbital_beta;
                                    std::pair<int, SpinType> finalOrbital_beta;

                                    // Add alpha transition
                                    initialOrbital_alpha.first = std::stoi(transitionRegexMatch[1]);
                                    initialOrbital_alpha.second = SpinType::ALPHA;

                                    finalOrbital_alpha.first = std::stoi(transitionRegexMatch[2]);
                                    finalOrbital_alpha.second = SpinType::ALPHA;

                                    double coefficient = std::stod(transitionRegexMatch[3]);

                                    excitedState.addTransition(initialOrbital_alpha, finalOrbital_alpha, coefficient);

                                    // Add beta transition
                                    initialOrbital_beta.first = initialOrbital_alpha.first;
                                    initialOrbital_beta.second = SpinType::BETA;

                                    finalOrbital_beta.first = finalOrbital_alpha.first;
                                    finalOrbital_beta.second = SpinType::BETA;

                                    excitedState.addTransition(initialOrbital_beta, finalOrbital_beta, coefficient);
                                }
                            }
                        }
                    } while (!logFile.eof() && !line.empty());

                    // Check that at least one transition was read
                    if (excitedState.getNumberOfTransitions() > 0)
                    {
                        // Add excited state to the list
                        excitedStates.push_back(excitedState);
                    }
                    else
                    {
                        ok = false;

                        std::stringstream errorMessage;
                        errorMessage << "Error: no transition found for excited state with energy " << excitedState.get_energy() << " in log file " << logFileName << '.' << std::endl;
                        errorMessage << "Please check the documentation for the format of the file.";

                        print_error(errorMessage.str());

                        std::exit(1);
                    }
                }
            }
        }
    }
    else
    {
        ok = false;

        std::stringstream errorMessage;
        errorMessage << "Error: could not open transitions file " << logFileName << '.' << std::endl;
        errorMessage << "Please check that the file exists and is readable.";

        print_error(errorMessage.str());

        std::exit(1);
    }

    return ok;
}

bool ExcitedState::readTransitions(const std::string& fileName, std::vector<ExcitedState>& excitedStates, const double groundStateEnergy)
{
    bool ok = true;

    if (fileName.substr(fileName.length() - 4) == ".log")
    {
        // Read as a .log file
        ok = readTransitionsFromLogFile(fileName, excitedStates, groundStateEnergy);
    }
    else
    {
        // Try to read as a transitions file
        ok = readTransitionsFile(fileName, excitedStates, groundStateEnergy);
    }

    return ok;
}

//----------------------------------------------------------------------------------------------------//
// OPERATOR OVERLOADS
//----------------------------------------------------------------------------------------------------//

std::ostream& operator<<(std::ostream& stream, const ExcitedState& excitedState)
{
    stream << "Excited State Energy: " << excitedState._energy << " Hartree." << std::endl;
    stream << "Transitions:" << std::endl;

    for (const auto& transition : excitedState._electronicTransitions)
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