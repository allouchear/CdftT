#include <cmath>
#include <cstdlib>
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
#include <Cube/Grid.h>
#include <Orbitals/ExcitedState.hpp>
#include <Utils/Enums.hpp>
#include <Utils/LOG.h>
#include <Utils/Utils.h>


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//----------------------------------------------------------------------------------------------------//

ExcitedState::ExcitedState(const double energy) :
    _energy(energy),
    _isGroundState(false),
    _electronicTransitions(),
    _slaterDeterminants()
{ }

ExcitedState::ExcitedState(const double energy, const SlaterDeterminant& slaterDeterminant) :
    _energy(energy),
    _isGroundState(true),
    _electronicTransitions(),
    _slaterDeterminants({ { slaterDeterminant, 1.0 } })
{ }


//----------------------------------------------------------------------------------------------------//
// GETTERS
//----------------------------------------------------------------------------------------------------//

double ExcitedState::get_energy() const
{
    return _energy;
}

const std::vector<std::pair<SlaterDeterminant, double>>& ExcitedState::get_slaterDeterminants() const
{
    return _slaterDeterminants;
}

bool ExcitedState::isGroundState() const
{
    return _isGroundState;
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
        double coefficient = std::get<2>(transition);

        // Update Slater determinant based on the transition
        slaterDeterminantTransition.updateFromTransition(initialOrbitalNumber, initialOrbitalSpin, finalOrbitalNumber, finalOrbitalSpin);

        // Store Slater determinant
        _slaterDeterminants.emplace_back(slaterDeterminantTransition, coefficient);
    }
}

int ExcitedState::getNumberOfTransitions() const
{
    return static_cast<int>(_electronicTransitions.size());
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

std::vector<ExcitedState> ExcitedState::buildPerturbedStates(const std::vector<ExcitedState>& unperturbedStates, const std::vector<double>& energies, const std::vector<std::vector<double>>& eigenvectors)
{
    const SlaterDeterminant& groundStateSlaterDeterminant = unperturbedStates[0].get_slaterDeterminants()[0].first;

    std::vector<ExcitedState> perturbedStates;
    perturbedStates.reserve(unperturbedStates.size());

    for (size_t i = 0; i < unperturbedStates.size(); ++i)
    {
        perturbedStates.emplace_back(energies[i]);
        ExcitedState& currentPerturbedState = perturbedStates.back();

        for (size_t k = 0; k < unperturbedStates.size(); ++k)
        {
            double c_k = eigenvectors[k][i];

            if (c_k != 0.0)
            {
                // First handle the case where the unperturbed state is the ground state (i.e., has no electronic transitions)
                if (unperturbedStates[k].isGroundState())
                {
                    // Add the ground state Slater determinant with its contribution
                    currentPerturbedState._slaterDeterminants.emplace_back(unperturbedStates[k]._slaterDeterminants[0].first, c_k);
                }
                else
                {
                    // Then handle the case where the unperturbed state is an excited state (i.e., has one or more electronic transitions)
                    for (const auto& electronicTransition : unperturbedStates[k]._electronicTransitions)
                    {
                        // Search for the electronic transition in the _electronicTransitions vector
                        auto it = std::find_if(currentPerturbedState._electronicTransitions.begin(),
                                            currentPerturbedState._electronicTransitions.end(),
                                            [&electronicTransition](const auto& element) { return (std::get<0>(element) == std::get<0>(electronicTransition) && std::get<1>(element) == std::get<1>(electronicTransition)); });

                        // If it is not found, add it to the electronic transitions with its contribution.
                        if (it == currentPerturbedState._electronicTransitions.end())
                        {
                            currentPerturbedState._electronicTransitions.emplace_back(std::get<0>(electronicTransition), std::get<1>(electronicTransition), c_k * std::get<2>(electronicTransition));
                        }
                        else
                        {
                            std::get<2>(*it) += c_k * std::get<2>(electronicTransition);
                        }
                    }
                }
            }
        }

        // Compute the Slater determinants associated with the perturbed state based on its electronic transitions
        currentPerturbedState.computeSlaterDeterminants(groundStateSlaterDeterminant);
    }

    return perturbedStates;
}

bool ExcitedState::readGroundStateEnergyFromOutFile(const std::string& orcaOutFileName, double& energy)
{
    bool ok = true;
    bool found = false;

    std::ifstream orcaOutFile(orcaOutFileName);
    if (orcaOutFile)
    {
        std::string line;
        while (!orcaOutFile.eof() && !found)
        {
            std::getline(orcaOutFile, line);
            line = trim_whitespaces(line, true, true);

            if (line.empty())
            {
                continue;
            }

            // Look for the ground state energy in the file
            std::regex energyRegex("Total Energy\\s+:\\s+(-?\\d*\\.?\\d+)\\s+Eh");
            std::smatch energyRegexMatch;
            if (std::regex_search(line, energyRegexMatch, energyRegex))
            {
                energy = std::stod(energyRegexMatch[1]);
                found = true;
            }
        }

        if (!found)
        {
            print_error("Error: could not read energy from OUT file.");
            std::exit(1);
        }
    }
    else
    {
        print_error("Error: could not read file " + orcaOutFileName + ".");
        std::exit(1);
    }

    return (ok && found);
}

bool ExcitedState::readGroundStateEnergyFromTransitionsFile(const std::string& transitionsFileName, double& groundStateEnergy)
{
    bool ok = true;
    bool found = false;

    std::ifstream transitionsFile(transitionsFileName);
    if (transitionsFile)
    {
        std::string line;
        while (!transitionsFile.eof() && !found)
        {
            std::getline(transitionsFile, line);
            line = trim_whitespaces(line, true, true);

            if (line.empty())
            {
                continue;
            }

            // Look for the ground state energy in the file
            std::regex groundStateEnergyRegex("Ground State Energy\\s+:\\s+(-?\\d*\\.?\\d+)\\s+(eV|H)");
            std::smatch groundStateEnergyRegexMatch;
            if (std::regex_search(line, groundStateEnergyRegexMatch, groundStateEnergyRegex))
            {
                groundStateEnergy = std::stod(groundStateEnergyRegexMatch[1]);
                found = true;
            }
        }
    }
    else
    {
        print_error("Error: could not read file " + transitionsFileName + ".");
        std::exit(1);
    }

    return (ok && found);
}

bool ExcitedState::readGroundStateEnergy(const std::string& fileName, double& groundStateEnergy)
{
    bool ok = true;

    if (to_lower(fileName.substr(fileName.length() - 5)) == ".fchk")
    {
        ok = FCHK::readGroundStateEnergy(fileName, groundStateEnergy);
    }
    else if (to_lower(fileName.substr(fileName.length() - 4)) == ".log")
    {
        ok = LOG::readGroundStateEnergy(fileName, groundStateEnergy);
    }
    else if (to_lower(fileName.substr(fileName.length() - 4)) == ".out")
    {
        ok = readGroundStateEnergyFromOutFile(fileName, groundStateEnergy);
    }
    else
    {
        // Try to read as a transitions file
        ok = readGroundStateEnergyFromTransitionsFile(fileName, groundStateEnergy);
    }
    
    return ok;
}

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
                std::regex energyRegex("(?:energy)\\s+(-?\\d*\\.?\\d+)\\s+(eV|H)", std::regex_constants::icase);
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
                            // First, consider the case where the spins are specified
                            std::regex transitionRegexAlphaBeta("(\\d+)\\s+([aAbB])\\s+(\\d+)\\s+([aAbB])\\s+(-?\\d*\\.?\\d+)");
                            std::smatch transitionRegexAlphaBetaMatch;
                            if (std::regex_search(line, transitionRegexAlphaBetaMatch, transitionRegexAlphaBeta))
                            {
                                std::pair<int, SpinType> initialOrbital;
                                std::pair<int, SpinType> finalOrbital;

                                initialOrbital.first = std::stoi(transitionRegexAlphaBetaMatch[1]);
                                initialOrbital.second = (transitionRegexAlphaBetaMatch[2] == "a" || transitionRegexAlphaBetaMatch[2] == "A") ? SpinType::ALPHA : SpinType::BETA;

                                finalOrbital.first = std::stoi(transitionRegexAlphaBetaMatch[3]);
                                finalOrbital.second = (transitionRegexAlphaBetaMatch[4] == "a" || transitionRegexAlphaBetaMatch[4] == "A") ? SpinType::ALPHA : SpinType::BETA;
                                double coefficient = std::stod(transitionRegexAlphaBetaMatch[5]);

                                excitedState.addTransition(initialOrbital, finalOrbital, coefficient);
                            }
                            else
                            {
                                // Then, consider the case where spins are not specified: both alpha and beta transitions are assumed
                                std::regex transitionRegex("(\\d+)\\s+(\\d+)\\s+(-?\\d*\\.?\\d+)");
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

        transitionsFile.close();
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

bool ExcitedState::readTransitionsFromOutFile(const std::string& orcaOutFileName, std::vector<ExcitedState>& excitedStates, const double groundStateEnergy)
{
    bool ok = true;
    bool hfTypeFound = false;
    HFType hfType = HFType::UNKNOWN;

    std::ifstream orcaOutFile(orcaOutFileName);
    if (orcaOutFile)
    {
        std::string line;
        while (!orcaOutFile.eof())
        {
            // Read line
            std::getline(orcaOutFile, line);
            line = trim_whitespaces(line, true, true);

            if (line.empty() || line[0] == '#')
            {
                continue;
            }
            else if (!hfTypeFound)
            {
                std::regex hfTypeRegex("Hartree-Fock type\\s+HFTyp\\s+\\.*\\s+(RHF|UHF)", std::regex_constants::icase);
                std::smatch hfTypeRegexMatch;
                if (std::regex_search(line, hfTypeRegexMatch, hfTypeRegex))
                {
                    hfTypeFound = true;
                    hfType = hfType_from_string(hfTypeRegexMatch[1]);
                }
            }
            else
            {
                // New excited state: read energy
                std::regex energyRegex("STATE\\s+\\d+:.*\\s+(?:(-?\\d*\\.?\\d+) au).*");
                std::smatch energyRegexMatch;
                if (std::regex_search(line, energyRegexMatch, energyRegex))
                {
                    double energy = std::stod(energyRegexMatch[1]);

                    ExcitedState excitedState(energy + groundStateEnergy);

                    do
                    {
                        std::getline(orcaOutFile, line);
                        line = trim_whitespaces(line, true, true);

                        if (!line.empty())
                        {
                            // First, consider the case where the spins are specified
                            std::regex transitionRegexAlphaBeta("(\\d+)(a|b)\\s+->\\s+(\\d+)(a|b).+c=\\s+(-?\\d*\\.?\\d+)");
                            std::smatch transitionRegexAlphaBetaMatch;
                            if (std::regex_search(line, transitionRegexAlphaBetaMatch, transitionRegexAlphaBeta))
                            {
                                std::pair<int, SpinType> initialOrbital;
                                std::pair<int, SpinType> finalOrbital;

                                initialOrbital.first = std::stoi(transitionRegexAlphaBetaMatch[1]) + 1; // +1 because Orca orbitals are 0-indexed, while we use 1-indexing for orbitals here
                                initialOrbital.second = (transitionRegexAlphaBetaMatch[2] == "a" ? SpinType::ALPHA : SpinType::BETA);

                                finalOrbital.first = std::stoi(transitionRegexAlphaBetaMatch[3]) + 1; // +1 because Orca orbitals are 0-indexed, while we use 1-indexing for orbitals here
                                finalOrbital.second = (transitionRegexAlphaBetaMatch[4] == "a" ? SpinType::ALPHA : SpinType::BETA);

                                double coefficient = std::stod(transitionRegexAlphaBetaMatch[5]);

                                // For RHF, we assume that alpha and beta transitions have the same coefficient
                                if (hfType == HFType::RHF)
                                {
                                    coefficient /= std::sqrt(2.0);
                                }

                                excitedState.addTransition(initialOrbital, finalOrbital, coefficient);

                                // For RHF, also add the beta transition
                                if (hfType == HFType::RHF)
                                {
                                    initialOrbital.second = SpinType::BETA;
                                    finalOrbital.second = SpinType::BETA;

                                    excitedState.addTransition(initialOrbital, finalOrbital, coefficient);
                                }
                            }
                        }
                    } while (!orcaOutFile.eof() && !line.empty());

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
                        errorMessage << "Error: no transition found for excited state with energy " << excitedState.get_energy() << " in Orca .out file " << orcaOutFileName << '.' << std::endl;
                        errorMessage << "Please check the documentation for the format of the file.";

                        print_error(errorMessage.str());

                        std::exit(1);
                    }
                }
            }
        }

        orcaOutFile.close();
    }
    else
    {
        ok = false;

        std::stringstream errorMessage;
        errorMessage << "Error: could not open transitions file " << orcaOutFileName << '.' << std::endl;
        errorMessage << "Please check that the file exists and is readable.";

        print_error(errorMessage.str());

        std::exit(1);
    }

    return ok;
}

bool ExcitedState::readTransitions(const std::string& fileName, std::vector<ExcitedState>& excitedStates, const double groundStateEnergy)
{
    bool ok = true;

    if (to_lower(fileName.substr(fileName.length() - 4)) == ".log")
    {
        ok = LOG::readTransitions(fileName, excitedStates, groundStateEnergy);
    }
    else if (to_lower(fileName.substr(fileName.length() - 4)) == ".out")
    {
        ok = readTransitionsFromOutFile(fileName, excitedStates, groundStateEnergy);
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
    stream << (excitedState._isGroundState ? "Ground" : "Excited") << " State Energy: " << excitedState._energy << " Hartree." << std::endl;

    if (!excitedState._isGroundState)
    {
        stream << "  Electronic transitions:" << std::endl;

        for (const auto& transition : excitedState._electronicTransitions)
        {
            const ExcitedState::OrbitalState& initialOrbital = std::get<0>(transition);
            const ExcitedState::OrbitalState& finalOrbital = std::get<1>(transition);
            const double& coefficient = std::get<2>(transition);

            stream << "    " << initialOrbital.first << to_char(initialOrbital.second)
                             << " -> " << finalOrbital.first << to_char(finalOrbital.second)
                             << ", coefficient: " << coefficient << std::endl;
        }

        stream << "  Slater determinants:" << std::endl;
        for (const auto& slaterCoeff : excitedState._slaterDeterminants)
        {
            stream << "    " << slaterCoeff.first << "; Coefficient: " << slaterCoeff.second << std::endl;
        }
    }

    return stream;
}