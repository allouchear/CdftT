#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

#include <Orbitals/Orbitals.h>
#include <Orbitals/SlaterDeterminant.hpp>
#include <Utils/Enums.hpp>

//----------------------------------------------------------------------------------------------------//
// STATIC FIELDS
//----------------------------------------------------------------------------------------------------//

bool SlaterDeterminant::_s_isOrbitalsSet_ = false;
Orbitals SlaterDeterminant::_s_orbitals_ = Orbitals();


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------//

SlaterDeterminant::SlaterDeterminant():
    _occupiedOrbitals(2, std::vector<std::pair<int, double>>())
{ }


SlaterDeterminant::SlaterDeterminant(const Orbitals& orbitals):
    _occupiedOrbitals(2, std::vector<std::pair<int, double>>())
{
    if (!_s_isOrbitalsSet_)
    {
        _s_orbitals_ = orbitals;
        _s_isOrbitalsSet_ = true;
    }

    // Get occupation numbers
    const std::vector<std::vector<double>>& occupationNumbers = orbitals.get_occupationNumber();

    // Get spin type int values
    const int ALPHA = static_cast<int>(SpinType::ALPHA);
    const int BETA = static_cast<int>(SpinType::BETA);

    // Populate the occupied orbitals and their occupation numbers based on the orbitals object
    for (int i = 0; i < orbitals.get_numberOfMo(); ++i)
    {
        // Alpha spin
        if (occupationNumbers[ALPHA][i] == 1)
        {
            _occupiedOrbitals[ALPHA].emplace_back(i + 1, 1.0);
        }
        else if (occupationNumbers[ALPHA][i] == 2) // Case where _alpha_and_beta = true
        {
            // Alpha spin
            _occupiedOrbitals[ALPHA].emplace_back(i + 1, 1.0);

            // Beta spin
            _occupiedOrbitals[BETA].emplace_back(i + 1, 1.0);
        }

        // Beta spin
        if (occupationNumbers[BETA][i] == 1)
        {
            _occupiedOrbitals[BETA].emplace_back(i + 1, 1.0);
        }
    }
}


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

bool SlaterDeterminant::updateFromTransition(int initialOrbitalNumber, SpinType initialSpin, int finalOrbitalNumber, SpinType finalSpin)
{
    // The transition will be tagged as valid if it has been done successfully, i.e. if the initial orbital was found in the occupied orbitals and replaced by the final orbital.
    // Otherwise (e.g., if the initial orbital was not found in the occupied orbitals), the transition is considered invalid.
    bool valid = false;

    // Get spin type int values
    const int ALPHA = static_cast<int>(SpinType::ALPHA);
    const int BETA = static_cast<int>(SpinType::BETA);

    // Update occupied orbitals based on the transition
    if (initialSpin == SpinType::ALPHA && finalSpin == SpinType::ALPHA)
    {
        for (size_t i = 0; i < _occupiedOrbitals[ALPHA].size(); ++i)
        {
            if (_occupiedOrbitals[ALPHA][i].first == initialOrbitalNumber)
            {
                _occupiedOrbitals[ALPHA][i].first = finalOrbitalNumber;
                valid = true;
                break;
            }
        }
    }
    else if (initialSpin == SpinType::BETA && finalSpin == SpinType::BETA)
    {
        for (size_t i = 0; i < _occupiedOrbitals[BETA].size(); ++i)
        {
            if (_occupiedOrbitals[BETA][i].first == initialOrbitalNumber)
            {
                _occupiedOrbitals[BETA][i].first = finalOrbitalNumber;
                valid = true;
                break;
            }
        }
    }
    else
    {
        std::stringstream errorMessage;
        errorMessage << "Error in SlaterDeterminant::updateFromTransition(): Mixed spin transitions are not supported." << std::endl;
        errorMessage << "For transition from " << initialOrbitalNumber << " " << to_char(initialSpin) << " to " << finalOrbitalNumber << " " << to_char(finalSpin) << '.' << std::endl;
        print_error(errorMessage.str());

        std::exit(1);
    }

    return valid;
}


//----------------------------------------------------------------------------------------------------//
// STATIC METHODS
//----------------------------------------------------------------------------------------------------//

std::vector<std::vector<std::pair<int, int>>> SlaterDeterminant::getDifferences(const SlaterDeterminant& d_i, const SlaterDeterminant& d_j)
{
    std::vector<std::vector<std::pair<int, int>>> differences(2, std::vector<std::pair<int, int>>());

    // Get spin type int values
    const int ALPHA = static_cast<int>(SpinType::ALPHA);
    const int BETA = static_cast<int>(SpinType::BETA);
    std::array<int, 2> spins = { ALPHA, BETA };

    // Check that the determinants have the same number of occupied orbitals
    if (d_i._occupiedOrbitals[ALPHA].size() != d_j._occupiedOrbitals[ALPHA].size()
        || d_i._occupiedOrbitals[BETA].size() != d_j._occupiedOrbitals[BETA].size())
    {
        std::stringstream errorMessage;
        errorMessage << "Error in SlaterDeterminant::getDifferences(): Slater determinants have different numbers of occupied orbitals." << std::endl;
        print_error(errorMessage.str());

        std::exit(1);
    }

    // Check for differences in occupied orbitals for each spin type
    for (int spin : spins)
    {
        for (size_t i = 0; i < d_i._occupiedOrbitals[spin].size(); ++i)
        {
            if (d_i._occupiedOrbitals[spin][i].first != d_j._occupiedOrbitals[spin][i].first)
            {
                differences[spin].emplace_back(d_i._occupiedOrbitals[spin][i].first, d_j._occupiedOrbitals[spin][i].first);
            }
        }
    }

    return differences;
}

bool SlaterDeterminant::equivalent(const SlaterDeterminant& d_i, const SlaterDeterminant& d_j)
{
    bool equivalent = true;

    // Get spin type int values
    const int ALPHA = static_cast<int>(SpinType::ALPHA);
    const int BETA = static_cast<int>(SpinType::BETA);
    std::array<int, 2> spins = { ALPHA, BETA };

    std::vector<std::vector<int>> occupiedOrbitalsNumbers_i(2, std::vector<int>(d_i._occupiedOrbitals[0].size()));
    std::vector<std::vector<int>> occupiedOrbitalsNumbers_j(2, std::vector<int>(d_j._occupiedOrbitals[0].size()));

    for (int spin : spins)
    {
        // Get and sort occupied orbitals numbers for first Slater determinant
        for (size_t i = 0; i < d_i._occupiedOrbitals[spin].size(); ++i)
        {
            occupiedOrbitalsNumbers_i[spin][i] = d_i._occupiedOrbitals[spin][i].first;
        }
        std::sort(occupiedOrbitalsNumbers_i[spin].begin(), occupiedOrbitalsNumbers_i[spin].end());

        // Get and sort occupied orbitals numbers for second Slater determinant
        for (size_t i = 0; i < d_j._occupiedOrbitals[spin].size(); ++i)
        {
            occupiedOrbitalsNumbers_j[spin][i] = d_j._occupiedOrbitals[spin][i].first;
        }
        std::sort(occupiedOrbitalsNumbers_j[spin].begin(), occupiedOrbitalsNumbers_j[spin].end());
        
        // Compare sorted occupied orbitals numbers for both Slater determinants
        if(occupiedOrbitalsNumbers_i[spin] != occupiedOrbitalsNumbers_j[spin])
        {
            equivalent = false;
            break;
        }
    }

    return equivalent;
}

double SlaterDeterminant::overlap(const SlaterDeterminant& d_i, const SlaterDeterminant& d_j)
{
    return d_i == d_j ? 1.0 : 0.0;
}

double SlaterDeterminant::ionicPotential(const SlaterDeterminant& d_i, const SlaterDeterminant& d_j, const std::vector<std::vector<std::vector<double>>>& ionicMatrix)
{
    double sum = 0.0;

    // Get spin type int values;
    const int ALPHA = static_cast<int>(SpinType::ALPHA);
    const int BETA = static_cast<int>(SpinType::BETA);
    std::array<int, 2> spins = {ALPHA, BETA};

    // Apply Slater-Condon rules
    if (d_i == d_j)
    {
        // Interaction with the electrons
        for (int spin : spins)
        {
            for (size_t i = 0; i < d_i._occupiedOrbitals[spin].size(); ++i)
            {
                int orbitalIndex = d_i._occupiedOrbitals[spin][i].first - 1;
                sum += ionicMatrix[spin][orbitalIndex][orbitalIndex];
            }
        }
    }
    else
    {
        std::vector<std::vector<std::pair<int, int>>> differences = getDifferences(d_i, d_j);

        // Check that there is only one difference in total
        if (differences[ALPHA].size() + differences[BETA].size() == 1)
        {
            // Get the spin for which there is a difference
            int spin = differences[ALPHA].size() == 1 ? ALPHA : BETA;

            int firstOrbitalIndex = differences[spin][0].first - 1;
            int secondOrbitalIndex = differences[spin][0].second - 1;
            
            sum += (firstOrbitalIndex > secondOrbitalIndex ? ionicMatrix[spin][firstOrbitalIndex][secondOrbitalIndex] : ionicMatrix[spin][secondOrbitalIndex][firstOrbitalIndex]);
        }
    }

    return sum;
}


//----------------------------------------------------------------------------------------------------//
// OPERATOR OVERLOADS
//----------------------------------------------------------------------------------------------------//

bool operator==(const SlaterDeterminant& lhs, const SlaterDeterminant& rhs)
{
    std::vector<std::vector<std::pair<int, int>>> differences = SlaterDeterminant::getDifferences(lhs, rhs);
    return (differences[0].empty() && differences[1].empty());
}

bool operator!=(const SlaterDeterminant& lhs, const SlaterDeterminant& rhs)
{
    return !(lhs == rhs);
}

std::ostream& operator<<(std::ostream& stream, const SlaterDeterminant& slaterDeterminant)
{
    // Get spin type int values
    const int ALPHA = static_cast<int>(SpinType::ALPHA);
    const int BETA = static_cast<int>(SpinType::BETA);

    for (size_t i = 0; i < slaterDeterminant._occupiedOrbitals[ALPHA].size(); ++i)
    {
        stream << slaterDeterminant._occupiedOrbitals[ALPHA][i].first << "A(" << slaterDeterminant._occupiedOrbitals[ALPHA][i].second << ") ";
    }

    stream << "; ";

    for (size_t i = 0; i < slaterDeterminant._occupiedOrbitals[BETA].size(); ++i)
    {
        stream << slaterDeterminant._occupiedOrbitals[BETA][i].first << "B(" << slaterDeterminant._occupiedOrbitals[BETA][i].second << ") ";
    }

    return stream;
}