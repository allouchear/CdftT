#include <algorithm>
#include <cstdlib>
#include <vector>

#include <Orbitals/Orbitals.h>
#include <Orbitals/SlaterDeterminant.hpp>

#include <iostream>

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

void SlaterDeterminant::updateFromTransition(int initialOrbitalNumber, SpinType initialSpin, int finalOrbitalNumber, SpinType finalSpin)
{
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

        exit(1);
    }
}


//----------------------------------------------------------------------------------------------------//
// STATIC METHODS
//----------------------------------------------------------------------------------------------------//

std::vector<std::vector<std::pair<int, int>>> SlaterDeterminant::getDifferences(const SlaterDeterminant& di, const SlaterDeterminant& dj)
{
    std::vector<std::vector<std::pair<int, int>>> differences(2, std::vector<std::pair<int, int>>());

    // Get spin type int values
    const int ALPHA = static_cast<int>(SpinType::ALPHA);
    const int BETA = static_cast<int>(SpinType::BETA);
    std::array<int, 2> spins = { ALPHA, BETA };

    // Check that the determinants have the same number of occupied orbitals
    if (di._occupiedOrbitals[ALPHA].size() != dj._occupiedOrbitals[ALPHA].size()
        || di._occupiedOrbitals[BETA].size() != dj._occupiedOrbitals[BETA].size())
    {
        std::stringstream errorMessage;
        errorMessage << "Error in SlaterDeterminant::getDifferences(): Slater determinants have different numbers of occupied orbitals." << std::endl;
        print_error(errorMessage.str());

        exit(1);
    }

    // Check for differences in occupied orbitals for each spin type
    for (int spin : spins)
    {
        for (size_t i = 0; i < di._occupiedOrbitals[spin].size(); ++i)
        {
            if (di._occupiedOrbitals[spin][i].first != dj._occupiedOrbitals[spin][i].first)
            {
                differences[spin].emplace_back(di._occupiedOrbitals[spin][i].first, dj._occupiedOrbitals[spin][i].first);
            }
        }
    }

    return differences;
}

double SlaterDeterminant::ionicPotential(const SlaterDeterminant& di, const SlaterDeterminant& dj, const std::vector<std::vector<std::vector<double>>>& ionicMatrix)
{
    double sum = 0.0;

    // Get spin type int values;
    std::array<int, 2> spins = {static_cast<int>(SpinType::ALPHA), static_cast<int>(SpinType::BETA)};

    // Apply Slater-Condon rules
    if (di == dj)
    {
        for (int spin : spins)
        {
            for (size_t i = 0; i < di._occupiedOrbitals[spin].size(); ++i)
            {
                int orbitalIndex = di._occupiedOrbitals[spin][i].first - 1;
                sum += ionicMatrix[spin][orbitalIndex][orbitalIndex];
            }
        }
    }
    else
    {
        std::vector<std::vector<std::pair<int, int>>> differences = getDifferences(di, dj);

        for (int spin : spins)
        {
            // Check that there is only one difference
            if (differences[spin].size() == 1)
            {
                int initialOrbitalIndex = differences[spin][0].first - 1;
                int finalOrbitalIndex = differences[spin][0].second - 1;

                sum += ionicMatrix[spin][initialOrbitalIndex][finalOrbitalIndex];
            }
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

    stream << std::endl;

    for (size_t i = 0; i < slaterDeterminant._occupiedOrbitals[BETA].size(); ++i)
    {
        stream << slaterDeterminant._occupiedOrbitals[BETA][i].first << "B(" << slaterDeterminant._occupiedOrbitals[BETA][i].second << ") ";
    }

    return stream;
}