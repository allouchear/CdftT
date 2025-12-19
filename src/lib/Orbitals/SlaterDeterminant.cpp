#include <algorithm>
#include <cstdlib>
#include <vector>

#include <Orbitals/Orbitals.h>
#include <Orbitals/SlaterDeterminant.hpp>

#include <iostream>

//----------------------------------------------------------------------------------------------------//
// STATIC FIELDS
//----------------------------------------------------------------------------------------------------//

std::vector<std::vector<std::vector<double>>> SlaterDeterminant::_s_ionicMatrix_ = std::vector<std::vector<std::vector<double>>>();
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

    // Check that the determinants have the same number of occupied orbitals
    if (di._occupiedOrbitals[ALPHA].size() != dj._occupiedOrbitals[ALPHA].size() ||
        di._occupiedOrbitals[BETA].size() != dj._occupiedOrbitals[BETA].size())
    {
        std::stringstream errorMessage;
        errorMessage << "Error in SlaterDeterminant::getDifferences(): Slater determinants have different numbers of occupied orbitals." << std::endl;
        print_error(errorMessage.str());

        exit(1);
    }

    // Check for differences in occupied orbitals for each spin type
    for (size_t i = 0; i < di._occupiedOrbitals[ALPHA].size(); ++i)
    {
        // Alpha spin
        if (di._occupiedOrbitals[ALPHA][i].first != dj._occupiedOrbitals[ALPHA][i].first)
        {
            differences[ALPHA].emplace_back(di._occupiedOrbitals[ALPHA][i].first, dj._occupiedOrbitals[ALPHA][i].first);
        }

        // Beta spin
        if (di._occupiedOrbitals[BETA][i].first != dj._occupiedOrbitals[BETA][i].first)
        {
            differences[BETA].emplace_back(di._occupiedOrbitals[BETA][i].first, dj._occupiedOrbitals[BETA][i].first);
        }
    }

    return differences;
}

double SlaterDeterminant::ionicPotential(const SlaterDeterminant& di, const SlaterDeterminant& dj, const std::array<double, 3>& position, double charge)
{
    if (_s_ionicMatrix_.empty())
    {
        if (_s_isOrbitalsSet_)
        {
            _s_ionicMatrix_ = _s_orbitals_.getIonicPotentialMatrix(position, charge);
        }
        else
        {
            std::stringstream errorMessage;
            errorMessage << "Error in SlaterDeterminant::ionicPotential(): shared Orbitals instance has not been set yet." << std::endl;
            errorMessage << "Please build at least the ground state Slater determinant." << std::endl;
            print_error(errorMessage.str());

            exit(1);
        }
    }

    double sum = 0.0;

    // Get spin type int values
    const int ALPHA = static_cast<int>(SpinType::ALPHA);
    const int BETA = static_cast<int>(SpinType::BETA);

    // Apply Slater-Condon rules
    if (di == dj)
    {
        for (size_t i = 0; i < di._occupiedOrbitals[ALPHA].size(); ++i)
        {
            int orbitalIndex = di._occupiedOrbitals[ALPHA][i].first - 1;
            sum += (_s_ionicMatrix_[ALPHA][orbitalIndex][orbitalIndex] + _s_ionicMatrix_[BETA][orbitalIndex][orbitalIndex]);
        }
    }
    else
    {
        std::vector<std::vector<std::pair<int, int>>> differences = getDifferences(di, dj);

        // Check that there is only one difference
        if (differences[ALPHA].size() + differences[BETA].size() == 1)
        {
            int initialOrbitalIndex, finalOrbitalIndex;
            
            if (!differences[ALPHA].empty())
            {
                initialOrbitalIndex = differences[ALPHA][0].first - 1;
                finalOrbitalIndex = differences[ALPHA][0].second - 1;
            }
            else
            {
                initialOrbitalIndex = differences[BETA][0].first - 1;
                finalOrbitalIndex = differences[BETA][0].second - 1;
            }

            sum += (_s_ionicMatrix_[ALPHA][initialOrbitalIndex][finalOrbitalIndex] + _s_ionicMatrix_[BETA][initialOrbitalIndex][finalOrbitalIndex]);
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