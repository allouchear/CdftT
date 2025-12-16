#include <algorithm>
#include <cstdlib>
#include <vector>

#include <Orbitals/Orbitals.h>
#include <Orbitals/SlaterDeterminant.hpp>


SlaterDeterminant::SlaterDeterminant():
    _alphaOccupation(),
    _betaOccupation()
{ }


SlaterDeterminant::SlaterDeterminant(const Orbitals& orbitals):
    _alphaOccupation(),
    _betaOccupation()
{
    // Get occupation numbers
    const std::vector<std::vector<double>>& occupationNumbers = orbitals.OccupationNumber();

    // Populate the occupied orbitals and their occupation numbers based on the orbitals object
    for (int i = 0; i < orbitals.NumberOfMo(); ++i)
    {
        // Alpha spin
        if (occupationNumbers[0][i] == 1)
        {
            _alphaOccupation.emplace_back(i + 1, 1.0);
        }
        else if (occupationNumbers[0][i] == 2) // Case where _alpha_and_beta = true
        {
            // Alpha spin
            _alphaOccupation.emplace_back(i + 1, 1.0);

            // Beta spin
            _betaOccupation.emplace_back(i + 1, 1.0);
        }

        // Beta spin
        if (occupationNumbers[1][i] == 1)
        {
            _betaOccupation.emplace_back(i + 1, 1.0);
        }
    }
}


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

void SlaterDeterminant::updateFromTransition(int initialOrbitalNumber, SpinType initialSpin, int finalOrbitalNumber, SpinType finalSpin)
{
    // Update occupied orbitals based on the transition
    if (initialSpin == SpinType::ALPHA && finalSpin == SpinType::ALPHA)
    {
        for (size_t i = 0; i < _alphaOccupation.size(); ++i)
        {
            if (_alphaOccupation[i].first == initialOrbitalNumber)
            {
                _alphaOccupation[i].first = finalOrbitalNumber;
                break;
            }
        }

        // Sort occupied orbitals to maintain order (facilitates comparison)
        std::sort(_alphaOccupation.begin(), _alphaOccupation.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    }
    else if (initialSpin == SpinType::BETA && finalSpin == SpinType::BETA)
    {
        for (size_t i = 0; i < _betaOccupation.size(); ++i)
        {
            if (_betaOccupation[i].first == initialOrbitalNumber)
            {
                _betaOccupation[i].first = finalOrbitalNumber;
                break;
            }
        }

        // Sort occupied orbitals to maintain order (facilitates comparison)
        std::sort(_betaOccupation.begin(), _betaOccupation.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
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
// FRIEND FUNCTIONS
//----------------------------------------------------------------------------------------------------//

bool operator==(const SlaterDeterminant& lhs, const SlaterDeterminant& rhs)
{
    return (lhs._alphaOccupation == rhs._alphaOccupation) && (lhs._betaOccupation == rhs._betaOccupation);
}

bool operator!=(const SlaterDeterminant& lhs, const SlaterDeterminant& rhs)
{
    return !(lhs == rhs);
}

std::ostream& operator<<(std::ostream& stream, const SlaterDeterminant& slaterDeterminant)
{
    for (size_t i = 0; i < slaterDeterminant._alphaOccupation.size(); ++i)
    {
        stream << slaterDeterminant._alphaOccupation[i].first << "A(" << slaterDeterminant._alphaOccupation[i].second << ") ";
    }

    stream << std::endl;

    for (size_t i = 0; i < slaterDeterminant._betaOccupation.size(); ++i)
    {
        stream << slaterDeterminant._betaOccupation[i].first << "B(" << slaterDeterminant._betaOccupation[i].second << ") ";
    }

    return stream;
}