#include <string>

#include <Cube/Domain.h>
#include <JobControl/Job.h>
#include <JobControl/Jobs/MakeOrbitalsCube.hpp>
#include <Orbitals/Orbitals.h>
#include <Utils/Enums.hpp>


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//----------------------------------------------------------------------------------------------------//

MakeOrbitalsCube::MakeOrbitalsCube(const std::string& inputFileName):
    Job(inputFileName)
{ }


//----------------------------------------------------------------------------------------------------//
// STATIC METHODS
//----------------------------------------------------------------------------------------------------//

void MakeOrbitalsCube::createCube(const std::string& analyticFileName, const std::string& outputFileName, OrbitalType orbitalType, SpinType spinType, GridSize gridSize, CustomSizeData customSizeData, bool showProgress, std::vector<int> orbitalsNumbers, std::vector<SpinType> spinList)
{
    // Loading orbitals
    Orbitals orbitals;
    computeOrbitalsOrBecke<Orbitals>(orbitals, analyticFileName);

    // Setting orbitals
    std::vector<SpinType> orbitalsSpins;
    setOrbitals(orbitals, orbitals.get_numberOfMo(), orbitalsNumbers, orbitalsSpins, orbitalType, spinType, spinList);

    // Building domain
    Domain domain = buildDomainForCube(orbitals, gridSize, customSizeData, orbitalsNumbers.size());

    // Saving orbitals cube file
    Job::createCube(orbitals, domain, outputFileName, CubeType::ORBITALS, showProgress, ELFMethod::UNKNOWN, orbitalsNumbers, orbitalsSpins);
}


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

void MakeOrbitalsCube::run()
{
    // Read analytic file names
    std::vector<std::string> analyticFilesNames;
    readAnalyticFilesNames(analyticFilesNames);

    // Check number of analytic files names
    if (analyticFilesNames.size() != 1)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect number of analytic files names (one file expected)." << std::endl;
        errorMessage << "Please check the documentation and the number of files specified in the \"AnalyticFiles\" parameter in " << _inputFileName << '.';

        print_error(errorMessage.str());

        std::exit(1);
    }

    // Read grid files name
    std::vector<std::string> gridFileName;
    readGridFilesNames(gridFileName);

    // Check number of grid files names
    if (gridFileName.size() != 1)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect number of grid files names (one file expected)." << std::endl;
        errorMessage << "Please check the documentation and the number of files specified in the \"GridFiles\" parameter in " << _inputFileName << '.';

        print_error(errorMessage.str());

        std::exit(1);
    }

    // Read grid parameters
    GridSize gridSize;
    CustomSizeData customSizeData;
    readSize(gridSize, customSizeData);
    
    // Read spin type
    SpinType spinType;
    readSpinType(spinType);

    // Read orbitals type
    OrbitalType orbitalType;
    readOrbitalType(orbitalType);

    std::vector<int> orbitalsNumbers;
    std::vector<SpinType> spinList;
    if (orbitalType == OrbitalType::CUSTOM)
    {
        // Read custom orbitals list
        readOrbitalsNumbers(orbitalsNumbers);
        readSpinList(spinList);
        

        // Check if the sizes of both lists match
        if (spinList.size() != orbitalsNumbers.size())
        {
            std::stringstream errorMessage;
            errorMessage << "Error: sizes of orbitals numbers list and spins list do not match." << std::endl;
            errorMessage << "Please check the documentation and the number of items specified in the \"OrbitalsList\" and \"SpinList\" parameters in " << _inputFileName << '.';

            print_error(errorMessage.str());

            std::exit(1);
        }
    }

    // Read progress bar display option
    bool showProgress;
    readShowProgress(showProgress);

    // Creating orbitals cube
    createCube(analyticFilesNames[0], gridFileName[0], orbitalType, spinType, gridSize, customSizeData, showProgress);
}