#include <string>

#include <Cube/Domain.h>
#include <JobControl/Job.h>
#include <JobControl/Jobs/MakeDensityCube.hpp>
#include <Orbitals/Orbitals.h>
#include <Utils/Enums.hpp>


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//----------------------------------------------------------------------------------------------------//

MakeDensityCube::MakeDensityCube(const std::string& inputFileName):
    Job(inputFileName)
{ }


//----------------------------------------------------------------------------------------------------//
// STATIC METHODS
//----------------------------------------------------------------------------------------------------//

void MakeDensityCube::createCube(const std::string& analyticFileName, const std::string& outputFileName, GridSize gridSize, CustomSizeData customSizeData, bool showProgress)
{
    // Loading orbitals
    Orbitals orbitals;
    computeOrbitalsOrBecke<Orbitals>(orbitals, analyticFileName);

    // Building domain
    Domain domain = buildDomainForCube(orbitals, gridSize, customSizeData, 1);

    // Saving density cube file
    Job::createCube(orbitals, domain, outputFileName, CubeType::DENSITY, showProgress);
}


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

void MakeDensityCube::run()
{
    // Read analytic file name
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

    // Read grid parameters
    GridSize gridSize;
    CustomSizeData customSizeData;
    readSize(gridSize, customSizeData);

    // Read grid file name
    std::vector<std::string> gridFilesName;
    readGridFilesNames(gridFilesName);

    // Check number of grid files names
    if (gridFilesName.size() != 1)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect number of grid files names (one file expected)." << std::endl;
        errorMessage << "Please check the documentation and the number of files specified in the \"GridFiles\" parameter in " << _inputFileName << '.';

        print_error(errorMessage.str());

        std::exit(1);
    }

    // Read progress bar display option
    bool showProgress;
    readShowProgress(showProgress);

    // Create density cube
    createCube(analyticFilesNames[0], gridFilesName[0], gridSize, customSizeData, showProgress);
}