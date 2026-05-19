#include <iostream>
#include <string>
#include <vector>

#include <Cube/Domain.h>
#include <Cube/Grid.h>
#include <JobControl/Job.h>
#include <JobControl/Jobs/LambdaDiagnostic.hpp>
#include <Orbitals/Orbitals.h>
#include <Orbitals/ExcitedState.hpp>
#include <Utils/Enums.hpp>


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//----------------------------------------------------------------------------------------------------//

LambdaDiagnostic::LambdaDiagnostic(const std::string& inputFileName):
    Job(inputFileName)
{ }


//----------------------------------------------------------------------------------------------------//
// STATIC METHODS
//----------------------------------------------------------------------------------------------------//

void LambdaDiagnostic::print(const std::string& analyticFileName, const std::string& transitionsFileName, GridSize gridSize, CustomSizeData customSizeData)
{
    // Loading orbitals
    Orbitals orbitals;
    computeOrbitalsOrBecke<Orbitals>(orbitals, analyticFileName);


    // Setting orbitals
    std::vector<SpinType> orbitalsSpins;
    std::vector<int> orbitalsNumbers;
    setOrbitals(orbitals, orbitals.get_numberOfMo(), orbitalsNumbers, orbitalsSpins, OrbitalType::ALL, SpinType::ALPHA_BETA);


    // Building domain and grid
    std::cout << "Building domain and grid, please wait..." << std::endl;

    Domain domain = buildDomainForCube(orbitals, gridSize, customSizeData, orbitalsNumbers.size());
    Grid orbitalsGrid = orbitals.makeOrbGrid(domain, orbitalsNumbers, orbitalsSpins);


    // Reading transitions file
    std::vector<ExcitedState> excitedStates;
    ExcitedState::readTransitions(transitionsFileName, excitedStates, orbitals.get_energy());
    std::cout << "Number of excited states read: " << excitedStates.size() << std::endl << std::endl;


    // Printing lambda diagnostic for each excited state
    for (const ExcitedState& excitedState : excitedStates)
    {
        std::cout << excitedState << std::endl;
        excitedState.printLambdaDiagnostic(orbitalsGrid);
        std::cout << "-----------------" << std::endl << std::endl;
    }
}


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

void LambdaDiagnostic::run()
{
    //Read analytic file name
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


    // Read grid size
    GridSize gridSize;
    CustomSizeData customSizeData;
    readSize(gridSize, customSizeData);


    // Read transitions file
    std::string transitionsFileName;
    readTransitionsFileName(transitionsFileName);


    // Compute and print lambda diagnostic
    print(analyticFilesNames[0], transitionsFileName, gridSize, customSizeData);
}