#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <JobControl/Job.h>
#include <JobControl/Jobs/ComputeElectronDensity.hpp>
#include <Orbitals/ExcitedState.hpp>
#include <Orbitals/Orbitals.h>
#include <Orbitals/SlaterDeterminant.hpp>
#include <Utils/Enums.hpp>
#include <Utils/Utils.h>


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//----------------------------------------------------------------------------------------------------//

ComputeElectronDensity::ComputeElectronDensity(const std::string& inputFileName):
    Job(inputFileName)
{ }


//----------------------------------------------------------------------------------------------------//
// STATIC METHODS
//----------------------------------------------------------------------------------------------------//

void ComputeElectronDensity::compute()
{

}


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

void ComputeElectronDensity::run()
{
    // Read output file prefix
    std::string outputPrefix;
    readOutputPrefix(outputPrefix);


    // Read verbose level and open log file if needed
    int verbose;
    readVerbose(verbose);

    std::stringstream logStream;

    std::ofstream logFile;
    if (verbose != 0)
    {
        logFile.open(outputPrefix + "_log.cdftt");
        if (!logFile)
        {
            std::cout << "Warning: could not open log file " << outputPrefix << "_log.cdftt for writing." << std::endl;
            std::cout << "The program will still display logging information on standard output." << std::endl << std::endl;
        }
    }
    std::ostream& outputStream = ((verbose != 0 && logFile) ? logFile : std::cout);


    // Read progress bar display option
    bool showProgress;
    readShowProgress(showProgress);


    // Read analytic files names
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


    // Read transitions file
    std::string transitionsFileName;
    readTransitionsFileName(transitionsFileName);


    // Read max number of transitions to read in transitions file
    int maxNbExcitedStates;
    readMaxNumberOfExcitedStates(maxNbExcitedStates);


    // Read states numbers (0-based) to keep
    std::vector<int> statesNumbers;
    readStatesNumbers(statesNumbers);


    // Load orbitals
    Orbitals orbitals;
    computeOrbitalsOrBecke<Orbitals>(orbitals, analyticFilesNames[0]);


    // Get Ground Slater Determinant
    SlaterDeterminant groundStateSlaterDeterminant(orbitals);
    ExcitedState groundState(orbitals.get_energy(), groundStateSlaterDeterminant);


    // Build states vector
    std::vector<ExcitedState> states;
    states.push_back(groundState);


    // Read transitions file
    if (!transitionsFileName.empty())
    {
        std::cout << "Reading transitions from file: " << transitionsFileName << ". Please wait..." << std::endl;
        ExcitedState::readTransitions(transitionsFileName, states, groundState.get_energy(), maxNbExcitedStates, statesNumbers);
    }
    else
    {
        std::cout << "Reading transitions from analytic file: " << analyticFilesNames[0] << ". Please wait..." << std::endl;
        ExcitedState::readTransitions(analyticFilesNames[0], states, groundState.get_energy(), maxNbExcitedStates, statesNumbers);
    }

    size_t nbStates = states.size();
    logStream << "Total number of states: " << nbStates << std::endl << std::endl;
    log(logStream, outputStream);


    // Compute Slater Determinants from electronic transitions for each state, 
    int state_index_ = 0;
    for (ExcitedState& state : states)
    {
        state.computeSlaterDeterminants(groundStateSlaterDeterminant);
        ++state_index_;

        if(showProgress)
        {
            std::cout << "computed SD for state " << state_index_ << "/" << nbStates << " of energy " << state.get_energy()  <<" Hartree (dE=" << state.get_energy()-states[0].get_energy() << "); state got " << state.get_slaterDeterminants().size() << " Slater determinants" << std::endl;
        }
    }
    std::cout << std::endl;


    // Build domain
    std::cout << "Building domain and grid, please wait..." << std::endl;
    Domain domain = buildDomainForCube(orbitals, gridSize, customSizeData, 1);


    // Compute grid values for each excited states
    for (size_t i = 0; i < nbStates; ++i)
    {
        Grid grid;
        grid.set_structure(orbitals.get_struct());
        grid.set_domain(domain);
        grid.reset();

        std::cout << "Computing grid for state #" << states[i].get_number() << " (" << i + 1 << " out of " << nbStates << "), please wait..." << std::endl;
        states[i].makeDensityGrid(orbitals, grid, showProgress);

        if (showProgress)
        {
            std::cout << std::endl;
        }

        // Save grid
        std::cout << "Writing cube file for state #" << states[i].get_number() << " (" << i + 1 << " out of " << nbStates << "), please wait..." << std::endl;
        std::ofstream out(outputPrefix + "_state" + std::to_string(states[i].get_number()) + ".cube");
        grid.save(out, showProgress);
        out.close();

        if (showProgress)
        {
            std::cout << std::endl;
        }
    }
}


