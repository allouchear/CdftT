#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <numeric>

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

void ComputeElectronDensity::computeStateDensities(const std::vector<ExcitedState>& states,std::vector<int>& excitedStatesNumbers, const Orbitals& orbitals, Grid& grid, const RDMMethod rdmMethod, const std::vector<int>& excludedOrbitalsNumbers, const std::string& outputPrefix, bool saveRDM, int outputPrecision, int verbose, std::ostream& logOutputStream, bool showProgress)
{
    size_t nbStates = states.size();

    std::stringstream logStream;
    std::ofstream outputFile;

    int iter = 0;
    for (size_t i : excitedStatesNumbers)
    {
        iter += 1;
        std::vector<std::vector<std::vector<double>>> reducedDensityMatrix;
        ExcitedState::reducedDensityMatrix(reducedDensityMatrix, rdmMethod, states[i], states[i], orbitals, excludedOrbitalsNumbers);

        if (saveRDM || verbose >= 1)
        {
            if (saveRDM)
            {
                std::string outputFileName = outputPrefix + "_state" + std::to_string(states[i].get_number()) + "_RDM.cdftt";
                outputFile.open(outputFileName);
                if (!outputFile)
                {
                    std::stringstream errorMessage;
                    errorMessage << "Error in ComputeElectronDensity::run(): could not open output file " << outputFileName << " for writing." << std::endl;

                    print_error(errorMessage.str(), logStream);
                    
                    std::exit(1);
                }

                outputFile << std::scientific;
                outputFile << std::setprecision(outputPrecision);
            }

            if (verbose >= 1)
            {
                logStream << std::scientific;
                logStream << std::setprecision(10);
            }

            for (int spin = 0; spin < 2; ++spin)
            {
                if (verbose >= 1)
                {
                    logStream << "Reduced Density Matrix for state #" << states[i].get_number() << " (spin " << ((spin == 0) ? "alpha" : "beta") << "):" << std::endl;
                }
                if (saveRDM)
                {
                    outputFile << "Spin " << ((spin == 0) ? "alpha" : "beta") << ':' << std::endl;
                }

                for(size_t i = 0; i < reducedDensityMatrix[spin].size(); ++i)
                {
                    for(size_t j = 0; j < reducedDensityMatrix[spin][i].size(); ++j)
                    {
                        if (verbose >= 1)
                        {
                            logStream << std::right << std::setw(17) << reducedDensityMatrix[spin][i][j] << '\t';
                        }
                        if (saveRDM)
                        {
                            outputFile << std::right << std::setw(17) << reducedDensityMatrix[spin][i][j] << '\t';
                        }
                    }
                    if (verbose >= 1)
                    {
                        logStream << std::endl;
                    }
                    if (saveRDM)
                    {
                        outputFile << std::endl;
                    }
                }
                if (verbose >= 1)
                {
                    logStream << std::endl;
                }
                if (saveRDM)
                {
                    outputFile << std::endl;
                }
            }
            if (verbose >= 1)
            {
                logStream << std::defaultfloat << std::endl;
                log(logStream, logOutputStream);
            }
            if (saveRDM)
            {
                outputFile.close();
            }
        }

        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Computing electronic density for state #" << states[i].get_number() << " (" << iter << " out of " << excitedStatesNumbers.size() << "), please wait..." << std::endl;
        grid.reset();
        orbitals.makeDensityGrid(grid, reducedDensityMatrix, showProgress);
        auto end = std::chrono::high_resolution_clock::now();
        if (showProgress)
        {
            std::cout << std::endl;
        }
        std::cout<<"compute time of grid : "<<std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()/1000.0<<"s"<<std::endl;
        start = std::chrono::high_resolution_clock::now();
        // Save grid
        std::cout << "Writing density cube file for state #" << states[i].get_number() << " (" << iter << " out of " << excitedStatesNumbers.size() << "), please wait..." << std::endl;
        std::ofstream out(outputPrefix + "_state" + std::to_string(states[i].get_number()) + ".cube");
        grid.save(out, showProgress, outputPrecision);
        out.close();
        end = std::chrono::high_resolution_clock::now();
        std::cout<<"compute time of writting grid : "<<std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()/1000.0<<"s"<<std::endl;
        if (showProgress)
        {
            std::cout << std::endl;
        }
    }
}

void ComputeElectronDensity::computeTransitionDensities(const std::vector<ExcitedState>& states, const std::vector<std::array<int, 2>>& transitionDensities, const Orbitals& orbitals, Grid& grid, const RDMMethod rdmMethod, const std::vector<int>& excludedOrbitalsNumbers, const std::string& outputPrefix, bool saveRDM, int outputPrecision, int verbose, std::ostream& logOutputStream, bool showProgress)
{
    size_t nbTransitionDensities = transitionDensities.size();

    std::stringstream logStream;
    std::ofstream outputFile;

    int currentTransitionDensity = 1;
    for (const std::array<int, 2>& transitionDensity : transitionDensities)
    {
        int i = transitionDensity[0];
        int j = transitionDensity[1];

        std::vector<std::vector<std::vector<double>>> reducedDensityMatrix;
        ExcitedState::reducedDensityMatrix(reducedDensityMatrix, rdmMethod, states[i], states[j], orbitals, excludedOrbitalsNumbers);

        if (saveRDM || verbose >= 1)
        {
            if (saveRDM)
            {
                std::string outputFileName = outputPrefix + "_transition_state" + std::to_string(states[i].get_number()) + "_state" + std::to_string(states[j].get_number()) + "_RDM.cdftt";
                outputFile.open(outputFileName);
                if (!outputFile)
                {
                    std::stringstream errorMessage;
                    errorMessage << "Error in ComputeElectronDensity::run(): could not open output file " << outputFileName << " for writing." << std::endl;

                    print_error(errorMessage.str(), logStream);
                    
                    std::exit(1);
                }

                outputFile << std::scientific;
                outputFile << std::setprecision(outputPrecision);
            }

            if (verbose >= 1)
            {
                logStream << std::scientific;
                logStream << std::setprecision(10);
            }

            for (int spin = 0; spin < 2; ++spin)
            {
                if (verbose >= 1)
                {            
                    logStream << "Reduced Density Matrix for transition between state #" << states[i].get_number() << " and state #" << states[j].get_number() << " (spin " << ((spin == 0) ? "alpha" : "beta") << "):" << std::endl;
                }

                if (saveRDM)
                {
                    outputFile << "Spin " << ((spin == 0) ? "alpha" : "beta") << ':' << std::endl;
                }

                for(size_t i = 0; i < reducedDensityMatrix[spin].size(); ++i)
                {
                    for(size_t j = 0; j < reducedDensityMatrix[spin][i].size(); ++j)
                    {
                        if (verbose >= 1)
                        {
                            logStream << std::right << std::setw(17) << reducedDensityMatrix[spin][i][j] << '\t';
                        }

                        if (saveRDM)
                        {
                            outputFile << std::right << std::setw(17) << reducedDensityMatrix[spin][i][j] << '\t';
                        }
                    }
                    if (verbose >= 1)
                    {
                        logStream << std::endl;
                    }

                    if (saveRDM)
                    {
                        outputFile << std::endl;
                    }
                }
                if (verbose >= 1)
                {
                    logStream << std::endl;
                }

                if (saveRDM)
                {
                    outputFile << std::endl;
                }
            }
            if (verbose >= 1)
            {
                logStream << std::defaultfloat << std::endl;
                log(logStream, logOutputStream);
            }
            
            if (saveRDM)
            {
                outputFile.close();
            }
        }
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Computing electronic density for transition between state #" << states[i].get_number() << " and state #" << states[j].get_number() << " (" << currentTransitionDensity << " out of " << nbTransitionDensities << "), please wait..." << std::endl;
        grid.reset();
        orbitals.makeDensityGrid(grid, reducedDensityMatrix, showProgress);
        auto end = std::chrono::high_resolution_clock::now();
        if (showProgress)
        {
            std::cout << std::endl;
        }
        std::cout<<"compute time of writting grid : "<<std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()/1000.0<<"s"<<std::endl;
        start = std::chrono::high_resolution_clock::now();
        // Save grid
        std::cout << "Writing density cube file for transition between state #" << states[i].get_number() << " and state #" << states[j].get_number() << " (" << currentTransitionDensity << " out of " << nbTransitionDensities << "), please wait..." << std::endl;
        std::ofstream out(outputPrefix + "_transition_state" + std::to_string(states[i].get_number()) + "_state" + std::to_string(states[j].get_number()) + ".cube");
        grid.save(out, showProgress, outputPrecision);
        out.close();
        end = std::chrono::high_resolution_clock::now();
        std::cout<<"compute time of writting grid : "<<std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()/1000.0<<"s"<<std::endl;

        if (showProgress)
        {
            std::cout << std::endl;
        }

        ++currentTransitionDensity;
    }
}


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

void ComputeElectronDensity::run()
{
    auto start = std::chrono::high_resolution_clock::now();

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

        print_error(errorMessage.str(), logStream);

        std::exit(1);
    }


    // Read grid parameters
    GridSize gridSize;
    CustomSizeData customSizeData;
    readSize(gridSize, customSizeData);


    // Read transitions file name
    std::string transitionsFileName;
    readTransitionsFileName(transitionsFileName);


    // Read max number of transitions to read in transitions file
    int maxNbExcitedStates;
    readMaxNumberOfExcitedStates(maxNbExcitedStates);


    // Read excited states numbers (1-based) to keep
    std::vector<int> excitedStatesNumbers;
    readExcitedStatesNumbers(excitedStatesNumbers);

    // Read transition densities to consider in the transitions reduced density matrix computation
    std::vector<std::array<int, 2>> transitionDensities;
    readTransitionDensities(transitionDensities);

    //get highest state number among excitedStatesNumbers and transitionDensities
    int highestState = 0;
    for (auto pair : transitionDensities)
    {
        highestState = std::max(highestState,std::max(pair[0],pair[1]));
    }
    int maxExcitedStates;
    if  (excitedStatesNumbers.size()>0)
        maxExcitedStates = *std::max_element(excitedStatesNumbers.begin(),excitedStatesNumbers.end());
    else
        maxExcitedStates = 0;
    highestState = std::max(highestState,maxExcitedStates);

    std::vector<int> statesToCompute;
    for (int i=0;i<=highestState;++i)
    {
        statesToCompute.push_back(i);
    }
    auto end = std::chrono::high_resolution_clock::now();
    double t_readFiles = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()/1000.0;

    // Load orbitals
    start = std::chrono::high_resolution_clock::now();
    Orbitals orbitals;
    computeOrbitalsOrBecke<Orbitals>(orbitals, analyticFilesNames[0]);
    end = std::chrono::high_resolution_clock::now();
    double t_loadOrbitals = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()/1000.0;

    // Get Ground Slater Determinant
    SlaterDeterminant groundStateSlaterDeterminant(orbitals);
    ExcitedState groundState(orbitals.get_energy(), groundStateSlaterDeterminant);

    // Build states vector
    std::vector<ExcitedState> states;
    states.push_back(groundState);

    // Read transitions file
    start = std::chrono::high_resolution_clock::now();
    if (!transitionsFileName.empty())
    {
        std::cout << "Reading transitions from file: " << transitionsFileName << ". Please wait..." << std::endl;
        ExcitedState::readTransitions(transitionsFileName, states, groundState.get_energy(), maxNbExcitedStates, statesToCompute);
    }
    else
    {
        std::cout << "Reading transitions from analytic file: " << analyticFilesNames[0] << ". Please wait..." << std::endl;
        ExcitedState::readTransitions(analyticFilesNames[0], states, groundState.get_energy(), maxNbExcitedStates, statesToCompute);
    }
    end = std::chrono::high_resolution_clock::now();
    double t_getTransitions = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()/1000.0;;

    size_t nbStates = states.size();
    logStream << "Total number of states: " << nbStates;
    if (verbose>=1) 
    {
        std::cout << std::endl << std::endl;
    }
    log(logStream, outputStream);

    double SDLimit;
    readSDLimit(SDLimit);

    // Compute Slater Determinants from electronic transitions for each state, and set the argsort array
    start = std::chrono::high_resolution_clock::now();
    for (ExcitedState& state : states)
    {
        state.computeSlaterDeterminants(groundStateSlaterDeterminant);
        state.set_argsortCoefs(SDLimit);

        if (verbose >= 1)
        {
            logStream << state << std::endl;
            log(logStream, outputStream);
        }
    }
    std::cout << std::endl;
    end = std::chrono::high_resolution_clock::now();
    double t_computeSD = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()/1000.0;;

    // Read orbitals numbers to exclude from the density computation
    std::vector<int> excludedOrbitals;
    readExcludedOrbitals(excludedOrbitals);
    std::cout<<"excluded orbitals : ";
    for (int i=0;i<excludedOrbitals.size();++i)
    {
        std::cout<<excludedOrbitals[i];
        if (i!=excludedOrbitals.size()-1) {std::cout<<",";}
    }
    std::cout<<std::endl<<std::endl<<std::endl;

    // Read method to use for Reduced Density Matrix computation
    RDMMethod rdmMethod;
    readRDMMethod(rdmMethod);


    // Read the option to save the reduced density matrix
    bool saveRDM;
    readSaveReducedDensityMatrix(saveRDM);


    // Read precision for numerical values in the output files
    int outputPrecision;
    readPrecision(outputPrecision);

    start = std::chrono::high_resolution_clock::now();
    // Build domain
    //std::cout << "Building domain and grid, please wait..." << std::endl;
    Domain domain = buildDomainForCube(orbitals, gridSize, customSizeData, 1);
    Grid grid;
    grid.set_structure(orbitals.get_struct());
    grid.set_domain(domain);
    end = std::chrono::high_resolution_clock::now();
    double t_buildDomaine = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()/1000.0;


    start = std::chrono::high_resolution_clock::now();
    logStream << "Total number of electronic densities to compute: " << excitedStatesNumbers.size() << std::endl << std::endl;
    log(logStream, outputStream);
    // Compute state electronic densities and save them in .cube files
    computeStateDensities(states, excitedStatesNumbers, orbitals, grid, rdmMethod, excludedOrbitals, outputPrefix, saveRDM, outputPrecision, verbose, outputStream, showProgress);
    end = std::chrono::high_resolution_clock::now();
    double t_computeDensity = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()/1000.0;


    start = std::chrono::high_resolution_clock::now();
    size_t nbTransitionDensities = transitionDensities.size();
    logStream << "Total number of transition densities to compute: " << nbTransitionDensities << std::endl << std::endl;
    log(logStream, outputStream);
    // Compute requested transition densities
    computeTransitionDensities(states, transitionDensities, orbitals, grid, rdmMethod, excludedOrbitals, outputPrefix, saveRDM, outputPrecision, verbose, outputStream, showProgress);
    end = std::chrono::high_resolution_clock::now();
    double t_computeTransitions = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()/1000.0;

    
    std::cout<<"time to read files : "<<t_readFiles<<"s"<<std::endl;
    std::cout<<"time to load orbitals : "<<t_loadOrbitals<<"s"<<std::endl;
    std::cout<<"time to get transitions : "<<t_getTransitions<<"s"<<std::endl;
    std::cout<<"time to compute SD + sort coefs: "<<t_computeSD<<"s"<<std::endl;
    std::cout<<"time to build domaine : "<<t_buildDomaine<<"s"<<std::endl;
    std::cout<<"time to compute densities : "<<t_computeDensity<<"s"<<std::endl;
    std::cout<<"time to compute transitions : "<<t_computeTransitions<<"s"<<std::endl;
}


