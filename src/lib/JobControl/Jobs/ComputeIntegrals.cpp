#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <Common/PeriodicTable.h>
#include <Cube/Grid.h>
#include <Cube/GridCP.h>
#include <JobControl/Job.h>
#include <JobControl/Jobs/ComputeIntegrals.hpp>
#include <Utils/Enums.hpp>
#include <Utils/Utils.h>


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//----------------------------------------------------------------------------------------------------//

ComputeIntegrals::ComputeIntegrals(const std::string& inputFileName):
    Job(inputFileName)
{ }


//----------------------------------------------------------------------------------------------------//
// STATIC METHODS
//----------------------------------------------------------------------------------------------------//

void ComputeIntegrals::buildBasins(GridCP& gridCP, const std::string& gridFileName, PartitionMethod partitionMethod)
{
    std::ifstream gridFile(gridFileName);
    if (!gridFile)
    {
        print_error("Error in ComputeIntegrals::buildBasins(): could not read file " + gridFileName + ".");
        std::exit(1);
    }

    Grid g(gridFile, PeriodicTable());
    gridFile.close();

    gridCP.buildBasins(g, partitionMethod);
}

void ComputeIntegrals::buildBasinsBySign(GridCP& gridCP, const std::string& gridFileName,  double cutoff, bool two) 
{
    std::ifstream gridFile(gridFileName);
    if (!gridFile)
    {
        print_error("Error in ComputeIntegrals::buildBasinsBySign(): could not read file " + gridFileName + ".");
        std::exit(1);
    }

    Grid g(gridFile, PeriodicTable());
    gridFile.close();

    if (two)
    {
        gridCP.build2BasinSign(g);
    }    
    else
    {
        gridCP.buildBasinsBySign(g, cutoff);
    }    
}

void ComputeIntegrals::computeLocalIntegrals(GridCP& gridCP, const std::vector<std::string>& gridFilesNames, PartitionMethod partitionMethod, double cutoff)
{
    // Check partition method validity
    if (partitionMethod == PartitionMethod::BECKE || partitionMethod == PartitionMethod::FD
                                                  || partitionMethod == PartitionMethod::FMO)
    {
        std::stringstream errorMessage;
        errorMessage << "Error in ComputeIntegrals::computeLocalIntegrals(): partitionMethod \"" << to_string(partitionMethod) << "\" invalid." << std::endl;
        errorMessage << "Please check the documentation of the \"ComputeIntegrals\" job.";

        print_error(errorMessage.str());

        std::exit(1);
    }

    // Print partition method
    std::cout << "Volume partition method: " << to_string(partitionMethod) << std::endl << std::endl;

    // Build basins
    std::cout << "Reading file " << gridFilesNames[0] << " to build basins." << std::endl << std::endl;

    if (partitionMethod == PartitionMethod::BBS)
    {
        buildBasinsBySign(gridCP, gridFilesNames[0], cutoff, false);
    }
    else if (partitionMethod == PartitionMethod::B2S)
    {
        buildBasinsBySign(gridCP, gridFilesNames[0], cutoff, true);
    }
    else
    {
        buildBasins(gridCP, gridFilesNames[0], partitionMethod);
    }

    
    // Compute local integrals
    for(size_t i = 0; i < gridFilesNames.size(); ++i)
    {
        std::ifstream gridFile(gridFilesNames[i]);
        if (!gridFile)
        {
            print_error("Error in ComputeIntegrals::computeLocalIntegrals(): could not read file " + gridFilesNames[i] + ".");
            std::exit(1);
        }

        Grid g(gridFile, PeriodicTable());
        gridFile.close();

        gridCP.computeIntegrals(g);
        gridCP.printCriticalPoints();
    }
}


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

void ComputeIntegrals::run()
{
    // Read grid files names
    std::vector<std::string> gridFilesNames;
    readGridFilesNames(gridFilesNames);


    // Read partition method
    PartitionMethod partitionMethod;
    readPartitionMethod(partitionMethod);


    // Read cutoff
    double cutoff = 0.0;
    readCutoff(cutoff);

    
    // Compute local integrals
    GridCP gridCP;
    computeLocalIntegrals(gridCP, gridFilesNames, partitionMethod, cutoff);
}
