#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <JobControl/Job.h>
#include <JobControl/Jobs/ComputePartialCharges.hpp>
#include <Utils/Enums.hpp>
#include <Utils/Utils.h>

//----------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//----------------------------------------------------------------------------------------------------//

ComputePartialCharges::ComputePartialCharges(const std::string& inputFileName):
    Job(inputFileName)
{ }


//----------------------------------------------------------------------------------------------------//
// STATIC METHODS
//----------------------------------------------------------------------------------------------------//

std::vector<double> ComputePartialCharges::compute(const std::string& gridFileName, PartitionMethod partitionMethod)
{
    // Check partition method validity
    if (partitionMethod == PartitionMethod::BBS || partitionMethod == PartitionMethod::B2S)
    {
        std::stringstream errorMessage;
        errorMessage << "Error in ComputePartialCharges::compute(): partitionMethod \""  << to_string(partitionMethod) << "\" invalid." << std::endl;
        errorMessage << "Please check the documentation of the \"ComputePartialCharges\" job.";

        print_error(errorMessage.str());

        std::exit(1);
    }


    // Print partition method
    std::cout << "Volume partition method: " << to_string(partitionMethod) << std::endl <<std::endl;


    std::vector<double> charges;

    std::ifstream gridFile(gridFileName);
    if (!gridFile)
    {
        print_error("Error in ComputePartialCharges::compute(): could not read file " + gridFileName + ".");
        std::exit(1);
    }

    Grid grid(gridFile, PeriodicTable());
    gridFile.close();

    if (partitionMethod == PartitionMethod::BECKE)
    {
        Becke B(grid);
        B.partial_charge(grid);
        B.printCharges();

        charges = B.get_partial_charge();
    }
    else
    {
        GridCP gridcp;
        gridcp.buildBasins(grid, partitionMethod);
        gridcp.computeIntegrals(grid);
        gridcp.printCriticalPoints();

        charges = gridcp.computeAIMCharges(grid);
    }

    return charges;
}


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

void ComputePartialCharges::run()
{
    // Read grid files names
    std::vector<std::string> gridFilesNames;
    readGridFilesNames(gridFilesNames);


    // Read partition method
    PartitionMethod partitionMethod;
    readPartitionMethod(partitionMethod);

    
    // Compute partial charges
    std::vector<double> charges = compute(gridFilesNames[0], partitionMethod);
}
