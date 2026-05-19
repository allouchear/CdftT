#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <Common/PeriodicTable.h>
#include <Cube/Grid.h>
#include <JobControl/Job.h>
#include <JobControl/Jobs/ComputeGridDifference.hpp>
#include <Utils/Utils.h>


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//----------------------------------------------------------------------------------------------------//

ComputeGridDifference::ComputeGridDifference(const std::string& inputFileName):
    Job(inputFileName)
{ }


//----------------------------------------------------------------------------------------------------//
// STATIC METHODS
//----------------------------------------------------------------------------------------------------//

void ComputeGridDifference::compute(const std::string& minuendGridFileName, const std::string& subtrahendGridFileName, Grid& differenceGrid)
{
    // Compute grid difference
    std::cout << "Computing difference between " << minuendGridFileName << " and " << subtrahendGridFileName << '.' << std::endl;
    
    std::ifstream minuendFile(minuendGridFileName);
    if (!minuendFile)
    {
        print_error("Error in ComputeGridDifference::compute(): could not read file " + minuendGridFileName + ".");
        std::exit(1);
    }

    Grid minuendGrid(minuendFile, PeriodicTable());
    minuendFile.close();

    std::ifstream subtrahendFile(subtrahendGridFileName);
    if (!subtrahendFile)
    {
        print_error("Error in ComputeGridDifference::compute(): could not read file " + subtrahendGridFileName + ".");
        std::exit(1);
    }

    Grid subtrahendGrid(subtrahendFile, PeriodicTable());
    subtrahendFile.close();

    differenceGrid = minuendGrid - subtrahendGrid;
}

void ComputeGridDifference::computeAndSave(const std::string& minuendGridFileName, const std::string& subtrahendGridFileName, const std::string& outputGridFileName)
{
    Grid differenceGrid;
    compute(minuendGridFileName, subtrahendGridFileName, differenceGrid);

    std::ofstream outputGridFile(outputGridFileName, std::ios::out);
    if (!outputGridFile)
    {
        print_error("Error in ComputeGridDifference::computeAndSave(): failed to write to file " + outputGridFileName + ".");
        std::exit(1);
    }

    differenceGrid.save(outputGridFile);
    outputGridFile.close();

    std::cout << "Difference grid saved to file " << outputGridFileName << '.' << std::endl;
}


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

void ComputeGridDifference::run()
{
    // Read grid files names
    std::vector<std::string> gridFilesNames;
    readGridFilesNames(gridFilesNames);


    // Check number of grid files names
    if (gridFilesNames.size() != 3)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect number of grid files names (three files expected)." << std::endl;
        errorMessage << "Please check the documentation and the number of files specified in the \"GridFiles\" parameter in " << _inputFileName << '.';

        print_error(errorMessage.str());

        std::exit(1);
    }


    // Compute and save grid difference
    computeAndSave(gridFilesNames[0], gridFilesNames[1], gridFilesNames[2]);
}