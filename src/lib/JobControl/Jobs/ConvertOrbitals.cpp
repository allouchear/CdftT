#include <regex>
#include <string>
#include <vector>

#include <JobControl/Job.h>
#include <JobControl/Jobs/ConvertOrbitals.hpp>
#include <Orbitals/Orbitals.h>


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//----------------------------------------------------------------------------------------------------//

ConvertOrbitals::ConvertOrbitals(const std::string& inputFileName):
    Job(inputFileName)
{ }


//----------------------------------------------------------------------------------------------------//
// STATIC METHODS
//----------------------------------------------------------------------------------------------------//


void ConvertOrbitals::convert(const std::string& inputFileName, const std::string& outputFileName)
{
    // Checking if the files have a different format
    std::regex fileExtensionRegex("\\.([a-zA-Z0-9]+)$");
    std::smatch matchInputFileExtension;
    std::smatch matchOutputFileExtension;
    std::string inputFileExtension;
    std::string outputFileExtension;
    std::regex_search(inputFileName, matchInputFileExtension, fileExtensionRegex);
    std::regex_search(outputFileName, matchOutputFileExtension, fileExtensionRegex);

    if (matchInputFileExtension[0] != matchOutputFileExtension[0])
    {
        // Loading orbitals and saving in the new format
        Orbitals o;
        computeOrbitalsOrBecke<Orbitals>(o, inputFileName);
        o.Save(outputFileName);
        std::cout << inputFileName << " has been converted to " << outputFileName << '.' << std::endl;
    }
    else
    {
        std::cout << "Input and output files have the same format (" << inputFileExtension << "). Nothing to be done." << std::endl;
    }
}


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

void ConvertOrbitals::run()
{
    // Read analytic files names
    std::vector<std::string> analyticFilesNames;
    readAnalyticFilesNames(analyticFilesNames);


    // Check number of analytic files names
    if (analyticFilesNames.size() != 2)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect number of analytic files names (two files expected)." << std::endl;
        errorMessage << "Please check the documentation and the number of files specified in the \"AnalyticFiles\" parameter in " << _inputFileName << '.';

        print_error(errorMessage.str());

        std::exit(1);
    }

    
    // Convert files
    convert(analyticFilesNames[0], analyticFilesNames[1]);
}