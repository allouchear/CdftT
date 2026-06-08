#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <JobControl/Job.h>
#include <Utils/FCHK.h>
#include <Utils/LOG.h>
#include <Utils/MOLDENGAB.h>
#include <Utils/WFX.h>


template<typename T, typename U>
U Job::computeOrbitalsOrBecke(const std::string& analyticFileName)
{
    std::ifstream analyticFile(analyticFileName);
    if (!analyticFile)
    {
        print_error("Error in Job::computeOrbitalsOrBecke(): could not read file " + analyticFileName + ".");
        std::exit(1);
    }

    T analyticFileParser(analyticFile);
    analyticFile.close();

    U analyticObject(analyticFileParser, Binomial(100), PeriodicTable());
    
    return analyticObject;
}

template<typename T>
void Job::computeOrbitalsOrBecke(T& analyticObject, const std::string& analyticFileName)
{
    if (analyticFileName.find(".wfx") != std::string::npos)
    {
        std::cout << "Reading data from " << analyticFileName << "... Please wait." << std::endl;
        analyticObject = computeOrbitalsOrBecke<WFX, T>(analyticFileName);
    }
    else if (analyticFileName.find(".fchk") != std::string::npos)
    {
        std::cout << "Reading data from " << analyticFileName << "... Please wait." << std::endl;
        analyticObject = computeOrbitalsOrBecke<FCHK, T>(analyticFileName);
    }
    else if (analyticFileName.find(".molden") != std::string::npos)
    {
        std::cout << "Reading data from " << analyticFileName << "... Please wait." << std::endl;
        analyticObject = computeOrbitalsOrBecke<MOLDENGAB, T>(analyticFileName);
    }
    else if (analyticFileName.find(".gab") != std::string::npos)
    {
        std::cout << "Reading data from " << analyticFileName << "... Please wait." << std::endl;
        analyticObject = computeOrbitalsOrBecke<MOLDENGAB, T>(analyticFileName);
    }
    else if (analyticFileName.find(".log") != std::string::npos)
    {
        std::cout << "Reading data from " << analyticFileName << "... Please wait." << std::endl;
        analyticObject = computeOrbitalsOrBecke<LOG, T>(analyticFileName);
    }
    else
    {
        std::stringstream errorMessage;
        errorMessage << "Error: invalid file format for file \"" << analyticFileName << "\"." << std::endl;
        errorMessage << "Please check the documentation and the \"AnalyticFiles\" parameter value in the provided input file.";

        print_error(errorMessage.str());

        std::exit(1);
    }
}