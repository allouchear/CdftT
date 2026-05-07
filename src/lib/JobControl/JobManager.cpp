#include <JobControl/ComputeEnergyWithPointCharges.hpp>
#include <JobControl/Job.h>
#include <JobControl/JobManager.hpp>
#include <JobControl/Help.hpp>
#include <Utils/Enums.hpp>

//----------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------//

JobManager::JobManager(const std::string& inputFileName):
    _timer()
{
    handleInputFile(inputFileName);
}

JobManager::JobManager(const std::vector<std::string>& inputFileNames):
    _timer()
{
    for (const std::string& inputFileName : inputFileNames)
    {
        handleInputFile(inputFileName);
    }
}


//----------------------------------------------------------------------------------------------------//
// PRIVATE METHODS
//----------------------------------------------------------------------------------------------------//

void JobManager::createJob(const std::string& inputFileName, RunType runtype)
{
    switch (runtype)
    {
        case RunType::COMPUTE_CONDENSED_LINEAR_RESPONSE:
        {
            _jobs.push_back(std::make_unique<Job>(inputFileName));
            break;
        }
        case RunType::COMPUTE_DESCRIPTORS:
        {
            _jobs.push_back(std::make_unique<Job>(inputFileName));
            break;
        }
        case RunType::COMPUTE_ENERGY_WITH_POINT_CHARGES:
        {
            _jobs.push_back(std::make_unique<ComputeEnergyWithPointCharges>(inputFileName));
            break;
        }
        case RunType::COMPUTE_GRID_DIFFERENCE:
        {
            _jobs.push_back(std::make_unique<Job>(inputFileName));
            break;
        }
        case RunType::COMPUTE_INTEGRALS:
        {
            _jobs.push_back(std::make_unique<Job>(inputFileName));
            break;
        }
        case RunType::COMPUTE_LINEAR_RESPONSE_WITH_POINT_CHARGES:
        {
            _jobs.push_back(std::make_unique<Job>(inputFileName));
            break;
        }
        case RunType::COMPUTE_PARTIAL_CHARGES:
        {
            _jobs.push_back(std::make_unique<Job>(inputFileName));
            break;
        }
        case RunType::CONVERT_ORBITALS:
        {
            _jobs.push_back(std::make_unique<Job>(inputFileName));
            break;
        }
        case RunType::HELP:
        {
            _jobs.push_back(std::make_unique<Help>());
            break;
        }
        case RunType::LAMBDA_DIAGNOSTIC:
        {
            _jobs.push_back(std::make_unique<Job>(inputFileName));
            break;
        }
        case RunType::MAKE_DENSITY_CUBE:
        {
            _jobs.push_back(std::make_unique<Job>(inputFileName));
            break;
        }
        case RunType::MAKE_ORBITALS_CUBE:
        {
            _jobs.push_back(std::make_unique<Job>(inputFileName));
            break;
        }
        case RunType::MAKE_ELF_CUBE:
        {
            _jobs.push_back(std::make_unique<Job>(inputFileName));
            break;
        }
        default:
        {
            // This should never happen because the run type is checked before calling this function, but we put it just in case
            std::cerr << "Error in JobManager::createJob(): unknown run type. No job will be created for input file \"" << inputFileName << "\"." << std::endl;
        }
    }
}

void JobManager::handleInputFile(const std::string& inputFileName)
{
    std::ifstream inputFile(inputFileName);
    if (!inputFile)
    {
        std::cerr << "Error in JobManager::handleInputFile(): could not open input file \"" << inputFileName << "\". This file will be skipped." << std::endl;
        std::cerr << "Please check that the file exists and that the program has the necessary permissions to read it." << std::endl;
    }
    else
    {
        std::string strRunType;
        bool read = readOneString(inputFile, "RunType", strRunType);
        inputFile.close();

        if (!read)
        {
            std::cerr << "Error in JobManager::handleInputFile(): the \"RunType\" parameter is not specified in the input file \"" << inputFileName << "\". This file will be skipped." << std::endl;
            std::cerr << "Please check the documentation for the format of the input file." << std::endl;
        }
        else
        {
            RunType runType = runType_from_string(strRunType);

            if (runType != RunType::UNKNOWN)
            {
                createJob(inputFileName, runType);
            }
            else
            {
                std::cerr << "Error in JobManager::handleInputFile(): the \"RunType\" parameter value \"" << strRunType << "\" specified in the input file \"" << inputFileName << "\" is unknown. This file will be skipped." << std::endl;
                std::cerr << "Please check the documentation and the value of the \"RunType\" parameter." << std::endl;
            }
        }
    }
}


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

void JobManager::runJobs()
{
    for (auto& job : _jobs)
    {
        job->run();
        std::cout << std::defaultfloat << "Time in ms: " << _timer.get() << std::endl;
    }
}