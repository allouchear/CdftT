#ifndef CDFTT_JOBMANAGER_HPP_INCLUDED
#define CDFTT_JOBMANAGER_HPP_INCLUDED

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include <Utils/Enums.hpp>
#include <JobControl/Job.h>
#include <Utils/Timer.hpp>


/** @brief Job manager class.
 * 
 * Manages the execution of jobs based on program input. Reads the "RunType" parameter from the input file(s) and runs the corresponding job(s).
 */
class JobManager
{
    private:
        /** @brief The timer for measuring execution time. */
        Timer _timer;

        /** @brief The vector of jobs to be executed. */
        std::vector<std::unique_ptr<Job>> _jobs;


        //----------------------------------------------------------------------------------------------------//
        // PRIVATE METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Creates a job based on the specified run type.
         * 
         * This function MUST be modified if a new Job is created!
         * 
         * @param[in] inputFileName The name of the input file associated with the job to create.
         * @param[in] runType The run type for which to create a job.
         */
        void createJob(const std::string& inputFileName, RunType runtype);

        /**
         * @brief Handles the input file: checks if it can be opened and creates the job based on the "RunType" parameter.
         * 
         * @param[in] inputFileName The name of the input file to handle.
         */
        void handleInputFile(const std::string& inputFileName);


    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTORS
        //----------------------------------------------------------------------------------------------------//
        
        /**
         * @brief Constructor for the JobManager class.
         * 
         * @param inputFileName[in] The name of the input file to be processed.
         */
        JobManager(const std::string& inputFileName);

        /** @brief Constructor for the JobManager class.
         * 
         * @param inputFileNames[in] A vector of input file names to be processed.
         */
        JobManager(const std::vector<std::string>& inputFileNames);


        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Runs all the jobs in the job manager.
         *
         * This function iterates through all the jobs stored in the job manager and calls their run() method to execute them.
         */
        void runJobs();
};


#endif /* CDFTT_JOBMANAGER_HPP_INCLUDED */