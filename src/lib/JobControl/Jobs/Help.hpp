#ifndef CDFTT_HELP_HPP_INCLUDED
#define CDFTT_HELP_HPP_INCLUDED

#include <JobControl/Job.h>
#include <Utils/Enums.hpp>


/** @brief Help job class.
 * 
 * This job is run when the "runType" parameter in the input file is set to "Help" or when the input file is not given.
 * It prints the list of available jobs and their descriptions.
 */
class Help : public Job
{
    private:
        //----------------------------------------------------------------------------------------------------//
        // STATIC FIELDS
        //----------------------------------------------------------------------------------------------------//

        /** @brief List of available job names (displayed in help). */
        static std::unordered_map<RunType, std::string> _s_availableJobs;
    
    
    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTOR
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Constructor for the Help class.
         */
        Help();


        //----------------------------------------------------------------------------------------------------//
        // STATIC METHODS
        //----------------------------------------------------------------------------------------------------//

        static void printAvailableJobs();
        

        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Prints the list of available jobs and their descriptions.
         */
        void run() override;
};

#endif // CDFTT_HELP_HPP_INCLUDED