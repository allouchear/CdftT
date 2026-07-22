#ifndef CDFTT_RUNLAMBDADIAGNOSTIC_HPP_INCLUDED
#define CDFTT_RUNLAMBDADIAGNOSTIC_HPP_INCLUDED

#include <string>

#include <JobControl/Job.h>
#include <Utils/Enums.hpp>


/** @brief RunLambdaDiagnostic job class.
 *
 *  This job is run when the "runType" parameter in the input file is set to "RunLambdaDiagnostic".
 */
class RunLambdaDiagnostic : public Job
{
    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTOR
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Constructor for the RunLambdaDiagnostic class.
         */
        RunLambdaDiagnostic(const std::string& inputFileName);


        //----------------------------------------------------------------------------------------------------//
        // STATIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Prints the lambda diagnostic computed with given electronic transitions and regular grid.
         *
         * @param[in] analyticFileName Input filename for the analytic files.
         * @param[in] transitionsFileName Input filename for the transitions file.
         * @param[in] maxNumberOfExcitedStates Maximum number of excited states to read in the transitions file.
         * @param[in] gridSize Grid size for the computation.
         * @param[in] customSizeData Custom size data for the computation.
         */
        static void print(const std::string& analyticFileName, const std::string& transitionsFileName, int maxNumberOfExcitedStates, GridSize gridSize, CustomSizeData customSizeData);


        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Runs the "RunLambdaDiagnostic" job: performs lambda diagnostic analysis.
         */
        void run() override;
};

#endif // CDFTT_RUNLAMBDADIAGNOSTIC_HPP_INCLUDED