#ifndef CDFTT_COMPUTECONDENSEDLINEARRESPONSE_HPP_INCLUDED
#define CDFTT_COMPUTECONDENSEDLINEARRESPONSE_HPP_INCLUDED

#include <JobControl/Job.h>
#include <Utils/Enums.hpp>


/** @brief ComputeCondensedLinearResponse job class.
 *
 *  This job is run when the "runType" parameter in the input file is set to "ComputeCondensedLinearResponse".
 */
class ComputeCondensedLinearResponse : public Job
{
    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTOR
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Constructor for the ComputeCondensedLinearResponse class.
         */
        ComputeCondensedLinearResponse(const std::string& inputFileName);

        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Runs the "ComputeCondensedLinearResponse" job: computes the linear responsed condensed on atoms.
         */
        void run() override;
};

#endif // CDFTT_COMPUTECONDENSEDLINEARRESPONSE_HPP_INCLUDED