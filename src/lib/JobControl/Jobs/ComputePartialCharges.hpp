#ifndef CDFTT_COMPUTEPARTIALCHARGES_HPP_INCLUDED
#define CDFTT_COMPUTEPARTIALCHARGES_HPP_INCLUDED

#include <string>

#include <JobControl/Job.h>


/** @brief ComputePartialCharges job class.
 *
 *  This job is run when the "runType" parameter in the input file is set to "ComputePartialCharges".
 */
class ComputePartialCharges : public Job
{
    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTOR
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Constructor for the ComputePartialCharges class.
         */
        ComputePartialCharges(const std::string& inputFileName);


        //----------------------------------------------------------------------------------------------------//
        // STATIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Computes partial atomic charges from a grid using the specified method.
         *
         * @param[in] gridFileName Input grid filename used for partitioning and integration.
         * @param[in] partitionMethod Partition method to use (AIM, VDD, Becke, ...).
         * @return Vector of computed partial charges.
         */
        static std::vector<double> compute(const std::string& gridFileName, PartitionMethod partitionMethod);


        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Runs the "ComputePartialCharges" job: computes the partial charges.
         */
        void run() override;
};

#endif // CDFTT_COMPUTEPARTIALCHARGES_HPP_INCLUDED