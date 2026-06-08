#ifndef CDFTT_COMPUTEGRIDDIFFERENCE_HPP_INCLUDED
#define CDFTT_COMPUTEGRIDDIFFERENCE_HPP_INCLUDED

#include <string>

#include <Cube/Grid.h>
#include <JobControl/Job.h>


/** @brief ComputeGridDifference job class.
 *
 *  This job is run when the "runType" parameter in the input file is set to "ComputeGridDifference".
 */
class ComputeGridDifference : public Job
{
    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTOR
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Constructor for the ComputeGridDifference class.
         * 
         * @param inputFileName[in] The name of the input file to be processed.
         */
        ComputeGridDifference(const std::string&  inputFileName);


        //----------------------------------------------------------------------------------------------------//
        // STATIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Computes the difference between two grids.
         *
         * @param[in] minuendGridFileName Left (minuend) input grid filename.
         * @param[in] subtrahendGridFileName Right (subtrahend) input grid filename.
         * @param[out] differenceGrid The grid containing the difference.
         */
        static void compute(const std::string& minuendGridFileName, const std::string& subtrahendGridFileName, Grid& differenceGrid);

        /**
         * @brief Computes the difference between two grids and saves the result to an output file.
         *
         * @param[in] minuendGridFileName Left (minuend) input grid filename.
         * @param[in] subtrahendGridFileName Right (subtrahend) input grid filename.
         * @param[in] outputGridFileName Output filename for the difference grid.
         */
        static void computeAndSave(const std::string& minuendGridFileName, const std::string& subtrahendGridFileName, const std::string& outputGridFileName);


        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Runs the "ComputeGridDifference" job: computes the grid difference.
         */
        void run() override;
};

#endif // CDFTT_COMPUTEGRIDDIFFERENCE_HPP_INCLUDED