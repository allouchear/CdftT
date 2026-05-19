#ifndef CDFTT_COMPUTEINTEGRALS_HPP_INCLUDED
#define CDFTT_COMPUTEINTEGRALS_HPP_INCLUDED

#include <string>
#include <vector>

#include <Cube/GridCP.h>
#include <JobControl/Job.h>


/** @brief ComputeIntegrals job class.
 *
 *  This job is run when the "runType" parameter in the input file is set to "ComputeIntegrals".
 */
class ComputeIntegrals : public Job
{
    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTOR
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Constructor for the ComputeIntegrals class.
         */
        ComputeIntegrals(const std::string& inputFileName);


        //----------------------------------------------------------------------------------------------------//
        // STATIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Builds basins from a grid using `partitionMethod`.
         *
         * @param[out] gridCP GridCP reference that will be populated with the constructed basins.
         * @param[in] gridFileName Path to the input grid (.cube file) used to build basins.
         * @param[in] partitionMethod Selected partition method.
         */
        static void buildBasins(GridCP& gridCP, const std::string& gridFileName, PartitionMethod partitionMethod);

        /**
         * @brief Builds basins by sign of grid values (BBS/B2S partition methods).
         *
         * @param[out] gridCP GridCP reference that will be populated with sign-based basins.
         * @param[in] gridFileName Path to the grid file used as input.
         * @param[in] cutoff Numerical cutoff below which values are considered zero.
         * @param[in] two If true, builds exactly two basins.
         */
        static void buildBasinsBySign(GridCP& gridCP, const std::string& gridFileName, double cutoff, bool two);

        /**
         * @brief Integrates the provided grids over the domain of a Critical Points grid (GridCP).
         *
         * @param[out] gridCP GridCP reference that stores the integration domain and the results.
         * @param[in] gridFilesNames Filenames of grids to integrate.
         * @param[in] partitionMethod The method used to partition the integration domain.
         * @param[in] cutoff The cutoff value for the integration.
         */
        static void computeLocalIntegrals(GridCP& gridCP, const std::vector<std::string>& gridFilesNames, PartitionMethod partitionMethod, double cutoff);


        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Runs the "ComputeIntegrals" job: computes the integrals.
         */
        void run() override;
};

#endif // CDFTT_COMPUTEINTEGRALS_HPP_INCLUDED