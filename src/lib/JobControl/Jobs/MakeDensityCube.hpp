#ifndef CDFTT_MAKEDENSITYCUBE_HPP_INCLUDED
#define CDFTT_MAKEDENSITYCUBE_HPP_INCLUDED

#include <string>

#include <JobControl/Job.h>
#include <Utils/Enums.hpp>


/** @brief MakeDensityCube job class.
 *
 *  This job is run when the "runType" parameter in the input file is set to "MakeDensityCube".
 */
class MakeDensityCube : public Job
{
    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTOR
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Constructor for the MakeDensityCube class.
         */
        MakeDensityCube(const std::string& inputFileName);


        //----------------------------------------------------------------------------------------------------//
        // STATIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Creates a density cube file from an analytic file.
         * 
         * @param[in] analyticFileName Name of the analytic file to read the orbitals from.
         * @param[in] outputFileName Name of the output cube file to create.
         * @param[in] gridSize Grid size (Coarse/Medium/Fine/Custom).
         * @param[in] customSizeData Custom size parameters (used for CUSTOM grid size).
         * @param[in] showProgress Whether to display a progress bar during cube creation.
         */
        static void createCube(const std::string& analyticFileName, const std::string& outputFileName, GridSize gridSize, CustomSizeData customSizeData, bool showProgress = false);


        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Runs the "MakeDensityCube" job: creates a density cube file.
         */
        void run() override;
};

#endif // CDFTT_MAKEDENSITYCUBE_HPP_INCLUDED