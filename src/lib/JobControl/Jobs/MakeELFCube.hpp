#ifndef CDFTT_MAKEELFCUBE_HPP_INCLUDED
#define CDFTT_MAKEELFCUBE_HPP_INCLUDED

#include <string>

#include <JobControl/Job.h>
#include <Utils/Enums.hpp>


/** @brief MakeELFCube job class.
 *
 *  This job is run when the "runType" parameter in the input file is set to "MakeELFCube".
 */
class MakeELFCube : public Job
{
    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTOR
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Constructor for the MakeELFCube class.
         */
        MakeELFCube(const std::string& inputFileName);


        //----------------------------------------------------------------------------------------------------//
        // STATIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Creates an Electron Localisation Function cube file from an analytic file.
         * 
         * @param[in] analyticFileName Name of the analytic file to read the orbitals from.
         * @param[in] outputFileName Name of the output cube file to create.
         * @param[in] elfMethod ELF method to use.
         * @param[in] gridSize Grid size (Coarse/Medium/Fine/Custom).
         * @param[in] customSizeData Custom size parameters (used for CUSTOM grid size).
         * @param[in] showProgress Whether to display a progress bar during cube creation.
         */
        static void createCube(const std::string& analyticFileName, const std::string& outputFileName, ELFMethod elfMethod, GridSize gridSize, CustomSizeData customSizeData, bool showProgress = false);


        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Runs the "MakeELFCube" job: creates an ELF cube file.
         */
        void run() override;
};

#endif // CDFTT_MAKEELFCUBE_HPP_INCLUDED