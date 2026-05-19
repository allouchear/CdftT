#ifndef CDFTT_MAKEORBITALSCUBE_HPP_INCLUDED
#define CDFTT_MAKEORBITALSCUBE_HPP_INCLUDED

#include <string>

#include <JobControl/Job.h>
#include <Utils/Enums.hpp>


/** @brief MakeOrbitalsCube job class.
 *
 *  This job is run when the "runType" parameter in the input file is set to "MakeOrbitalsCube".
 */
class MakeOrbitalsCube : public Job
{
    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTOR
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Constructor for the MakeOrbitalsCube class.
         */
        MakeOrbitalsCube(const std::string& inputFileName);


        //----------------------------------------------------------------------------------------------------//
        // STATIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Creates an orbitals cube file from an analytic file.
         * 
         * @param[in] analyticFileName Name of the analytic file to read the orbitals from.
         * @param[in] outputFileName Name of the output cube file to create.
         * @param[in] orbitalType Type of orbitals to include in the cube (All/Occupied/Virtual/Custom).
         * @param[in] spinType Spin type of orbitals to include in the cube (Alpha/Beta/Alpha-Beta).
         * @param[in] gridSize Grid size (Coarse/Medium/Fine/Custom).
         * @param[in] customSizeData Custom size parameters (used for CUSTOM grid size).
         * @param[in] showProgress Whether to display a progress bar during cube creation.
         */
        static void createCube(const std::string& analyticFileName, const std::string& outputFileName, OrbitalType orbitalType, SpinType spinType, GridSize gridSize, CustomSizeData customSizeData, bool showProgress = false, std::vector<int> orbitalsNumbers = { 0 }, std::vector<SpinType> spinList = { SpinType::ALPHA });


        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Runs the "MakeOrbitalsCube" job: creates an orbitals cube file.
         */
        void run() override;
};

#endif // CDFTT_MAKEORBITALSCUBE_HPP_INCLUDED