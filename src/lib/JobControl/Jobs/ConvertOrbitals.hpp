#ifndef CDFTT_CONVERTORBITALS_HPP_INCLUDED
#define CDFTT_CONVERTORBITALS_HPP_INCLUDED

#include <string>

#include <JobControl/Job.h>


/** @brief ConvertOrbitals job class.
 *
 *  This job is run when the "runType" parameter in the input file is set to "ConvertOrbitals".
 */
class ConvertOrbitals : public Job
{
    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTOR
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Constructor for the ConvertOrbitals class.
         */
        ConvertOrbitals(const std::string& inputFileName);


        //----------------------------------------------------------------------------------------------------//
        // STATIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Converts orbitals from one format to another.
         *
         * @param[in] inputFileName Input filename for the orbitals to be converted.
         * @param[in] outputFileName Output filename for the converted orbitals.
         */
        static void convert(const std::string& inputFileName, const std::string& outputFileName);


        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Runs the "ConvertOrbitals" job: converts the analytic files.
         */
        void run() override;
};

#endif // CDFTT_CONVERTORBITALS_HPP_INCLUDED