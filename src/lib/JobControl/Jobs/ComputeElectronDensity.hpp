#ifndef CDFTT_COMPUTEELECTRONDENSITY_HPP_INCLUDED
#define CDFTT_COMPUTEELECTRONDENSITY_HPP_INCLUDED

#include <JobControl/Job.h>
#include <Utils/Enums.hpp>


/** @brief Compute electron density job class.
 * 
 * This job is run when the "runType" parameter in the input file is set to "ComputeElectronDensity".
 * It computes the electron density for the given system.
 */
class ComputeElectronDensity : public Job
{   
    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTOR
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Constructor for the ComputeElectronDensity class.
         */
        ComputeElectronDensity(const std::string& inputFileName);


        //----------------------------------------------------------------------------------------------------//
        // STATIC METHODS
        //----------------------------------------------------------------------------------------------------//

        static void compute();
        

        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Computes the electron density for the given system.
         */
        void run() override;
};

#endif // CDFTT_COMPUTEELECTRONDENSITY_HPP_INCLUDED