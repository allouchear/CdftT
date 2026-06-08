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

        /**
         * @brief Computes electronic densities for a set of excited states and saves them in .cube files.
         * 
         * @param[in] states Vector of excited states (index 0 corresponds to the Ground State)
         * @param[in] excitedStatesNumbers List of states for which to compute the density.
         * @param[in] orbitals Orbitals of the system.
         * @param[in] grid Grid on which the densities will be computed.
         * @param[in] rdmMethod Method to use for the reduced density matrix computation.
         * @param[in] excludedOrbitalsNumbers List of orbital numbers to exclude from the density computation.
         * @param[in] outputPrefix Output filename prefix for saving results.
         * @param[in] saveRDM Flag indicating whether the computed reduced density matrix should be saved in .cdftt files.
         * @param[in] outputPrecision Number of digits to use for numerical values in the output files.
         * @param[in] verbose Verbosity level for outputting intermediate values during computation (default 0).
         * @param[in] logOutputStream Output stream for logging information during the computation.
         * @param[in] showProgress Whether to show a progress bar during the computation.
         */
        static void computeStateDensities(const std::vector<ExcitedState>& states, std::vector<int>& excitedStatesNumbers, const Orbitals& orbitals, Grid& grid, const RDMMethod rdmMethod, const std::vector<int>& excludedOrbitalsNumbers, const std::string& outputPrefix, bool saveRDM, int outputPrecision, int verbose, std::ostream& logOutputStream, bool showProgress);

        /**
         * @brief Computes electronic densities for a set of transitions between excited states and saves them in .cube files.
         * 
         * @param[in] states Vector of excited states for which the transition densities will be computed.
         * @param[in] transitionDensities List of pairs of excited state numbers for which the transition density will be computed.
         * @param[in] orbitals Orbitals of the system.
         * @param[in] grid Grid on which the densities will be computed.
         * @param[in] rdmMethod Method to use for the reduced density matrix computation.
         * @param[in] excludedOrbitalsNumbers List of orbital numbers to exclude from the density computation.
         * @param[in] outputPrefix Output filename prefix for saving results.
         * @param[in] saveRDM Flag indicating whether the computed reduced density matrix should be saved in .cdftt files.
         * @param[in] outputPrecision Number of digits to use for numerical values in the output files.
         * @param[in] verbose Verbosity level for outputting intermediate values during computation (default 0).
         * @param[in] logOutputStream Output stream for logging information during the computation.
         * @param[in] showProgress Whether to show a progress bar during the computation.
         */
        static void computeTransitionDensities(const std::vector<ExcitedState>& states, const std::vector<std::array<int, 2>>& transitionDensities, const Orbitals& orbitals, Grid& grid, const RDMMethod rdmMethod, const std::vector<int>& excludedOrbitalsNumbers, const std::string& outputPrefix, bool saveRDM, int outputPrecision, int verbose, std::ostream& logOutputStream, bool showProgress);
        

        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Computes the electron density for the given system.
         */
        void run() override;
};

#endif // CDFTT_COMPUTEELECTRONDENSITY_HPP_INCLUDED