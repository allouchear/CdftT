#ifndef CDFTT_COMPUTEDESCRIPTORS_HPP_INCLUDED
#define CDFTT_COMPUTEDESCRIPTORS_HPP_INCLUDED

#include <string>
#include <vector>

#include <JobControl/Job.h>
#include <Utils/Enums.hpp>


/** @brief ComputeDescriptors job class.
 *
 *  This job is run when the "runType" parameter in the input file is set to "ComputeDescriptors".
 */
class ComputeDescriptors : public Job
{
    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTOR
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Constructor for the ComputeDescriptors class.
         * 
         * @param inputFileName[in] The name of the input file to be processed.
         */
        ComputeDescriptors(const std::string& inputFileName);


        //----------------------------------------------------------------------------------------------------//
        // STATIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Computes descriptors from three grid files using ionization energy and electron affinity.
         *
         * @param[in] gridFileName1 Name of the electrophilic grid file.
         * @param[in] gridFileName2 Name of the nucleophilic grid file.
         * @param[in] gridFileName3 Name of the radical grid file.
         * @param[in] ionisationEnergy Ionization energy.
         * @param[in] electronAffinity Electron affinity.
         * @param[in] partitionMethod Selected partition method.
         * @return Descriptors object containing computed descriptors.
         */
        static Descriptors computeWithIonisationEnergyAndElectronAffinity(const std::string& gridFileName1, const std::string& gridFileName2, const std::string& gridFileName3, double ionisationEnergy, double electronAffinity, PartitionMethod partitionMethod);

        /**
         * @brief Computes descriptors from three grid files using these file energies.
         *
         * @param[in] gridFileName1 Name of the first grid file.
         * @param[in] gridFileName2 Name of the second grid file.
         * @param[in] gridFileName3 Name of the third grid file.
         * @param[in] energies Vector of energies respectively corresponding to the files.
         * @param[in] partitionMethod Selected partition method.
         * @return Descriptors object containing computed descriptors.
         */
        static Descriptors computeWithEnergies(const std::string& gridFileName1, const std::string& gridFileName2, const std::string& gridFileName3, const std::vector<double>& energies, PartitionMethod partitionMethod);

        /**
         * @brief Computes descriptors with the finite-difference Becke method from analytic files.
         *
         * @param[in] ANAFileName1 First analytic input file.
         * @param[in] ANAFileName2 Second analytic input file.
         * @param[in] ANAFileName3 Third analytic input file.
         * @param[in] kmax Fuzzyness of the Voronoi polyhedrons (default 3).
         * @param[in] lebedev_order Lebedev order for angular quadrature (default 41).
         * @param[in] radial_grid_factor Radial grid multiplicative factor (default 5).
         * 
         * @return Descriptors object containing computed descriptors.
         */
        static Descriptors computeWithFiniteDifferenceBecke(const std::string& ANAFileName1, const std::string& ANAFileName2, const std::string& ANAFileName3, int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5);


        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Runs the "ComputeDescriptors" job: computes the descriptors.
         */
        void run() override;
};

#endif // CDFTT_COMPUTEDESCRIPTORS_HPP_INCLUDED