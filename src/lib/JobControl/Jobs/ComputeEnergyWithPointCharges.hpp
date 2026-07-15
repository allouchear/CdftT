#ifndef CDFTT_COMPUTEENERGYWITHPOINTCHARGES_HPP_INCLUDED
#define CDFTT_COMPUTEENERGYWITHPOINTCHARGES_HPP_INCLUDED

#include <iostream>
#include <string>
#include <vector>

#include <Common/Atom.h>
#include <JobControl/Job.h>
#include <Orbitals/ExcitedState.hpp>


/**
 * @brief ComputeEnergyWithPointCharges job class.
 * 
 * This job is run when the "runType" parameter in the input file is set to "ComputeEnergyWithPointCharges".
 */
class ComputeEnergyWithPointCharges : public Job
{
    private:
        //----------------------------------------------------------------------------------------------------//
        // TYPE DEFINITIONS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Represents a point charge.
         */
        struct PointCharge
        {
            /* @brief The charge of the point charge. */
            double charge;

            /* @brief The position of the point charge. */
            std::array<double, 3> position;

            /* @brief A description of the point charge (its charge and position, and the name of the atom it is placed on if this is the case). */
            std::string description;
        };

        /** @brief Represents a run with a list of PointCharges. */
        typedef std::vector<PointCharge> Run;
        

        //----------------------------------------------------------------------------------------------------//
        // PRIVATE METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Computes < psi_i | H_0 | psi_j >, < psi_i | H_1 | psi_j > and < psi_i | H | psi_j > (where H = H_0 + H_1) for a set of excited states and one or many point charge(s).
         * 
         * @param[in] states Vector of excited states for which the computations will be performed.
         * @param[in] chargesNucleiContributions Values of the < psi_i | V_ion/nuclei | psi_i > contributions (in the order of the charges).
         * @param[in] ionicMatrixes Matrixes of the < phi_i | V_ion/electrons | phi_j > contributions for the point charges.
         * @param[out] psi_i_H_0_psi_j Output lower triangular matrix where the computed < psi_i | H_0 | psi_j > values will be stored.
         * @param[out] psi_i_H_1_psi_j Output lower triangular matrix where the computed < psi_i | H_1 | psi_j > values will be stored.
         * @param[out] psi_i_H_psi_j Output lower triangular matrix where the computed < psi_i | H | psi_j > values will be stored.
         * @param[in,out] outputStream Stream where information will be logged during the computation.
         * @param[in] verbose Verbosity level for outputting intermediate values during computation (default 0).
         */
        void computeHamiltonianMatrixes(const std::vector<ExcitedState>& excitedStates, const std::vector<double>& chargesNucleiContributions, const std::vector<std::vector<std::vector<std::vector<double>>>>& ionicMatrixes, std::vector<std::vector<double>>& psi_i_H_0_psi_j, std::vector<std::vector<double>>& psi_i_H_1_psi_j, std::vector<std::vector<double>>& psi_i_H_psi_j, std::ostream& outputStream, int verbose = 0);

        /**
         * @brief Computes and prints the results for the variational approach.
         *
         * @param[in] states Vector of excited states of the system.
         * @param[in] psi_i_H_0_psi_j Matrix of the < psi_i | H_0 | psi_j > values.
         * @param[in] psi_i_H_1_psi_j Matrix of the < psi_i | H_1 | psi_j > values.
         * @param[in] outputFilePrefix Output filename prefix for saving results.
         * @param[in,out] outputStream Output stream for printing results.
         * @param[in] verbose Verbosity level for outputting intermediate values during computation (default 0).
         */
        void computeResults_perturbative(const std::vector<ExcitedState>& states, const std::vector<std::vector<double>>& psi_i_H_0_psi_j, const std::vector<std::vector<double>>& psi_i_H_1_psi_j, const std::string& outputFilePrefix, std::ostream& outputStream, int verbose = 0);

        /**
         * @brief Computes and prints the results for the variational approach.
         *
         * @param[in] states Vector of excited states of the system.
         * @param[in] psi_i_H_psi_j Matrix of the < psi_i | H | psi_j > values.
         * @param[in] outputFilePrefix Output filename prefix for saving results.
         * @param[in,out] outputStream Output stream for printing results.
         * @param[in] verbose Verbosity level for outputting intermediate values during computation (default 0).
         */
        void computeResults_variational(const std::vector<ExcitedState>& states, const std::vector<std::vector<double>>& psi_i_H_psi_j, const std::string& outputFilePrefix, std::ostream& outputStream, int verbose = 0);

        /**
         * @brief TODO
         */
        void computeResults_linearResponse(const std::vector<std::vector<double>>& lrfMatrixEigenvalues, const std::vector<std::vector<std::vector<double>>>& ionicPotentialVectors, const std::vector<std::vector<double>>& psi_i_H_0_psi_j, const std::vector<std::vector<double>>& psi_i_H_1_psi_j, const std::string& outputFilePrefix, std::ostream& outputStream, int verbose = 0);

        /**
         * @brief Prints the results.
         *
         * @param[in] methods Vector of methods used to compute the energy with point charges.
         * @param[in] orbitals Orbital information for the system.
         * @param[in] states Vector of excited states of the system.
         * @param[in] ionicPotentialMatrixes Ionic potential matrixes (values of the < psi_i | V_ion/electrons | psi_i >). Used for perturbative and variational approaches only.
         * @param[in] chargeNucleiContributions Values of the < psi_i | V_ion/nuclei | psi_i > contributions. Used for perturbative and variational approaches only.
         * @param[in] lrfMatrixEigenvalues Eigenvalues of the LRF matrix. Used for the linear response approach only.
         * @param[in] ionicPotentialVectors Ionic potential vectors. Used for the linear response approach only.
         * @param[in] runs Vector of runs, where each run is a vector of PointCharges and their descriptions.
         * @param[in] outputPrefix Output filename prefix for saving results.
         * @param[in,out] outputStream Output stream for printing results.
         * @param[in] verbose Verbosity level for outputting intermediate values during computation (default 0).
         */
        void printResults(const std::vector<EnergyPointChargeMethod>& methods, const std::vector<ExcitedState>& states, const std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& ionicPotentialMatrixes, const std::vector<std::vector<double>>& chargeNucleiContributions, const std::vector<std::vector<double>>& lrfMatrixEigenvalues, const std::vector<std::vector<std::vector<std::vector<double>>>>& ionicPotentialVectors, const std::vector<Run>& runs, const std::string& outputPrefix, std::ostream& outputStream, int verbose = 0);

        /**
         * TODO
         */
        void readChargesAndPositions(std::vector<Run>& runs, const std::vector<Atom>& atoms, std::ostream& outputStream);

        /**
         * TODO
         */
        void useLinearResponseApproach(Orbitals& orbitals, const std::vector<Run>& runs, const std::string& analyticFileName, const std::string& outputPrefix, const std::vector<int>& beckeParams, bool savePseudoOrbitals, bool showProgress, int verbose, std::ostream& outputStream);


    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTOR
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Constructor for the ComputeEnergyWithPointCharges class.
         * 
         * @param inputFileName[in] The name of the input file to be processed.
         */
        ComputeEnergyWithPointCharges(const std::string& inputFileName);


        //----------------------------------------------------------------------------------------------------//
        // STATIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Computes the contributions to the energy of the interactions between point charges and nuclei, i.e. the < psi_i | V_ion/nuclei | psi_i > values.
         * 
         * @param[in] atoms Vector of atoms in the system.
         * @param[out] chargeNucleiContributions Output vector where the computed contributions will be stored. The first dimension corresponds to the run index, and the second dimension corresponds to the PointCharge index in the run.
         * @param[in] runs Vector of runs, where each run is a vector of PointCharges and their descriptions.
         * @param[in] nuclearCutoff Distance cutoff under which the contribution of a nucleus is not taken into account (to avoid divergences when the point charge is very close to a nucleus).
         */
        static void computeChargeNucleiContributions(const std::vector<Atom>& atoms, std::vector<std::vector<double>>& chargeNucleiContributions, const std::vector<Run>& runs, double nuclearCutoff);

        /**
         * @brief Computes the contributions to the energy of the interactions between point charges and electrons, i.e. the < psi_i | V_ion/electrons | psi_i > values.
         *
         * @param[in] becke Becke instance used for the computation.
         * @param[out] ionicPotentialMatrixes Output matrix where the computed contributions will be stored. The first dimension corresponds to the run index, and the second dimension corresponds to the PointCharge index in the run.
         * @param[in] runs Vector of runs, where each run is a vector of PointCharges and their descriptions.
         * @param[in] kmax Fuzzyness of the Voronoi polyhedrons (default 3).
         * @param[in] lebedev_order Lebedev order for angular quadrature (default 41).
         * @param[in] radial_grid_factor Radial grid multiplicative factor (default 5).
         */
        static void computeIonicPotentialMatrixesFromBecke(Becke& becke, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& ionicPotentialMatrixes, const std::vector<Run>& runs, int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5);

        /**
         * @brief Computes the contributions to the energy of the interactions between point charges and electrons, i.e. the < psi_i | V_ion/electrons | psi_i > values.
         * 
         * @param[in] orbitals Orbitals instance used for the computation.
         * @param[in] grid Grid instance used for the computation.
         * @param[out] ionicPotentialMatrixes Output matrix where the computed contributions will be stored. The first dimension corresponds to the run index, and the second dimension corresponds to the PointCharge index in the run.
         * @param[in] runs Vector of runs, where each run is a vector of PointCharges and their descriptions.
         */
        static void computeIonicPotentialMatrixesFromGrid(const Orbitals& orbitals, Grid& grid, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& ionicPotentialMatrixes, const std::vector<Run>& runs);

        /**
         * @brief Computes the contributions to the energy of the interactions between point charges and electrons, i.e. the < psi_i | V_ion/electrons | psi_i > values.
         * 
         * @param[in] orbitals Orbitals instance used for the computation.
         * @param[out] ionicPotentialMatrixes Output matrix where the computed contributions will be stored. The first dimension corresponds to the run index, and the second dimension corresponds to the PointCharge index in the run.
         * @param[in] runs Vector of runs, where each run is a vector of PointCharges and their descriptions.
         */
        static void computeIonicPotentialMatrixesFromOrbitals(Orbitals& orbitals, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& ionicPotentialMatrixes, const std::vector<Run>& runs);

        /**
         * @brief Computes the contributions to the energy of the interactions between point charges and electrons, using a Becke-grid-based approach.
         * 
         * @param[in] becke Becke instance used for the computation.
         * @param[out] ionicPotentialVectors Output matrix where the computed contributions will be stored. The first dimension corresponds to the run index, the second dimension corresponds to the point charge index in the run, and the third dimension corresponds to the spin (0 for alpha, 1 for beta).
         * @param[in] runs Vector of runs, where each run is a vector of PointCharges and their descriptions.
         * @param[in] kmax Fuzzyness of the Voronoi polyhedrons (default 3).
         * @param[in] lebedev_order Lebedev order for angular quadrature (default 41).
         * @param[in] radial_grid_factor Radial grid multiplicative factor (default 5).
         */
        static void computeIonicPotentialVectorsFromBecke(Becke& becke, std::vector<std::vector<std::vector<std::vector<double>>>>& ionicPotentialVectors, const std::vector<Run>& runs, int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5);

        /**
         * @brief Computes the contributions to the energy of the interactions between point charges and electrons, using a grid-based approach.
         * 
         * @param[in] orbitals Orbitals instance used for the computation.
         * @param[in] grid Grid instance used for the computation.
         * @param[out] ionicPotentialVectors Output matrix where the computed contributions will be stored. The first dimension corresponds to the run index, the second dimension corresponds to the point charge index in the run, and the third dimension corresponds to the spin (0 for alpha, 1 for beta).
         * @param[in] runs Vector of runs, where each run is a vector of PointCharges and their descriptions.
         */
        static void computeIonicPotentialVectorsFromOrbitals(Orbitals& orbitals, std::vector<std::vector<std::vector<std::vector<double>>>>& ionicPotentialVectors, const std::vector<Run>& runs);

        /**
         * TODO
         */
        static void computeLinearResponseFunctionMatrix(const Orbitals& orbitals, const std::vector<std::vector<std::vector<std::vector<double>>>>& tripleOrbitalIntegralMatrix, std::vector<std::vector<std::vector<double>>>& lrfMatrix);


        //----------------------------------------------------------------------------------------------------//
        // PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Runs the "ComputeEnergyWithPointCharges" job: computes the energy of excited states in the presence of one or many point charge(s).
         */
        void run() override;


        /**
         * @brief Runs the job associated with the "ComputeLinearResponseWithPointCharges" runtype.
         */
        void run_computeLinearResponseWithPointCharges();
};

#endif // CDFTT_COMPUTEENERGYWITHPOINTCHARGES_HPP_INCLUDED
