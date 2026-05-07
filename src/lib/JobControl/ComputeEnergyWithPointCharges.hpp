#ifndef CDFTT_COMPUTEENERGYWITHPOINTCHARGES_HPP_INCLUDED
#define CDFTT_COMPUTEENERGYWITHPOINTCHARGES_HPP_INCLUDED

#include <iostream>
#include <string>
#include <vector>

#include <Common/Atom.h>
#include <JobControl/Job.h>
#include <Orbitals/ExcitedState.hpp>


class ComputeEnergyWithPointCharges : public Job
{
    private:
        //----------------------------------------------------------------------------------------------------//
        // PRIVATE METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Computes < psi_i | H | psi_j > and < psi_i | H-H0 | psi_j > for a set of excited states and one or many point charge(s).
         * 
         * @param[in] states Vector of excited states for which the computations will be performed.
         * @param[in] chargesNucleiContributions Values of the < psi_i | V_ion/nuclei | psi_i > contributions (in the order of the charges).
         * @param[in] ionicMatrixes Matrixes of the < phi_i | V_ion/electrons | phi_j > contributions for the point charges.
         * @param[out] psi_i_H_psi_j Output lower triangular matrix where the computed < psi_i | H | psi_j > values will be stored.
         * @param[out] psi_i_HminusH0_psi_j Output lower triangular matrix where the computed < psi_i | H - H0 | psi_j > values will be stored.
         * @param[in,out] outputStream Stream where information will be logged during the computation.
         * @param[in] verbose Verbosity level for outputting intermediate values during computation (default 0).
         */
        void computeHamiltonianMatrix(const std::vector<ExcitedState>& excitedStates, const std::vector<double>& chargesNucleiContributions, const std::vector<std::vector<std::vector<std::vector<double>>>>& ionicMatrixes, std::vector<std::vector<double>>& psi_i_H_psi_j, std::vector<std::vector<double>>& psi_i_HminusH0_psi_j, std::ostream& outputStream, int verbose = 0);

        /**
         * @brief Computes and prints the energy results from Hamiltonian matrix elements with point charges.
         *
         * @param[in] states Vector of excited states of the system.
         * @param[in] psi_i_H_psi_j Matrix of the < psi_i | H | psi_j > values.
         * @param[in] psi_i_HminusH0_psi_j Matrix of the < psi_i | H - H0 | psi_j > values.
         * @param[in] outputFilePrefix Output filename prefix for saving results.
         * @param[in,out] outputStream Output stream for printing results.
         * @param[in] verbose Verbosity level for outputting intermediate values during computation (default 0).
         */
        void computeResults(const std::vector<ExcitedState>& states, const std::vector<std::vector<double>>& psi_i_H_psi_j, const std::vector<std::vector<double>>& psi_i_HminusH0_psi_j, const std::string& outputFilePrefix, std::ostream& outputStream, int verbose = 0);

        /**
         * @brief Computes and prints the results with point charges.
         */
        void computeResultsLinearResponse(const std::vector<std::vector<double>>& eigenvalues, const std::vector<std::vector<std::vector<double>>>& ionicPotentialVectors, std::ostream& outputStream, int verbose);

        /**
         * @brief Prints the results.
         *
         * @param[in] states Vector of excited states of the system.
         * @param[in] ionicPotentialMatrixes Ionic potential matrixes (values of the < psi_i | V_ion/electrons | psi_i >).
         * @param[in] chargeNucleiContributions Values of the < psi_i | V_ion/nuclei | psi_i > contributions.
         * @param[in] charges Charges of the point charges.
         * @param[in] chargesPositions Positions of the point charges.
         * @param[in] loopOnAtoms Whether to loop on atoms (in that case, each charge uses every chargePosition).
         * @param[in] atoms Vector of atoms in the system (used if `loopOnAtoms` is true).
         * @param[in] outputPrefix Output filename prefix for saving results.
         * @param[in,out] outputStream Output stream for printing results.
         * @param[in] verbose Verbosity level for outputting intermediate values during computation (default 0).
         */
        void printResults(const std::vector<ExcitedState>& states, const std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& ionicPotentialMatrixes, const std::vector<std::vector<double>>& chargeNucleiContributions, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms, const std::vector<Atom>& atoms, const std::string& outputPrefix, std::ostream& outputStream, int verbose = 0);

        /**
         * @brief Prints the results for the Linear Response approach.
         *
         * @param[in] eigenvalues Eigenvalues of the system.
         * @param[in] ionicPotentialVectors Ionic potential vectors.
         * @param[in] charges Charges of the point charges.
         * @param[in] chargesPositions Positions of the point charges.
         * @param[in] loopOnAtoms Whether to loop on atoms (in that case, each charge uses every chargePosition).
         * @param[in,out] outputStream Output stream for printing results.
         * @param[in] verbose Verbosity level for outputting intermediate values during computation (default 0).
         */
        void printResultsLinearResponse(const std::vector<std::vector<double>>& eigenvalues, const std::vector<std::vector<std::vector<std::vector<double>>>>& ionicPotentialVectors, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms, const std::vector<Atom>& atoms, std::ostream& outputStream, int verbose = 0);


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
         * @param[out] chargeNucleiContributions Output vector where the computed contributions will be stored. The first dimension corresponds to the charges, and the second dimension corresponds to the charge positions.
         * @param[in] charges Vector of the charges of the point charges.
         * @param[in] chargesPositions Vector of the positions of the point charges.
         * @param[in] loopOnAtoms Whether to loop on atoms (in that case, each charge uses every chargePosition).
         * @param[in] nuclearCutoff Distance cutoff under which the contribution of a nucleus is not taken into account (to avoid divergences when the point charge is very close to a nucleus).
         */
        static void computeChargeNucleiContributions(const std::vector<Atom>& atoms, std::vector<std::vector<double>>& chargeNucleiContributions, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms, double nuclearCutoff);

        /**
         * @brief Computes the contributions to the energy of the interactions between point charges and electrons, i.e. the < psi_i | V_ion/electrons | psi_i > values.
         *
         * @param[in] becke Becke instance used for the computation.
         * @param[out] ionicPotentialMatrixes Output matrix where the computed contributions will be stored.
         * @param[in] charges Vector of the charges of the point charges.
         * @param[in] chargesPositions Vector of the positions of the point charges.
         * @param[in] loopOnAtoms Whether to loop on atoms (in that case, each charge uses every chargePosition).
         * @param[in] kmax Fuzzyness of the Voronoi polyhedrons (default 3).
         * @param[in] lebedev_order Lebedev order for angular quadrature (default 41).
         * @param[in] radial_grid_factor Radial grid multiplicative factor (default 5).
         */
        static void computeIonicPotentialMatrixesFromBecke(Becke& becke, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& ionicPotentialMatrixes, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms, int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5);

        /**
         * @brief Computes the contributions to the energy of the interactions between point charges and electrons, i.e. the < psi_i | V_ion/electrons | psi_i > values.
         * 
         * @param[in] orbitals Orbitals instance used for the computation.
         * @param[in] grid Grid instance used for the computation.
         * @param[out] ionicPotentialMatrixes Output matrix where the computed contributions will be stored.
         * @param[in] charges Vector of the charges of the point charges.
         * @param[in] chargesPositions Vector of the positions of the point charges.
         * @param[in] loopOnAtoms Whether to loop on atoms (in that case, each charge uses every chargePosition).
         */
        static void computeIonicPotentialMatrixesFromGrid(const Orbitals& orbitals, Grid& grid, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& ionicPotentialMatrixes, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms);

        /**
         * @brief Computes the contributions to the energy of the interactions between point charges and electrons, i.e. the < psi_i | V_ion/electrons | psi_i > values.
         * 
         * @param[in] orbitals Orbitals instance used for the computation.
         * @param[out] ionicPotentialMatrixes Output matrix where the computed contributions will be stored.
         * @param[in] charges Vector of the charges of the point charges.
         * @param[in] chargesPositions Vector of the positions of the point charges.
         * @param[in] loopOnAtoms Whether to loop on atoms (in that case, each charge uses every chargePosition).
         */
        static void computeIonicPotentialMatrixesFromOrbitals(Orbitals& orbitals, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& ionicPotentialMatrixes, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms);

        /**
         * @brief Computes the contributions to the energy of the interactions between point charges and electrons, using a .
         */
        static void computeIonicPotentialVectorsFromBecke(Becke& becke, std::vector<std::vector<std::vector<std::vector<double>>>>& ionicPotentialVectors, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms, int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5);

        /**
         * TODO
         */
        static void computeIonicPotentialVectorsFromOrbitals(Orbitals& orbitals, std::vector<std::vector<std::vector<std::vector<double>>>>& ionicPotentialVectors, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms);

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
