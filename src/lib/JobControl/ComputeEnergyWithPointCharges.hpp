
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
         */
        void computeResults(const std::vector<ExcitedState>& states, const std::vector<std::vector<double>>& psi_i_H_psi_j, const std::vector<std::vector<double>>& psi_i_HminusH0_psi_j, const std::string& outputFilePrefix, std::ostream& outputStream, int verbose);

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
         * @param[in] verbose Verbosity level.
         */
        void printResults(const std::vector<ExcitedState>& states, const std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& ionicPotentialMatrixes, const std::vector<std::vector<double>>& chargeNucleiContributions, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms, const std::vector<Atom>& atoms, const std::string& outputPrefix, std::ostream& outputStream, int verbose);


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
        // PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Runs the "ComputeEnergyWithPointCharges" job: computes the energy of excited states in the presence of one or many point charge(s).
         */
        void run() override;
};