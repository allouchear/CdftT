#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include <Common/Atom.h>
#include <Common/Constants.h>
#include <JobControl/Job.h>
#include <JobControl/Jobs/ComputeEnergyWithPointCharges.hpp>
#include <Orbitals/ExcitedState.hpp>
#include <Orbitals/Orbitals.h>


//----------------------------------------------------------------------------------------------------//
// PRIVATE METHODS
//----------------------------------------------------------------------------------------------------//

void ComputeEnergyWithPointCharges::computeHamiltonianMatrix(const std::vector<ExcitedState>& states, const std::vector<double>& chargesNucleiContributions, const std::vector<std::vector<std::vector<std::vector<double>>>>& ionicMatrixes, std::vector<std::vector<double>>& psi_i_H_psi_j, std::vector<std::vector<double>>& psi_i_HminusH0_psi_j, std::ostream& outputStream, int verbose)
{
    std::stringstream logStream;


    // Build and initialise two lower triangular matrixes:
    //     (*) < psi_i | H | psi_j > for the variational approach,
    //     (*) < psi_i | H - H_0 | psi_j > for the perturbative approach.
    size_t nbStates = states.size();

    psi_i_H_psi_j.resize(nbStates, std::vector<double>());
    psi_i_HminusH0_psi_j.resize(nbStates, std::vector<double>());
    
    for (size_t i = 0; i < nbStates; ++i)
    {
        psi_i_H_psi_j[i].resize(i + 1, 0.0);
        psi_i_HminusH0_psi_j[i].resize(i + 1, 0.0);
    }

    // Determine the number of charges
    size_t nbCharges = chargesNucleiContributions.size();
    if (ionicMatrixes.size() != nbCharges)
    {
        std::string errorMessage = "Error in ComputeEnergyWithPointCharges::computeHamiltonianMatrix(): the first dimension of ionicMatrixes does not match the dimension of chargesNucleiContributions.";
        print_error(errorMessage);

        std::exit(1);
    }

    
    // Compute matrix elements < psi_i | H | psi_j > and < psi_i | H - H_0 | psi_j >
    size_t i, j;
    if (verbose >= 1)
    {
        logStream << "Matrix elements < psi_i | H | psi_j > and < psi_i | H - H_0 | psi_j >:" << std::endl;
        log(logStream, outputStream);
    }
    else
    {
        #ifdef ENABLE_OPENMP
        #pragma omp parallel for private(i, j)
        #endif
    }
    for (i = 0; i < nbStates; ++i)
    {
        const ExcitedState& psi_i = states[i];

        for (j = 0; j <= i; ++j)
        {
            const ExcitedState& psi_j = states[j];

            // Initialize < psi_i | H | psi_j > matrix element
            double matrixElement = 0.0;

            // Compute < psi_i | H_0 | psi_j >
            double h0Contribution = (i == j ? states[i].get_energy() : 0.0);
            matrixElement += h0Contribution;
            if (verbose >= 2)
            {
                logStream << "< " << i << " | H_0 | " << j << " > = " << std::setprecision(12) << h0Contribution << std::endl;
                log(logStream, outputStream);
            }

            // Compute < psi_i | V_ions/nuclei | psi_j >
            double nuclearContribution = 0.0;
            for (size_t chargeIndex = 0; chargeIndex < nbCharges; ++chargeIndex)
            {
                nuclearContribution += (i == j ? chargesNucleiContributions[chargeIndex] : 0.0);

                if (verbose >= 3)
                {
                    logStream << "< " << i << " | V" << chargeIndex + 1 << "_nuclei | " << j << " > = " << std::setprecision(12) << nuclearContribution << std::endl;
                    log(logStream, outputStream);
                }
            }
            matrixElement += nuclearContribution;
            if (verbose >= 2)
            {
                logStream << "< " << i << " | V_ions/nuclei | " << j << " > = " << std::setprecision(12) << nuclearContribution << std::endl;
                log(logStream, outputStream);
            }

            // Compute < psi_i | V_ion/electrons | psi_j >
            double chargeContribution = 0.0;
            for (size_t chargeIndex = 0; chargeIndex < nbCharges; ++chargeIndex)
            {
                double currentChargeContribution = ExcitedState::ionicPotential(psi_i, psi_j, ionicMatrixes[chargeIndex]);
                chargeContribution += currentChargeContribution;
                if (verbose >= 3)
                {
                    logStream << "< " << i << " | V" << chargeIndex + 1 << "_electrons | " << j << " > = " << std::setprecision(12) << currentChargeContribution << std::endl;
                    log(logStream, outputStream);
                }
            }
            matrixElement += chargeContribution;
            if (verbose >= 2)
            {
                logStream << "< " << i << " | V_ions/electrons | " << j << " > = " << chargeContribution << std::endl;
                log(logStream, outputStream);
            }

            // Store < psi_i | H | psi_j > matrix element
            psi_i_H_psi_j[i][j] = matrixElement;
            psi_i_HminusH0_psi_j[i][j] = psi_i_H_psi_j[i][j] - h0Contribution;

            if (verbose >= 1)
            {
                logStream << "< " << i << " | H | " << j << " > = " << std::setprecision(12) << psi_i_H_psi_j[i][j] << std::endl;
                logStream << "< " << i << " | H - H0| " << j << " > = " << std::setprecision(12) << psi_i_HminusH0_psi_j[i][j] << std::endl;
                log(logStream, outputStream);
            }
            if (verbose >= 2 && j != i)
            {
                logStream << std::endl;
                log(logStream, outputStream);
            }
        }

        if (verbose >= 1)
        {
            logStream << std::endl;
            log(logStream, outputStream);
        }
    }
}

void ComputeEnergyWithPointCharges::computeResults(const std::vector<ExcitedState>& states, const std::vector<std::vector<double>>& psi_i_H_psi_j, const std::vector<std::vector<double>>& psi_i_HminusH0_psi_j, const std::string& outputFilePrefix, std::ostream& outputStream, int verbose)
{
    std::stringstream logStream;
    size_t nbStates = states.size();


    // Diagonalize < psi_i | H | psi_j > matrix
    std::vector<double> eigenvalues;
    std::vector<std::vector<double>> eigenvectors;
    findEigenValuesAndEigenVectorsOfSymmetricalMatrix(psi_i_H_psi_j, eigenvalues, eigenvectors);

    if (verbose >= 3)
    {
        logStream << "Eigenvalues:" << std::endl;
        log(logStream, outputStream);
        logStream << std::scientific;
        logStream << std::setprecision(10);
        for (size_t k = 0; k < eigenvalues.size(); ++k)
        {
            logStream << eigenvalues[k] << ' ';
        }
        logStream << std::endl << std::endl;
        log(logStream, outputStream);

        logStream << "Eigenvectors (columns): " << std::endl;
        log(logStream, outputStream);
        logStream << std::scientific;
        logStream << std::setprecision(10);
        for (size_t i = 0; i < eigenvectors.size(); ++i)
        {
            for (size_t j = 0; j < eigenvectors[i].size(); ++j)
            {
                logStream << std::right << std::setw(17) << eigenvectors[i][j] << '\t';
            }

            logStream << std::endl;
        }
        logStream << std::defaultfloat << std::endl;
        log(logStream, outputStream);
    }
    

    // Sort eigenvalues and eigenvectors
    sortEigenValuesAndEigenVectors(eigenvalues, eigenvectors);

    std::ofstream outputFile(outputFilePrefix + "_energies.cdftt");
    if (!outputFile)
    {
        std::stringstream errorMessage;
        errorMessage << "Error in Job::computeResultsEnergyWithPointCharges(): could not open output file " << outputFilePrefix << "_energies.cdftt for writing." << std::endl;
        
        print_error(errorMessage.str());
        
        std::exit(1);
    }

    logStream << "Sorted Eigenvalues:" << std::endl;
    log(logStream, outputStream);
    logStream << std::scientific;
    logStream << std::setprecision(10);
    outputFile << std::scientific;
    outputFile << std::setprecision(10);
    for (size_t k = 0; k < eigenvalues.size(); ++k)
    {
        logStream << eigenvalues[k] << ' ';
        outputFile << eigenvalues[k] << std::endl;
    }
    logStream << std::endl << std::endl;
    log(logStream, outputStream);

    outputFile.close();
    outputFile.open(outputFilePrefix + "_eigenvectors.cdftt");
    if (!outputFile)
    {
        std::stringstream errorMessage;
        errorMessage << "Error in Job::computeResultsEnergyWithPointCharges(): could not open output file " << outputFilePrefix << "_eigenvectors.cdftt for writing." << std::endl;
        
        print_error(errorMessage.str());
        
        std::exit(1);
    }

    logStream << "Sorted Eigenvectors (columns): " << std::endl;
    log(logStream, outputStream);
    logStream << std::scientific;
    logStream << std::setprecision(10);
    outputFile << std::scientific;
    outputFile << std::setprecision(10);
    for (size_t i = 0; i < eigenvectors.size(); ++i)
    {
        for (size_t j = 0; j < eigenvectors[i].size(); ++j)
        {
            logStream << std::right << std::setw(17) << eigenvectors[i][j] << '\t';
            outputFile << std::right << std::setw(17) << eigenvectors[i][j] << ' ';
        }

        logStream << std::endl;
        outputFile << std::endl;
    }
    logStream << std::defaultfloat << std::endl;
    log(logStream, outputStream);

    outputFile.close();


    if (verbose >= 3)
    {
        // Compute projections of perturbed states onto unperturbed basis and in terms of Slater determinants
        logStream << "Projection onto unperturbed basis:" << std::endl;
        log(logStream, outputStream);

        for (size_t i = 0; i < nbStates; ++i)
        {
            std::vector<std::pair<double, int>> contributions;
            contributions.reserve(nbStates);

            std::vector<std::pair<double, SlaterDeterminant>> contributions_SD;

            logStream << "Perturbed state " << i << " (E = " << std::setprecision(10) << eigenvalues[i] << " H):" << std::endl;
            logStream << "  | " << i << "' > = ";
            log(logStream, outputStream);
            
            bool firstTerm = true;
            for (size_t k = 0; k < nbStates; ++k)
            {
                double c_k = eigenvectors[k][i];
                double c_k_squared = c_k * c_k;

                if (c_k != 0.0)
                {
                    contributions.emplace_back(c_k_squared, k);

                    if (!firstTerm && c_k > 0)
                    {
                        logStream << " + ";
                    }
                    else if (c_k < 0)
                    {
                        logStream << " - ";
                    }

                    logStream << std::setprecision(6) << std::abs(c_k) << " | " << k << " >";
                    firstTerm = false;

                    for (const std::pair<SlaterDeterminant, double>& slaterCoef : states[k].get_slaterDeterminants())
                    {
                        // Search for the Slater determinant in the contributions_SD vector
                        auto it = std::find_if(contributions_SD.begin(),
                                               contributions_SD.end(),
                                               [&slaterCoef](const std::pair<double, SlaterDeterminant>& element)
                                               { return element.second == slaterCoef.first; });

                        // If it is not found, add it to the contributions with its contribution.
                        if (it == contributions_SD.end())
                        {
                            contributions_SD.emplace_back(c_k * slaterCoef.second, slaterCoef.first);
                        }
                        else
                        {
                            it->first += c_k * slaterCoef.second;
                        }
                    }
                }
            }
            logStream << std::endl;
            log(logStream, outputStream);
            
            // Show dominant contributions
            std::sort(contributions.begin(), contributions.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
            logStream << "  Main contributions:" << std::endl;
            log(logStream, outputStream);
            for (size_t ii = 0; ii < std::min(size_t(5), contributions.size()); ++ii)
            {
                if (contributions[ii].first > 1e-6)
                {
                    size_t k = contributions[ii].second;
                    logStream << "    State " << k << ": "
                              << std::setprecision(6) << std::setw(10) << contributions[ii].first * 100 << " %"
                              << "  (c_" << k << " = " << std::setprecision(8) << eigenvectors[k][i] << ")" << std::endl;
                }
            }
            logStream << std::endl;
            log(logStream, outputStream);

            logStream << "Expansion in terms of Slater determinants:" << std::endl;
            logStream << "  | " << i << "' > = ";
            log(logStream, outputStream);

            firstTerm = true;
            for (size_t l = 0; l < contributions_SD.size(); ++l)
            {
                double c_l = contributions_SD[l].first;
                //double c_l_squared = c_l * c_l;

                if (!firstTerm && c_l >= 0)
                {
                    logStream << " + ";
                }
                else if (c_l < 0)
                {
                    logStream << " - ";
                }

                logStream << std::setprecision(6) << std::abs(c_l) << " | D_" << l << " >";
                firstTerm = false;
            }
            logStream << std::endl << "  where:" << std::endl;
            for (size_t l = 0; l < contributions_SD.size(); ++l)
            {
                logStream << "    | D_" << l << " > = " << contributions_SD[l].second << std::endl;
            }
            
            logStream << std::endl;
            log(logStream, outputStream);
        }
    }



    logStream << "------ Perturbative approach (cf. Guégan et al., PCCP 2020) ------" << std::endl << std::endl;
    log(logStream, outputStream);

    bool warningPrinted = false;

    // Compute dp_k for each state using Eq. (27) in Guégan et al., PCCP 2020
    std::vector<double> dpk_perturb_state0_withoutRenormalisation(nbStates, 0.0);
    std::vector<double> normalisationFactors(nbStates, 0.0);
    std::vector<std::vector<double>> dpk_perturb(nbStates, std::vector<double>(nbStates, 0.0));

    // Compute extra-diagonal dp_k coefficients
    for (size_t i = 0; i < nbStates; ++i)
    {
        for (size_t j = 0; j < nbStates; ++j)
        {
            if (i != j)
            {
                double Ei_minus_Ej = states[i].get_energy() - states[j].get_energy();

                // Check for degeneracy to avoid division by zero
                if (std::abs(Ei_minus_Ej) >= 1e-10)
                {
                    // psi_i_HminusH0_psi_j is a lower triangular matrix
                    dpk_perturb[i][j] = (j <= i ? psi_i_HminusH0_psi_j[i][j] : psi_i_HminusH0_psi_j[j][i]) / Ei_minus_Ej;
                    dpk_perturb[i][j] *= dpk_perturb[i][j];

                    normalisationFactors[i] += dpk_perturb[i][j];

                    // For the first state, we keep the unrenormalised dp_k coefficients so we can compare them with the paper and the variational approach later.
                    if (i == 0)
                    {
                        dpk_perturb_state0_withoutRenormalisation[j] = dpk_perturb[i][j];
                    }

                    if (dpk_perturb[i][j] > 1.0)
                    {
                        dpk_perturb[i][j] = 0.0;
                        if (i < j) // to avoid printing twice the same warning for the pair (i, j) and (j, i)
                        {
                            warningPrinted = true;

                            logStream << "Warning: the dp_" << j << " coefficient for the state " << i << " and the dp_" << i << " coefficient for the state " << j << " are greater than 1 (dp_" << j << " = " << dpk_perturb[i][j] << ")." << std::endl;
                            logStream << "They will be set to 0 to maintain the normalisation condition on dp_k (limitation of the perturbative approach)." << std::endl;
                            log(logStream, outputStream);
                        }
                    }
                }
                else
                {
                    dpk_perturb[i][j] = 0.0;

                    if (i < j) // to avoid printing twice the same warning for the pair (i, j) and (j, i)
                    {
                        warningPrinted = true;

                        logStream << "Warning: degeneracy detected between states " << i << " and " << j << " (|E_i - E_j| < 1e-10)." << std::endl;
                        logStream << "The dp_" << j << "coefficient for the state " << i << "and the dp_" << i << " coefficient for the state " << j << " will be set to zero to avoid division by zero." << std::endl;
                        log(logStream, outputStream);
                    }
                }
            }
        }
    }

    // Compute dp_0 without renormalisation
    double sum = 0.0;
    for (size_t i = 1; i < nbStates; ++i)
    {
        sum += dpk_perturb_state0_withoutRenormalisation[i];
    }
    dpk_perturb_state0_withoutRenormalisation[0] = 1.0 -  sum;

    // Renormalization of dp_k coefficients to ensure that their sum is equal to 1 for each state (normalisation condition)
    for (size_t i = 0; i < nbStates; ++i)
    {
        if (normalisationFactors[i] > 1.0)
        {
            for (size_t j = 0; j < nbStates; ++j)
            {
                if (i != j)
                {
                    dpk_perturb[i][j] = dpk_perturb[i][j] / (1.0 + normalisationFactors[i]);
                }
            }
        }
    }

    // Compute diagonal dp_k coefficients using the normalisation condition
    for (size_t i = 0; i < nbStates; ++i)
    {
        double sumExtraDiagonal = 0.0;
        for (size_t j = 0; j < nbStates; ++j)
        {
            sumExtraDiagonal += (i != j ? dpk_perturb[i][j] : 0.0);
        }

        dpk_perturb[i][i] = 1.0 / (1.0 + sumExtraDiagonal);
    }

    if (warningPrinted)
    {
        logStream << std::endl;
        log(logStream, outputStream);
    }

    logStream << "dp_k values for ground state, using Eq. (27). Excited states are on the columns:" << std::endl;
    log(logStream, outputStream);
    logStream << std::scientific;
    logStream << std::setprecision(10);
    for (size_t i = 0; i < nbStates; ++i)
    {
        logStream << std::right << std::setw(17) << dpk_perturb_state0_withoutRenormalisation[i] << '\t';
    }
    logStream << std::defaultfloat << std::endl << std::endl;
    log(logStream, outputStream);

    logStream << "Renormalized dp_k values. Excited states are on the columns:" << std::endl;
    log(logStream, outputStream);
    logStream << std::scientific;
    logStream << std::setprecision(10);
    for (size_t i = 0; i < nbStates; ++i)
    {
        for (size_t j = 0; j < nbStates; ++j)
        {
            logStream << std::right << std::setw(17) << dpk_perturb[i][j] << '\t';
        }

        logStream << std::endl;
    }
    logStream << std::defaultfloat << std::endl;
    log(logStream, outputStream);

    // Compute E_polarisation and dS for each state using respectively Eq. (26) and Eq. (32) in Guégan et al., PCCP 2020
    std::vector<double> dS_perturb(nbStates, 0.0);
    std::vector<double> E_pola_perturb(nbStates, 0.0);
    double dS_perturb_state0_withoutRenormalisation = 0.0;
    double E_pola_perturb_state0_withoutRenormalisation = 0.0;

    logStream << "E_polarisation and dS using respectively Eq. (26) and Eq. (32) in Guégan et al., PCCP 2020:" << std::endl;
    log(logStream, outputStream);

    for (size_t i = 0; i < nbStates; ++i)
    {
        double sum_dS = 0.0;
        double sum_Epola = 0.0;
        for (size_t j = 0; j < nbStates; ++j)
        {
            if (dpk_perturb[i][j] != 0)
            {
                sum_dS -= dpk_perturb[i][j] * std::log(dpk_perturb[i][j]);

                if (i != j)
                {
                    // Note : degeneracy is already handled in the computation of dp_k coefficients
                    // So we can safely compute the energy difference here without checking for division by zero again.
                    sum_Epola -= dpk_perturb[i][j] * (states[j].get_energy() - states[i].get_energy());
                }
            }
        }

        dS_perturb[i] = sum_dS * Constants::BOLTZMANN_CONSTANT * Constants::AVOGADRO_CONSTANT;
        E_pola_perturb[i] = sum_Epola * Constants::HARTREE_TO_JOULE * Constants::AVOGADRO_CONSTANT;

        // Treating the ground state without renormalisation separately
        if (dpk_perturb_state0_withoutRenormalisation[i] != 0)
        {
            dS_perturb_state0_withoutRenormalisation -= dpk_perturb_state0_withoutRenormalisation[i] * std::log(dpk_perturb_state0_withoutRenormalisation[i]);

            if (i != 0)
            {
                E_pola_perturb_state0_withoutRenormalisation -= dpk_perturb_state0_withoutRenormalisation[i] * (states[0].get_energy() - states[i].get_energy());
            }
        }
    }

    dS_perturb_state0_withoutRenormalisation *= Constants::BOLTZMANN_CONSTANT * Constants::AVOGADRO_CONSTANT;
    E_pola_perturb_state0_withoutRenormalisation *= Constants::HARTREE_TO_JOULE * Constants::AVOGADRO_CONSTANT;

    logStream << "dS (J/mol/K) and |E_polarisation| (J/mol) for ground state without renormalisation:" << std::endl;
    log(logStream, outputStream);
    logStream << std::scientific;
    logStream << std::setprecision(10);
    logStream << std::right << std::setw(17) << dS_perturb_state0_withoutRenormalisation << '\t';
    logStream << std::right << std::setw(17) << std::abs(E_pola_perturb_state0_withoutRenormalisation) << std::endl << std::endl;

    logStream << "dS (J/mol/K) for each state:" << std::endl;
    for (size_t i = 0; i < nbStates; ++i)
    {
        logStream << std::right << std::setw(17) << dS_perturb[i] << '\t';
    }
    logStream << std::endl << std::endl;

    logStream << "|E_polarisation| (J/mol) for each state:" << std::endl;
    for (size_t i = 0; i < nbStates; ++i)
    {
        logStream << std::right << std::setw(17) << std::abs(E_pola_perturb[i]) << '\t';
    }
    logStream << std::defaultfloat << std::endl;
    log(logStream, outputStream);

    
    logStream << std::endl << "------ Variational approach ------" << std::endl << std::endl;

    // Compute dp_k
    std::vector<std::vector<double>> dpk_varia(nbStates, std::vector<double>(nbStates, 0.0));

    logStream << "dp_k values:" << std::endl;
    log(logStream, outputStream);
    logStream << std::scientific;
    logStream << std::setprecision(10);
    for (size_t i = 0; i < eigenvectors.size(); ++i)
    {
        for (size_t j = 0; j < eigenvectors[i].size(); ++j)
        {
            dpk_varia[i][j] = eigenvectors[i][j] * eigenvectors[i][j];
            logStream << std::right << std::setw(17) << dpk_varia[i][j] << '\t';
        }

        logStream << std::endl;
    }
    logStream << std::defaultfloat << std::endl;
    log(logStream, outputStream);

    // Search for the excited state that has the maximum contribution from the ground state to compare with the perturbative approach.
    size_t maxGroundContributionExcitedState = 0;
    double maxContribution = dpk_varia[0][0];
    for (size_t i = 1; i < nbStates; ++i)
    {
        if (dpk_varia[0][i] > maxContribution)
        {
            maxContribution = dpk_varia[0][i];
            maxGroundContributionExcitedState = i;
        }
    }

    // Compute E_polarisation and dS for each state
    std::vector<double> dS_varia(nbStates, 0.0);
    std::vector<double> E_pola_varia(nbStates, 0.0);

    for (size_t i = 0; i < nbStates; ++i)
    {
        double sum_dS = 0.0;
        double sum_Epola = 0.0;
        for (size_t j = 0; j < nbStates; ++j)
        {
            if (dpk_varia[j][i] != 0)
            {
                sum_dS -= dpk_varia[j][i] * std::log(dpk_varia[j][i]);

                if (i != j)
                {
                    // Note : degeneracy is already handled in the computation of dp_k coefficients
                    // So we can safely compute the energy difference here without checking for division by zero again.
                    sum_Epola -= dpk_varia[j][i] * (states[i].get_energy() - states[j].get_energy());
                }
            }
        }

        dS_varia[i] = sum_dS * Constants::BOLTZMANN_CONSTANT * Constants::AVOGADRO_CONSTANT;
        E_pola_varia[i] = sum_Epola * Constants::HARTREE_TO_JOULE * Constants::AVOGADRO_CONSTANT;
    }

    logStream << "dS (J/mol/K) and |E_polarisation| (J/mol) for the excited state with maximum contribution from the ground state | 0 > (| " << maxGroundContributionExcitedState << "' >, with dp_" << maxGroundContributionExcitedState << " = " << std::setprecision(10) << dpk_varia[0][maxGroundContributionExcitedState] << ")." << std::endl;
    logStream << std::scientific;
    logStream << std::setprecision(10);
    logStream << std::right << std::setw(17) << dS_varia[maxGroundContributionExcitedState] << '\t';
    logStream << std::right << std::setw(17) << std::abs(E_pola_varia[maxGroundContributionExcitedState]) << std::endl << std::endl;
    logStream<< "dS (J/mol/K) for each state:" << std::endl;
    for (size_t i = 0; i < nbStates; ++i)
    {
        logStream << std::right << std::setw(17) << dS_varia[i] << '\t';
    }
    logStream << std::endl << std::endl;

    logStream << "|E_polarisation| (J/mol) for each state:" << std::endl;
    for (size_t i = 0; i < nbStates; ++i)
    {
        logStream << std::right << std::setw(17) << std::abs(E_pola_varia[i]) << '\t';
    }
    logStream << std::defaultfloat << std::endl;
    log(logStream, outputStream);
}

void ComputeEnergyWithPointCharges::computeResultsLinearResponse(const std::vector<std::vector<double>>& eigenvalues, const std::vector<std::vector<std::vector<double>>>& ionicPotentialVectors, std::ostream& outputStream, int verbose)
{
    std::stringstream logStream;

    double energy_pseudoOrbitals = 0.0;
    for (size_t i = 0; i < ionicPotentialVectors.size(); ++i)
    {
        for (int spin = 0; spin < 2; ++spin)
        {
            for (size_t j = 0; j < eigenvalues[spin].size(); ++j)
            {
                energy_pseudoOrbitals += eigenvalues[spin][j] * ionicPotentialVectors[i][spin][j] * ionicPotentialVectors[i][spin][j];
            }
        }
    }
    energy_pseudoOrbitals *= 0.5;

    logStream << std::scientific;
    logStream << std::setprecision(10);
    logStream << "|E_polarisation| = " << std::abs(energy_pseudoOrbitals) * Constants::HARTREE_TO_JOULE * Constants::AVOGADRO_CONSTANT << " J/mol." << std::endl;
    log(logStream, outputStream);
}

void ComputeEnergyWithPointCharges::printResults(const std::vector<ExcitedState>& states, const std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& ionicPotentialMatrixes, const std::vector<std::vector<double>>& chargeNucleiContributions, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms, const std::vector<Atom>& atoms, const std::string& outputPrefix, std::ostream& outputStream, int verbose)
{
    std::stringstream logStream;

    size_t nbCharges = charges.size();
    size_t nbChargePositions = chargesPositions.size();


    // Compute Hamiltonian matrixes
    std::vector<std::vector<std::vector<double>>> psi_i_H_psi_j_matrixes;
    std::vector<std::vector<std::vector<double>>> psi_i_HminusH0_psi_j_matrixes;

    if (loopOnAtoms)
    {
        const int maxRunNumber = nbCharges * nbChargePositions;

        for (size_t chargeIndex = 0; chargeIndex < nbCharges; ++chargeIndex)
        {
            for (size_t atomIndex = 0; atomIndex < atoms.size(); ++atomIndex)
            {
                int runNumber = chargeIndex * nbChargePositions + atomIndex + 1;
                std::string runNumberStr = int_to_string_withLeadingZeros(runNumber, maxRunNumber);

                logStream << "====================== RUN #" << runNumberStr
                          << ": charge of " << charges[chargeIndex]
                          << " e on " << atoms[atomIndex].get_name()
                          << " at position (" << std::setprecision(10) << chargesPositions[atomIndex][0] << ", " << chargesPositions[atomIndex][1] << ", " << chargesPositions[atomIndex][2]
                          << ") ======================"
                          << std::defaultfloat << std::endl
                          << std::endl;
                log(logStream, outputStream);

                std::vector<std::vector<std::vector<std::vector<double>>>> currentIonicMatrix(1, ionicPotentialMatrixes[chargeIndex][atomIndex]);
                std::vector<double> currentChargeNucleiContribution = std::vector<double>(1, chargeNucleiContributions[chargeIndex][atomIndex]);

                std::vector<std::vector<double>> psi_i_H_psi_j;
                std::vector<std::vector<double>> psi_i_HminusH0_psi_j;

                computeHamiltonianMatrix(states, currentChargeNucleiContribution, currentIonicMatrix, psi_i_H_psi_j, psi_i_HminusH0_psi_j, outputStream, verbose);
                computeResults(states, psi_i_H_psi_j, psi_i_HminusH0_psi_j, outputPrefix + "_run" + runNumberStr, outputStream, verbose);

                logStream << std::endl;
                log(logStream, outputStream);
            }
        }
    }
    else
    {
        std::vector<std::vector<std::vector<std::vector<double>>>> currentIonicMatrixes(nbCharges);
        std::vector<double> currentChargeNucleiContribution(nbCharges);
        for (size_t i = 0; i < nbCharges; ++i)
        {
            currentIonicMatrixes[i] = ionicPotentialMatrixes[i][0];
            currentChargeNucleiContribution[i] = chargeNucleiContributions[i][0];
        }

        std::vector<std::vector<double>> psi_i_H_psi_j;
        std::vector<std::vector<double>> psi_i_HminusH0_psi_j;

        computeHamiltonianMatrix(states, currentChargeNucleiContribution, currentIonicMatrixes, psi_i_H_psi_j, psi_i_HminusH0_psi_j, outputStream, verbose);
        computeResults(states, psi_i_H_psi_j, psi_i_HminusH0_psi_j, outputPrefix, outputStream, verbose);
    }
}

void ComputeEnergyWithPointCharges::printResultsLinearResponse(const std::vector<std::vector<double>>& eigenvalues, const std::vector<std::vector<std::vector<std::vector<double>>>>& ionicPotentialVectors, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms, const std::vector<Atom>& atoms, std::ostream& outputStream, int verbose)
{
    size_t nbCharges = charges.size();
    size_t nbChargePositions = chargesPositions.size();

    std::stringstream logStream;

    // Compute polarization energy
    if (loopOnAtoms)
    {
        const int maxRunNumber = nbCharges * nbChargePositions;

        for (size_t chargeIndex = 0; chargeIndex < nbCharges; ++chargeIndex)
        {
            for (size_t atomIndex = 0; atomIndex < atoms.size(); ++atomIndex)
            {
                int runNumber = chargeIndex * nbChargePositions + atomIndex + 1;
                std::string runNumberStr = int_to_string_withLeadingZeros(runNumber, maxRunNumber);

                logStream << "====================== RUN #" << runNumberStr
                          << ": charge of " << charges[chargeIndex]
                          << " e on " << atoms[atomIndex].get_name()
                          << " at position (" << std::setprecision(10) << chargesPositions[atomIndex][0] << ", " << chargesPositions[atomIndex][1] << ", " << chargesPositions[atomIndex][2]
                          << ") ======================"
                          << std::defaultfloat << std::endl
                          << std::endl;
                log(logStream, outputStream);

                std::vector<std::vector<std::vector<double>>> currentIonicPotentialVectors(1, ionicPotentialVectors[chargeIndex][atomIndex]);

                computeResultsLinearResponse(eigenvalues, currentIonicPotentialVectors, outputStream, verbose);

                logStream << std::defaultfloat << std::endl;
                log(logStream, outputStream);
            }
        }
    }
    else
    {
        std::vector<std::vector<std::vector<double>>> currentIonicPotentialVectors(nbCharges);

        for (size_t i = 0; i < nbCharges; ++i)
        {
            currentIonicPotentialVectors[i] = ionicPotentialVectors[i][0];
        }
        
        computeResultsLinearResponse(eigenvalues, currentIonicPotentialVectors, outputStream, verbose);

        logStream << std::defaultfloat << std::endl;
        log(logStream, outputStream);
    }
}


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//----------------------------------------------------------------------------------------------------//

ComputeEnergyWithPointCharges::ComputeEnergyWithPointCharges(const std::string& inputFileName):
    Job(inputFileName)
{ }


//----------------------------------------------------------------------------------------------------//
// STATIC METHODS
//----------------------------------------------------------------------------------------------------//

void ComputeEnergyWithPointCharges::computeChargeNucleiContributions(const std::vector<Atom>& atoms, std::vector<std::vector<double>>& chargeNucleiContributions, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms, double nuclearCutoff)
{
    size_t nbCharges = charges.size();
    size_t nbChargePositions = chargesPositions.size();

    chargeNucleiContributions.resize(nbCharges);

    for (size_t i = 0; i < nbCharges; ++i)
    {
        if (loopOnAtoms)
        {
            chargeNucleiContributions[i].resize(nbChargePositions, 0.0);

            for (size_t j = 0; j < nbChargePositions; ++j)
            {
                for (const Atom& atom : atoms)
                {
                    double distance = std::sqrt((atom.get_coordinates()[0] - chargesPositions[j][0]) * (atom.get_coordinates()[0] - chargesPositions[j][0])
                                                + (atom.get_coordinates()[1] - chargesPositions[j][1]) * (atom.get_coordinates()[1] - chargesPositions[j][1])
                                                + (atom.get_coordinates()[2] - chargesPositions[j][2]) * (atom.get_coordinates()[2] - chargesPositions[j][2]));

                    if (distance > nuclearCutoff)
                    {
                        chargeNucleiContributions[i][j] += charges[i] * atom.get_atomicNumber() / distance;
                    }
                }
            }
        }
        else
        {
            chargeNucleiContributions[i].resize(1, 0.0);

            for (const Atom& atom : atoms)
            {
                double distance = std::sqrt((atom.get_coordinates()[0] - chargesPositions[i][0]) * (atom.get_coordinates()[0] - chargesPositions[i][0])
                                            + (atom.get_coordinates()[1] - chargesPositions[i][1]) * (atom.get_coordinates()[1] - chargesPositions[i][1])
                                            + (atom.get_coordinates()[2] - chargesPositions[i][2]) * (atom.get_coordinates()[2] - chargesPositions[i][2]));

                if (distance > nuclearCutoff)
                {
                    chargeNucleiContributions[i][0] += charges[i] * atom.get_atomicNumber() / distance;
                }
            }
        }
    }
}

void ComputeEnergyWithPointCharges::computeIonicPotentialMatrixesFromBecke(Becke& becke, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& ionicPotentialMatrixes, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms, int kmax, int lebedev_order, int radial_grid_factor)
{
    size_t nbCharges = charges.size();
    size_t nbChargePositions = chargesPositions.size();

    ionicPotentialMatrixes.resize(charges.size());

    if (loopOnAtoms)
    {
        for (size_t i = 0; i < nbCharges; ++i)
        {
            // In the looping case, each charge has multiple positions (one for each atom).
            // So we need to compute the ionic matrixes for each position of the charge.
            ionicPotentialMatrixes[i].resize(nbChargePositions);

            for (size_t j = 0; j < nbChargePositions; ++j)
            {
                ionicPotentialMatrixes[i][j] = becke.getIonicPotentialMatrix(chargesPositions[j], charges[i], kmax, lebedev_order, radial_grid_factor);
            }
        }
    }
    else
    {
        for (size_t i = 0; i < nbCharges; ++i)
        {
            // In the non-looping case, each charge has only one position.
            // So we compute the ionic matrixes only once for each charge and store it in the first position of the second dimension of the vector.
            ionicPotentialMatrixes[i].resize(1);
            ionicPotentialMatrixes[i][0] = becke.getIonicPotentialMatrix(chargesPositions[i], charges[i], kmax, lebedev_order, radial_grid_factor);
        }
    }
}

void ComputeEnergyWithPointCharges::computeIonicPotentialMatrixesFromGrid(const Orbitals& orbitals, Grid& grid, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& ionicPotentialMatrixes, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms)
{
    size_t nbCharges = charges.size();
    size_t nbChargePositions = chargesPositions.size();

    ionicPotentialMatrixes.resize(charges.size());

    if (loopOnAtoms)
    {
        for (size_t i = 0; i < nbCharges; ++i)
        {
            // In the looping case, each charge has multiple positions (one for each atom).
            // So we need to compute the ionic matrixes for each position of the charge.
            ionicPotentialMatrixes[i].resize(nbChargePositions);

            for (size_t j = 0; j < nbChargePositions; ++j)
            {
                ionicPotentialMatrixes[i][j] = grid.getIonicPotentialMatrix(orbitals, chargesPositions[j], charges[i]);
            }
        }
    }
    else
    {
        for (size_t i = 0; i < nbCharges; ++i)
        {
            // In the non-looping case, each charge has only one position.
            // So we compute the ionic matrixes only once for each charge and store it in the first position of the second dimension of the vector.
            ionicPotentialMatrixes[i].resize(1);
            ionicPotentialMatrixes[i][0] = grid.getIonicPotentialMatrix(orbitals, chargesPositions[i], charges[i]);
        }
    }
}

void ComputeEnergyWithPointCharges::computeIonicPotentialMatrixesFromOrbitals(Orbitals& orbitals, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& ionicPotentialMatrixes, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms)
{
    size_t nbCharges = charges.size();
    size_t nbChargePositions = chargesPositions.size();

    ionicPotentialMatrixes.resize(charges.size());

    if (loopOnAtoms)
    {
        for (size_t i = 0; i < nbCharges; ++i)
        {
            // In the looping case, each charge has multiple positions (one for each atom).
            // So we need to compute the ionic matrixes for each position of the charge.
            ionicPotentialMatrixes[i].resize(nbChargePositions);

            for (size_t j = 0; j < nbChargePositions; ++j)
            {
                ionicPotentialMatrixes[i][j] = orbitals.getIonicPotentialMatrix(chargesPositions[j], charges[i]);
            }
        }
    }
    else
    {
        for (size_t i = 0; i < nbCharges; ++i)
        {
            // In the non-looping case, each charge has only one position.
            // So we compute the ionic matrixes only once for each charge and store it in the first position of the second dimension of the vector.
            ionicPotentialMatrixes[i].resize(1);
            ionicPotentialMatrixes[i][0] = orbitals.getIonicPotentialMatrix(chargesPositions[i], charges[i]);
        }
    }
}

void ComputeEnergyWithPointCharges::computeIonicPotentialVectorsFromBecke(Becke& becke, std::vector<std::vector<std::vector<std::vector<double>>>>& ionicPotentialVectors, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms, int kmax, int lebedev_order, int radial_grid_factor)
{
    size_t nbCharges = charges.size();
    size_t nbChargePositions = chargesPositions.size();

    ionicPotentialVectors.resize(charges.size());

    if (loopOnAtoms)
    {
        for (size_t i = 0; i < nbCharges; ++i)
        {
            // In the looping case, each charge has multiple positions (one for each atom).
            // So we need to compute the ionic matrixes for each position of the charge.
            ionicPotentialVectors[i].resize(nbChargePositions);

            for (size_t j = 0; j < nbChargePositions; ++j)
            {
                ionicPotentialVectors[i][j] = becke.getIonicPotentialVector(chargesPositions[j], charges[i], kmax, lebedev_order, radial_grid_factor);
            }
        }
    }
    else
    {
        for (size_t i = 0; i < nbCharges; ++i)
        {
            // In the non-looping case, each charge has only one position.
            // So we compute the ionic matrixes only once for each charge and store it in the first position of the second dimension of the vector.
            ionicPotentialVectors[i].resize(1);
            ionicPotentialVectors[i][0] = becke.getIonicPotentialVector(chargesPositions[i], charges[i], kmax, lebedev_order, radial_grid_factor);
        }
    }
}

void ComputeEnergyWithPointCharges::computeIonicPotentialVectorsFromOrbitals(Orbitals& orbitals, std::vector<std::vector<std::vector<std::vector<double>>>>& ionicPotentialVectors, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms)
{
    size_t nbCharges = charges.size();
    size_t nbChargePositions = chargesPositions.size();

    ionicPotentialVectors.resize(charges.size());

    if (loopOnAtoms)
    {
        for (size_t i = 0; i < nbCharges; ++i)
        {
            // In the looping case, each charge has multiple positions (one for each atom).
            // So we need to compute the ionic matrixes for each position of the charge.
            ionicPotentialVectors[i].resize(nbChargePositions);

            for (size_t j = 0; j < nbChargePositions; ++j)
            {
                ionicPotentialVectors[i][j] = orbitals.getIonicPotentialVector_unitPseudoCgtf(chargesPositions[j], charges[i]);
            }
        }
    }
    else
    {
        for (size_t i = 0; i < nbCharges; ++i)
        {
            // In the non-looping case, each charge has only one position.
            // So we compute the ionic matrixes only once for each charge and store it in the first position of the second dimension of the vector.
            ionicPotentialVectors[i].resize(1);
            ionicPotentialVectors[i][0] = orbitals.getIonicPotentialVector_unitPseudoCgtf(chargesPositions[i], charges[i]);
        }
    }
}

void ComputeEnergyWithPointCharges::computeLinearResponseFunctionMatrix(const Orbitals& orbitals, const std::vector<std::vector<std::vector<std::vector<double>>>>& tripleOrbitalIntegralMatrix, std::vector<std::vector<std::vector<double>>>& lrfMatrix)
{
    // Get number of MOs
    int numberOfMo = orbitals.get_numberOfMo();

    // Get occupied and virtual orbitals numbers
    std::vector<std::vector<int>> occupiedOrbitalsNumbers;
    std::vector<std::vector<int>> virtualOrbitalsNumbers;
    orbitals.getOccupiedAndVirtualOrbitalNumbers(occupiedOrbitalsNumbers, virtualOrbitalsNumbers);

    // Get orbital energies
    const std::vector<std::vector<double>>& orbitalEnergies = orbitals.get_orbitalEnergy();


    // Build and initialise the lower triangular LRF matrix for each spin
    lrfMatrix.resize(2, std::vector<std::vector<double>>(numberOfMo, std::vector<double>()));
    for (int spin = 0; spin < 2; ++spin)
    {
        for (int i = 0; i < numberOfMo; ++i)
        {
            lrfMatrix[spin][i].resize(i + 1, 0.0);
        }
    }

    std::cout << std::scientific;
    std::cout << std::setprecision(10);
    for (int spin = 0; spin < 2; ++spin)
    {
        std::cout << "Computing LRF matrix for " << (spin == static_cast<int>(SpinType::ALPHA) ? "Alpha" : "Beta") << " spin (analytical):" << std::endl;

        for (int i = 0; i < numberOfMo; ++i)
        {
            for (int j = 0; j <= i; ++j)
            {
                for (int occupiedOrbital : occupiedOrbitalsNumbers[spin])
                {
                    int occupiedOrbitalIndex = occupiedOrbital - 1; // because occupiedOrbitalsNumbers are 1-based

                    for (int virtualOrbital : virtualOrbitalsNumbers[spin])
                    {
                        int virtualOrbitalIndex = virtualOrbital - 1; // because virtualOrbitalsNumbers are 1-based

                        std::array<int, 3> indices_i = {i, occupiedOrbitalIndex, virtualOrbitalIndex};
                        std::array<int, 3> indices_j = {j, occupiedOrbitalIndex, virtualOrbitalIndex};

                        std::sort(indices_i.begin(), indices_i.end(), std::greater<size_t>());
                        std::sort(indices_j.begin(), indices_j.end(), std::greater<size_t>());

                        lrfMatrix[spin][i][j] += tripleOrbitalIntegralMatrix[spin][indices_i[0]][indices_i[1]][indices_i[2]]
                                                  * tripleOrbitalIntegralMatrix[spin][indices_j[0]][indices_j[1]][indices_j[2]]
                                                  / (orbitalEnergies[spin][occupiedOrbitalIndex] - orbitalEnergies[spin][virtualOrbitalIndex]);
                    }
                }

                lrfMatrix[spin][i][j] *= 2.0;
                std::cout << "< phi_" << i + 1 << " | Xi | phi_" << j + 1<< " > = " << lrfMatrix[spin][i][j] << std::endl;
            }

            std::cout << std::endl;
        }

        std::cout << std::endl;
    }
}


//----------------------------------------------------------------------------------------------------//
// PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

void ComputeEnergyWithPointCharges::run()
{
    // Read output file prefix
    std::string outputPrefix;
    readOutputPrefix(outputPrefix);


    // Read progress bar display option
    bool showProgress;
    readShowProgress(showProgress);

    
    // Read verbose level and open log file if needed
    int verbose;
    readVerbose(verbose);

    std::stringstream logStream;

    std::ofstream logFile;
    if (verbose != 0)
    {
        logFile.open(outputPrefix + "_log.cdftt");
        if (!logFile)
        {
            std::cout << "Warning: could not open log file " << outputPrefix << "_log.cdftt for writing." << std::endl;
            std::cout << "The program will still display logging information on standard output." << std::endl << std::endl;
        }
    }
    std::ostream& outputStream = ((verbose != 0 && logFile) ? logFile : std::cout);


    // Read analytic file name
    std::vector<std::string> analyticFilesNames;
    readAnalyticFilesNames(analyticFilesNames);


    // Check number of analytic files names
    if (analyticFilesNames.size() != 1)
    {
        if (analyticFilesNames[0].substr(analyticFilesNames[0].length() - 4) == ".log"
            || analyticFilesNames[0].substr(analyticFilesNames[0].length() - 5) == ".fchk")
        {
            std::stringstream errorMessage;
            errorMessage << "Error: incorrect number of analytic files names (one file expected)." << std::endl;
            errorMessage << "Please check the documentation and the number of files specified in the \"AnalyticFiles\" parameter in " << _inputFileName << '.';

            print_error(errorMessage.str(), outputStream);

            std::exit(1);
        }
        else if (analyticFilesNames.size() != 2)
        {
            std::stringstream errorMessage;
            errorMessage << "Error: incorrect number of analytic files names (one or two file(s) expected)." << std::endl;
            errorMessage << "Please check the documentation and the number of files specified in the \"AnalyticFiles\" parameter in " << _inputFileName << '.';

            print_error(errorMessage.str(), outputStream);

            std::exit(1);
        }
    }


    // Read cutoff distance for nuclear contribution
    double nuclearCutoff;
    readNuclearCutoff(nuclearCutoff);
    
    
    // Read transitions file
    std::string transitionsFileName;
    readTransitionsFileName(transitionsFileName);


    // Read max number of transitions to read in transitions file
    int maxNbExcitedStates;
    readMaxNumberOfExcitedStates(maxNbExcitedStates);


    // Loading orbitals
    Orbitals orbitals;
    computeOrbitalsOrBecke<Orbitals>(orbitals, analyticFilesNames[0]);
    std::cout << std::endl;

    double groundStateEnergy = orbitals.get_energy();
    std::string groundStateEnergySource = "analytic file " + analyticFilesNames[0];

    // The "GroundStateEnergy" parameter in the input file has priority over the other ways of reading the ground state energy.
    // This allows to use a different ground state energy than the one in the orbitals file if needed.
    if (readGroundStateEnergy(groundStateEnergy))
    {
        groundStateEnergySource = "input file " + _inputFileName;
        orbitals.set_energy(groundStateEnergy);
    }
    else if (orbitals.get_energy() == 0.0)
    {
        if (analyticFilesNames.size() == 2)
        {
            ExcitedState::readGroundStateEnergy(analyticFilesNames[1], groundStateEnergy);
            groundStateEnergySource = "analytic file " + analyticFilesNames[1];
        }
        else
        {
            if (ExcitedState::readGroundStateEnergy(transitionsFileName, groundStateEnergy))
            {
                groundStateEnergySource = "transitions file " + transitionsFileName;
            }
            else
            {
                std::stringstream errorMessage;
                errorMessage << "Error: cannot determine the ground state energy." << std::endl;
                errorMessage << "Please check the documentation and the \"GroundStateEnergy\" parameter value or the \"AnalyticalFiles\" parameter value in the provided input file (" << _inputFileName << ")." << std::endl << std::endl;

                print_error(errorMessage.str(), outputStream);

                std::exit(1);
            }
        }

        orbitals.set_energy(groundStateEnergy);
    }

    logStream << "Ground State Energy read from " << groundStateEnergySource << ": " << std::setprecision(10) << groundStateEnergy << " H." << std::endl;
    log(logStream, outputStream);

    // Keep a const reference on orbitals' atoms
    const std::vector<Atom>& atoms = orbitals.get_struct().get_atoms();


    // Read point charges
    std::vector<double> charges;
    readCharges(charges);
    size_t nbCharges = charges.size();


    // Read point charges positions
    bool loopOnAtoms = false;
    std::vector<std::array<double, 3>> chargesPositions;
    readPositions(chargesPositions);

    if (chargesPositions.empty())
    {
        logStream << "Note: the \"Positions\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        logStream << "The program will place the point charge" << (nbCharges > 1 ? "s" : "") << " on each atom successively." << std::endl << std::endl;
        log(logStream, outputStream);

        loopOnAtoms = true;
        for (const Atom& atom : atoms)
        {
            chargesPositions.push_back(atom.get_coordinates());
        }
    }
    size_t nbChargePositions = chargesPositions.size();


    // Check number of charges positions
    if (!loopOnAtoms && nbChargePositions != nbCharges)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect number of point charges positions." << std::endl;
        errorMessage << "Please check the documentation and the positions specified in the \"ChargesPositions\" parameter in " << _inputFileName << '.';

        print_error(errorMessage.str(), outputStream);

        std::exit(1);
    }


    // Print charges information
    logStream << "Number of point charges: " << nbCharges << std::endl;
    log(logStream, outputStream);
    if (!loopOnAtoms)
    {
        for (size_t i = 0; i < nbCharges; ++i)
        {
            logStream << "Point charge #" << i + 1 << ": " << charges[i] << " e at position (" << std::setprecision(10) << chargesPositions[i][0] << ", " << chargesPositions[i][1] << ", " << chargesPositions[i][2] << ")." << std::defaultfloat << std::endl;
        }
    }
    else
    {
        for (size_t i = 0; i < nbCharges; ++i)
        {
            for (size_t j = 0; j < nbChargePositions; ++j)
            {
                logStream << "Run #" << i * nbChargePositions + j + 1 << ": point charge #" << i + 1 << " of " << charges[i] << " e, on " << atoms[j].get_name() << " at position (" << std::setprecision(10) << chargesPositions[j][0] << ", " << chargesPositions[j][1] << ", " << chargesPositions[j][2] << ")." << std::defaultfloat << std::endl;
            }
        }
    }
    logStream << std::endl;
    log(logStream, outputStream);


    if (verbose >= 3)
    {
        logStream << "Molecular orbitals:" << std::endl;
        logStream << orbitals << std::endl;
        log(logStream, outputStream);
    }
    

    // Get Ground Slater Determinant
    SlaterDeterminant groundStateSlaterDeterminant(orbitals);
    ExcitedState groundState(orbitals.get_energy(), groundStateSlaterDeterminant);


    // Building states vector
    std::vector<ExcitedState> states;
    states.push_back(groundState);
    

    // Reading transitions file
    if (!transitionsFileName.empty())
    {
        std::cout << "Reading transitions from file: " << transitionsFileName << ". Please wait..." << std::endl;
        ExcitedState::readTransitions(transitionsFileName, states, groundState.get_energy(), maxNbExcitedStates);
    }
    else
    {
        std::cout << "Reading transitions from analytic file: " << analyticFilesNames[0] << ". Please wait..." << std::endl;
        ExcitedState::readTransitions(analyticFilesNames[0], states, groundState.get_energy(), maxNbExcitedStates);
    }

    size_t nbStates = states.size();
    logStream << "Total number of states: " << nbStates << std::endl << std::endl;
    log(logStream, outputStream);


    // Compute Slater Determinants from electronic transitions for each state
    for (ExcitedState& state : states)
    {
        state.computeSlaterDeterminants(groundStateSlaterDeterminant);

        if (verbose >= 1)
        {
            logStream << state << std::endl;
            log(logStream, outputStream);
        }
    }


    // Compute ions-nuclei interactions only once
    std::vector<std::vector<double>> chargeNucleiContributions;
    computeChargeNucleiContributions(atoms, chargeNucleiContributions, charges, chargesPositions, loopOnAtoms, nuclearCutoff);
    

    /************/
    /* ANALYTIC */
    /************/
    
    logStream << std::endl << std::endl
              << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl
              << "|||||||||||||||||||||||||||||||||||||||||||||||                          |||||||||||||||||||||||||||||||||||||||||||||||" << std::endl
              << "|||||||||||||||||||||||||||||||||||||||||||||||   ANALYTIC COMPUTATION   |||||||||||||||||||||||||||||||||||||||||||||||" << std::endl
              << "|||||||||||||||||||||||||||||||||||||||||||||||                          |||||||||||||||||||||||||||||||||||||||||||||||" << std::endl
              << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl << std::endl;
    log(logStream, outputStream);


    // Compute ionic matrixes for each charge and position
    // We build a 5D vector of dimensions [charge][position][spin][MO][MO] to store the ionic matrixes for each charge and position.
    // The first dimension corresponds to the point charge index.
    // The second dimension corresponds to the charge position index (in case the program has to loop over atom positions).
    // The third dimension corresponds to the spin (0 for alpha, 1 for beta).
    // The fourth and fifth dimensions correspond to the MO indexes (matrix elements i and j, with j <= i: lower triangular matrix).
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> ionicMatrixes;
    computeIonicPotentialMatrixesFromOrbitals(orbitals, ionicMatrixes, charges, chargesPositions, loopOnAtoms);


    // Print results
    printResults(states, ionicMatrixes, chargeNucleiContributions, charges, chargesPositions, loopOnAtoms, atoms, outputPrefix, outputStream, verbose);
    

    /****************/
    /* REGULAR GRID */
    /****************/
    
    // Read grid size
    GridSize gridSize;
    CustomSizeData customSizeData;
    if (readSize(gridSize, customSizeData))
    {
        logStream << std::endl << std::endl
                  << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl
                  << "|||||||||||||||||||||||||||||||||||||||||||||                              |||||||||||||||||||||||||||||||||||||||||||||" << std::endl
                  << "|||||||||||||||||||||||||||||||||||||||||||||   REGULAR GRID COMPUTATION   |||||||||||||||||||||||||||||||||||||||||||||" << std::endl
                  << "|||||||||||||||||||||||||||||||||||||||||||||                              |||||||||||||||||||||||||||||||||||||||||||||" << std::endl
                  << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl << std::endl;
        log(logStream, outputStream);


        // Setting orbitals
        std::vector<SpinType> orbitalsSpins;
        std::vector<int> orbitalsNumbers;
        setOrbitals(orbitals, orbitals.get_numberOfMo(), orbitalsNumbers, orbitalsSpins, OrbitalType::ALL, SpinType::ALPHA_BETA);


        // Building domain and grid
        std::cout << "Building domain and grid, please wait..." << std::endl;
        Domain domain = buildDomainForCube(orbitals, gridSize, customSizeData, orbitalsNumbers.size());
        Grid orbitalsGrid = orbitals.makeOrbGrid(domain, orbitalsNumbers, orbitalsSpins, showProgress);
        std::cout << std::endl;


        // Compute ionic matrixes for each charge and position
        std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> ionicMatrixes_regularGrid;
        computeIonicPotentialMatrixesFromGrid(orbitals, orbitalsGrid, ionicMatrixes_regularGrid, charges, chargesPositions, loopOnAtoms);


        // Print results
        printResults(states, ionicMatrixes_regularGrid, chargeNucleiContributions, charges, chargesPositions, loopOnAtoms, atoms, outputPrefix + "_regularGrid", outputStream, verbose);
    }
    

    /****************/
    /* BECKE GRID */
    /****************/
    
    // Read Becke grid parameters
    std::vector<int> beckeParams;
    if (readBecke(beckeParams))
    {
        if (beckeParams.size() == 0)
        {
            beckeParams = { 3, 41, 5 };
        }

        Becke becke;
        computeOrbitalsOrBecke<Becke>(becke, analyticFilesNames[0]);

        
        logStream << std::endl << std::endl
                  << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl
                  << "||||||||||||||||||||||||||||||||||||||||||||||                             |||||||||||||||||||||||||||||||||||||||||||||" << std::endl
                  << "||||||||||||||||||||||||||||||||||||||||||||||   BECKE GRID COMPUTATION    |||||||||||||||||||||||||||||||||||||||||||||" << std::endl
                  << "||||||||||||||||||||||||||||||||||||||||||||||                             |||||||||||||||||||||||||||||||||||||||||||||" << std::endl
                  << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl << std::endl;
        log(logStream, outputStream);


        // Compute ionic matrixes for each charge and position
        std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> ionicMatrixes_becke;
        computeIonicPotentialMatrixesFromBecke(becke, ionicMatrixes_becke, charges, chargesPositions, loopOnAtoms);


        // Print results
        printResults(states, ionicMatrixes_becke, chargeNucleiContributions, charges, chargesPositions, loopOnAtoms, atoms, outputPrefix + "_becke", outputStream, verbose);
    }
    

    /*
    ////////////////////////////////////
    // Debug - V_nuclear calculations //
    ////////////////////////////////////
    logStream << std::endl << std::endl << std::endl;
    logStream << "==============================================================================================" << std::endl;
    logStream << "====================== DEBUG - Computation of < ϕ_0 | V_nuclear | ϕ_0 > ======================" << std::endl;
    logStream << "==============================================================================================" << std::endl << std::endl;
    
    logStream << "====================== ANALYTIC COMPUTATION ======================" << std::endl << std::endl;
    log(logStream, outputStream);

    // Compute Nuclear Matrices for each atom (print AO matrix elements)
    logStream << "AO / MO basis:" <<std::endl;
    logStream << "--------------" << std::endl << std::endl;
    log(logStream, outputStream);

    std::vector<std::vector<std::vector<std::vector<double>>>> nuclearMatrices(orbitals.get_numberOfAtoms());
    int atomIndex = 0;
    int nbAtoms = static_cast<int>(orbitals.get_struct().get_atoms().size());
    for (const Atom& atom : orbitals.get_struct().get_atoms())
    {
        logStream << "Computing nuclear matrix for atom " << atom.get_name() << ": Z = " << atom.get_atomicNumber() << " ; position = (" << atom.get_coordinates()[0] << ", " << atom.get_coordinates()[1] << ", " << atom.get_coordinates()[2] << ")..." << std::endl;

        nuclearMatrices[atomIndex] = orbitals.getIonicPotentialMatrix(atom.get_coordinates(), atom.get_atomicNumber(), true, atomIndex == nbAtoms - 1, atomIndex == nbAtoms - 1);
        ++atomIndex;
    }
    logStream << std::endl;
    log(logStream, outputStream);

    double sum_phi_i_Vnuclear_phi_j = 0.0;
    for (int atomIndex = 0; atomIndex < orbitals.get_numberOfAtoms(); ++atomIndex)
    {
        for (size_t spin = 0; spin < nuclearMatrices[atomIndex].size(); ++spin)
        {
            for (size_t i = 0; i < nuclearMatrices[atomIndex][spin].size(); ++i)
            {
                for (size_t j = 0; j <= i; ++j)
                {
                    sum_phi_i_Vnuclear_phi_j += (i == j ? nuclearMatrices[atomIndex][spin][i][j] : 2.0 * nuclearMatrices[atomIndex][spin][i][j]);
                }
            }
        }
    }
    logStream << "Total sum of MO matrix elements for Alpha and Beta spins (analytic): " << std::setprecision(10) << sum_phi_i_Vnuclear_phi_j << std::endl << std::endl;
    log(logStream, outputStream);
    
    double sum_psi_i_Vnuclear_psi_j = 0.0;
    double V_ij = 0.0;
    
    // Print Slater determinant matrix elements
    logStream << "-------------------------------------------------------------------------------------------" << std::endl << std::endl;
    logStream << "States basis:" << std::endl;
    logStream << "-------------" << std::endl << std::endl;
    log(logStream, outputStream);

    for (size_t i = 0; i < nbStates; ++i)
    {
        logStream << "ψ_" << i << ": energy = " << std::setprecision(10) << states[i].get_energy() << " Hartree" << std::endl;
        for (const auto& slaterCoeff : states[i].get_slaterDeterminants())
        {
            logStream << "    " << slaterCoeff.first << "; Coefficient: " << slaterCoeff.second << std::endl;
        }
        logStream << std::endl;
    }
    log(logStream, outputStream);

    for (size_t i = 0; i < nbStates; ++i)
    {
        // Get Slater Determinants for state i
        std::vector<std::pair<SlaterDeterminant, double>> slaterDeterminants_i(states[i].get_slaterDeterminants());

        for (size_t j = 0; j <= i; ++j)
        {
            logStream << "Computing < ψ_" << i << " | V_nuclear | ψ_" << j << " > matrix element..." << std::endl;

            // Initialize < D_k | V_nuclear | D_l > matrix element
            V_ij = 0.0;

            // Get Slater Determinants for state j
            std::vector<std::pair<SlaterDeterminant, double>> slaterDeterminants_j(states[j].get_slaterDeterminants());

            logStream << "    Computing < D_k | V_nuclear | D_l > matrix elements:" << std::endl;
            // Compute < D_k | V_nuclear | D_l >
            for (const std::pair<SlaterDeterminant, double>& slaterCoeff_k : slaterDeterminants_i)
            {
                for (const std::pair<SlaterDeterminant, double>& slaterCoeff_l : slaterDeterminants_j)
                {
                    double D_kl = 0.0;

                    for (int atomIndex = 0; atomIndex < orbitals.get_numberOfAtoms(); ++atomIndex)
                    {
                        D_kl += SlaterDeterminant::ionicPotential(slaterCoeff_k.first, slaterCoeff_l.first, nuclearMatrices[atomIndex]) * (slaterCoeff_k.second * slaterCoeff_l.second);
                    }

                    V_ij += D_kl;
                    logStream << std::right << std::setw(17) << D_kl << '\t';
                }
                logStream << std::endl;
            }

            sum_psi_i_Vnuclear_psi_j += (i == j ? V_ij : 2.0 * V_ij);
            logStream << std::endl;
        }
    }
    log(logStream, outputStream);
    

    // COMPARAISON AVEC CALCUL SUR GRILLE
    GridSize gridSize;
    CustomSizeData customSizeData;
    if (readSize(gridSize, customSizeData))
    {
        logStream << std::endl << std::endl << "====================== REGULAR GRID COMPUTATION ======================" << std::endl << std::endl;
        log(logStream, outputStream);

        // Setting orbitals
        std::vector<SpinType> orbitalsSpins;
        std::vector<int> orbitalsNumbers;
        setOrbitals(orbitals, orbitals.get_numberOfMo(), orbitalsNumbers, orbitalsSpins, OrbitalType::ALL, SpinType::ALPHA_BETA);

        // Building domain and grid
        std::cout << "Building domain and grid, please wait..." << std::endl;
        Domain domain = buildDomainForCube(orbitals, gridSize, customSizeData, orbitalsNumbers.size());
        Grid orbitalsGrid = orbitals.makeOrbGrid(domain, orbitalsNumbers, orbitalsSpins, showProgress);
        std::cout << std::endl;

        double sum_phi_i_Vnuclear_phi_j_alpha = 0.0;
        double sum_phi_i_Vnuclear_phi_j_beta = 0.0;
        double V_ij = 0.0;

        logStream << "Computing ionic potential matrix in MO basis for Alpha spin..." << std::endl;
        log(logStream, outputStream);

        int nbStepsTotal = orbitals.get_numberOfMo() * (orbitals.get_numberOfMo() + 1) / 2;
        int currentStep = 0;
        int lastProgress = -1;
        if (showProgress)
        {
            print_progressBar(0, nbStepsTotal, lastProgress);
        }

        std::cout << std::scientific;
        std::cout << std::setprecision(10);
        logStream << std::scientific;
        logStream << std::setprecision(10);
        for (int i = 0; i < orbitals.get_numberOfMo(); ++i)
        {
            for (int j = 0; j <= i; ++j)
            {
                V_ij = 0.0;
                for (const Atom& atom : orbitals.get_struct().get_atoms())
                {
                    V_ij += orbitalsGrid.phiStarVionicStarPhi(i, j, SpinType::ALPHA, atom.get_coordinates(), atom.get_atomicNumber());
                }

                sum_phi_i_Vnuclear_phi_j_alpha += (i == j ? V_ij : 2.0 * V_ij);
                logStream << std::right << std::setw(17) << V_ij << '\t';

                currentStep++;
            }
            logStream << std::endl;

            if (showProgress)
            {
                print_progressBar(currentStep, nbStepsTotal, lastProgress);
            }
        }
        if (showProgress)
        {
            std::cout << std::endl;
        }
        logStream << std::defaultfloat << "Total sum of MO matrix elements for Alpha spin: " << std::setprecision(10) << sum_phi_i_Vnuclear_phi_j_alpha << std::endl;
        log(logStream, outputStream);

        logStream << "Ionic potential matrix in MO basis for Beta spin:" << std::endl;
        log(logStream, outputStream);

        currentStep = 0;
        lastProgress = -1;
        if (showProgress)
        {
            print_progressBar(0, nbStepsTotal, lastProgress);
        }

        std::cout << std::scientific;
        std::cout << std::setprecision(10);
        logStream << std::scientific;
        logStream << std::setprecision(10);

        for (int i = 0; i < orbitals.get_numberOfMo(); ++i)
        {
            for (int j = 0; j <= i; ++j)
            {
                V_ij = 0.0;
                for (const Atom& atom : orbitals.get_struct().get_atoms())
                {
                    V_ij += orbitalsGrid.phiStarVionicStarPhi(i, j, SpinType::BETA, atom.get_coordinates(), atom.get_atomicNumber());
                }

                sum_phi_i_Vnuclear_phi_j_beta += (i == j ? V_ij : 2.0 * V_ij);
                logStream << std::right << std::setw(17) << V_ij << '\t';

                currentStep++;
            }
            logStream << std::endl;

            if (showProgress)
            {
                print_progressBar(currentStep, nbStepsTotal, lastProgress);
            }
        }
        if (showProgress)
        {
            std::cout << std::endl;
        }
        logStream << std::defaultfloat << "Total sum of MO matrix elements for Beta spin: " << std::setprecision(10) << sum_phi_i_Vnuclear_phi_j_beta << std::endl << std::endl;
        log(logStream, outputStream);

        double sum_phi_i_Vnuclear_phi_j = sum_phi_i_Vnuclear_phi_j_alpha + sum_phi_i_Vnuclear_phi_j_beta;
        logStream << "Total sum of MO matrix elements for Alpha and Beta spins (regular Grid): " << std::setprecision(10) << sum_phi_i_Vnuclear_phi_j << std::endl;
        log(logStream, outputStream);
    }

    
    // COMPARAISON AVEC BECKE
    logStream << std::endl << std::endl << "====================== BECKE GRID COMPUTATION ======================" << std::endl << std::endl;
    log(logStream, outputStream);

    // Compute Nuclear Matrices for each atom (print AO matrix elements)
    logStream << "AO / MO basis:" << std::endl;
    logStream << "--------------" << std::endl << std::endl;
    log(logStream, outputStream);

    std::vector<std::vector<std::vector<std::vector<double>>>> nuclearMatricesBecke(orbitals.get_numberOfAtoms());
    atomIndex = 0;
    for (const Atom& atom : orbitals.get_struct().get_atoms())
    {
        logStream << "Computing nuclear matrix for atom " << atom.get_name() << ": Z = " << atom.get_atomicNumber() << " ; position = (" << atom.get_coordinates()[0] << ", " << atom.get_coordinates()[1] << ", " << atom.get_coordinates()[2] << ")..." << std::endl;
        log(logStream, outputStream);

        nuclearMatricesBecke[atomIndex] = becke.getIonicPotentialMatrix(atom.get_coordinates(), atom.get_atomicNumber(), beckeParams[0], beckeParams[1], beckeParams[2]);
        ++atomIndex;
    }
    logStream << std::endl;
    log(logStream, outputStream);

    sum_phi_i_Vnuclear_phi_j = 0.0;
    for (atomIndex = 0; atomIndex < orbitals.get_numberOfAtoms(); ++atomIndex)
    {
        for (size_t spin = 0; spin < nuclearMatricesBecke[atomIndex].size(); ++spin)
        {
            for (size_t i = 0; i < nuclearMatricesBecke[atomIndex][spin].size(); ++i)
            {
                for (size_t j = 0; j <= i; ++j)
                {
                    sum_phi_i_Vnuclear_phi_j += (i == j ? nuclearMatricesBecke[atomIndex][spin][i][j] : 2.0 * nuclearMatricesBecke[atomIndex][spin][i][j]);
                }
            }
        }
    }
    logStream << "Total sum of MO matrix elements for Alpha and Beta spins (Becke): " << std::setprecision(10) << sum_phi_i_Vnuclear_phi_j << std::endl << std::endl;
    log(logStream, outputStream);
    */

    if (verbose != 0 && logFile)
    {
        logFile.close();
    }
}


void ComputeEnergyWithPointCharges::run_computeLinearResponseWithPointCharges()
{
    // Read output file prefix
    std::string outputPrefix;
    readOutputPrefix(outputPrefix);


    // Read option to save pseudo orbitals in cube format
    bool savePseudoOrbitals;
    readSavePseudoOrbitals(savePseudoOrbitals);


    // Read progress bar display option
    bool showProgress;
    readShowProgress(showProgress);


    // Read verbose level and open log file if needed
    int verbose;
    readVerbose(verbose);

    std::stringstream logStream;
    
    std::ofstream logFile;
    if (verbose != 0)
    {
        logFile.open(outputPrefix + "_log.cdftt");
        if (!logFile)
        {
            std::cout << "Warning: could not open log file " << outputPrefix << "_log.cdftt for writing." << std::endl;
            std::cout << "The program will still display logging information on standard output." << std::endl << std::endl;
        }
    }
    std::ostream& outputStream = ((verbose != 0 && logFile) ? logFile : std::cout);

    
    //Read analytic file name
    std::vector<std::string> analyticFilesNames;
    readAnalyticFilesNames(analyticFilesNames);


    // TODO: check number of analytic files


    // Loading orbitals
    std::cout << "Building Orbitals object... ";
    Orbitals orbitals;
    computeOrbitalsOrBecke<Orbitals>(orbitals, analyticFilesNames[0]);
    std::cout << std::endl;

    // Keep a const reference on orbitals' atoms
    const std::vector<Atom>& atoms = orbitals.get_struct().get_atoms();

    
    // Read point charges
    std::vector<double> charges;
    readCharges(charges);
    size_t nbCharges = charges.size();


    // Read point charges positions
    bool loopOnAtoms = false;
    std::vector<std::array<double, 3>> chargesPositions;
    readPositions(chargesPositions);

    if (chargesPositions.empty())
    {
        logStream << "Note: the \"Positions\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        logStream << "The program will place the point charge" << (nbCharges > 1 ? "s" : "") << " on each atom successively." << std::endl << std::endl;
        log(logStream, outputStream);

        loopOnAtoms = true;
        for (const Atom& atom : atoms)
        {
            chargesPositions.push_back(atom.get_coordinates());
        }
    }
    size_t nbChargePositions = chargesPositions.size();


    // Check number of charges positions
    if (!loopOnAtoms && nbChargePositions != nbCharges)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect number of point charges positions." << std::endl;
        errorMessage << "Please check the documentation and the positions specified in the \"ChargesPositions\" parameter in " << _inputFileName << '.';

        print_error(errorMessage.str(), outputStream);

        std::exit(1);
    }


    // Print charges information
    logStream << "Number of point charges: " << nbCharges << std::endl;
    log(logStream, outputStream);
    if (!loopOnAtoms)
    {
        for (size_t i = 0; i < nbCharges; ++i)
        {
            logStream << "Point charge #" << i + 1 << ": " << charges[i] << " e at position (" << std::setprecision(10) << chargesPositions[i][0] << ", " << chargesPositions[i][1] << ", " << chargesPositions[i][2] << ")." << std::defaultfloat << std::endl;
        }
    }
    else
    {
        for (size_t i = 0; i < nbCharges; ++i)
        {
            for (size_t j = 0; j < nbChargePositions; ++j)
            {
                logStream << "Run #" << i * nbChargePositions + j + 1 << ": point charge #" << i + 1 << " of " << charges[i] << " e, on " << atoms[j].get_name() << " at position (" << std::setprecision(10) << chargesPositions[j][0] << ", " << chargesPositions[j][1] << ", " << chargesPositions[j][2] << ")." << std::defaultfloat << std::endl;
            }
        }
    }
    logStream << std::endl;
    log(logStream, outputStream);
    

    /************/
    /* ANALYTIC */
    /************/

    logStream << std::endl << std::endl
              << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl
              << "|||||||||||||||||||||||||||||||||||||||||||||||                          |||||||||||||||||||||||||||||||||||||||||||||||" << std::endl
              << "|||||||||||||||||||||||||||||||||||||||||||||||   ANALYTIC COMPUTATION   |||||||||||||||||||||||||||||||||||||||||||||||" << std::endl
              << "|||||||||||||||||||||||||||||||||||||||||||||||                          |||||||||||||||||||||||||||||||||||||||||||||||" << std::endl
              << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl << std::endl;
    log(logStream, outputStream);

    
    // Compute triple-orbital-integral matrix only once.
    // We build a 4D vector of dimensions [spin][MO_i][MO_j][MO_k] to store the triple-orbital-integral matrix.
    // The first dimension corresponds to the spin (0 for alpha, 1 for beta).
    // The second, third and fourth dimensions correspond to the indices of the matrix elements (i, j and k) with k <= j <= i (lower triangular matrixes).
    std::vector<std::vector<std::vector<std::vector<double>>>> tripleOrbitalIntegralMatrix = orbitals.getTripleOrbitalIntegralMatrix(showProgress);


    // Compute linear response function (LRF) matrix
    std::vector<std::vector<std::vector<double>>> lrfMatrix;
    computeLinearResponseFunctionMatrix(orbitals, tripleOrbitalIntegralMatrix, lrfMatrix);


    // Diagonalize LRF matrix to get pseudo orbitals from eigenvectors
    std::vector<std::vector<double>> eigenvalues;
    std::vector<std::vector<std::vector<double>>> eigenvectors;
    Orbitals pseudoOrbitals = computePseudoOrbitalsFromLrfMatrix(orbitals, lrfMatrix, eigenvalues, eigenvectors, outputPrefix, savePseudoOrbitals, outputStream, verbose, showProgress);


    // Compute ionic vectors obtained with a pseudo CGTF made from a unit pseudo GTF (exponent = 0, coefficient = 1) only once.
    // We build a 4D of dimensions [charge][position][spin][MO] to store the ionic potential vectors for each charge and position.
    // The first dimension corresponds to the charge index.
    // The second dimension corresponds to the charge position index (in case the program has to loop over atom positions).
    // The third dimension corresponds to the spin (0 for alpha, 1 for beta).
    // The fourth dimension corresponds to the MO index.
    std::vector<std::vector<std::vector<std::vector<double>>>> ionicPotentialVectors;
    computeIonicPotentialVectorsFromOrbitals(pseudoOrbitals, ionicPotentialVectors, charges, chargesPositions, loopOnAtoms);


    // Print results
    printResultsLinearResponse(eigenvalues, ionicPotentialVectors, charges, chargesPositions, loopOnAtoms, atoms, outputStream, verbose);
    

    /**************/
    /* BECKE GRID */
    /**************/

    // Read Becke grid parameters
    std::vector<int> beckeParams;
    readBecke(beckeParams);

    if (beckeParams.size() != 0)
    {
        logStream << std::endl << std::endl
                  << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl
                  << "||||||||||||||||||||||||||||||||||||||||||||||                            ||||||||||||||||||||||||||||||||||||||||||||||" << std::endl
                  << "||||||||||||||||||||||||||||||||||||||||||||||   BECKE GRID COMPUTATION   ||||||||||||||||||||||||||||||||||||||||||||||" << std::endl
                  << "||||||||||||||||||||||||||||||||||||||||||||||                            ||||||||||||||||||||||||||||||||||||||||||||||" << std::endl
                  << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl << std::endl;
        log(logStream, outputStream);


        // Build Becke grid
        std::cout << "Building Becke object... ";
        Becke becke;
        computeOrbitalsOrBecke<Becke>(becke, analyticFilesNames[0]);


        // Get triple-orbital-integral matrix for Becke grid
        std::vector<std::vector<std::vector<std::vector<double>>>> tripleOrbitalIntegralMatrix_becke = becke.getTripleOrbitalIntegralMatrix(beckeParams[0], beckeParams[1], beckeParams[2], showProgress);


        // Compute LRF Matrix for Becke grid
        std::vector<std::vector<std::vector<double>>> lrfMatrix_becke;
        computeLinearResponseFunctionMatrix(becke.get_orbitals(), tripleOrbitalIntegralMatrix_becke, lrfMatrix_becke);


        // Diagonalize LRF matrix to get pseudo orbitals from eigenvectors
        std::vector<std::vector<double>> eigenvalues_becke(2);
        std::vector<std::vector<std::vector<double>>> eigenvectors_becke(2);
        Orbitals pseudoOrbitals_becke = computePseudoOrbitalsFromLrfMatrix(becke.get_orbitals(), lrfMatrix_becke, eigenvalues_becke, eigenvectors_becke, outputPrefix + "_becke", savePseudoOrbitals, outputStream, verbose, showProgress);
        
        
        // Compute ionic vectors obtained with a pseudo CGTF made from a unit pseudo GTF (exponent = 0, coefficient = 1) only once.
        std::vector<std::vector<std::vector<std::vector<double>>>> ionicPotentialVectors_becke;
        computeIonicPotentialVectorsFromOrbitals(pseudoOrbitals_becke, ionicPotentialVectors_becke, charges, chargesPositions, loopOnAtoms);


        // Print results
        printResultsLinearResponse(eigenvalues_becke, ionicPotentialVectors_becke, charges, chargesPositions, loopOnAtoms, atoms, outputStream, verbose);



        // DEBUG - Manually compute sigma vectors (i.e. obtaining the same values than the pseudo Orbitals coefficients.)
        /*
        // Compute ionic vectors obtained with a pseudo CGTF made from a unit pseudo GTF (exponent = 0, coefficient = 1) only once.
        std::vector<std::vector<std::vector<std::vector<double>>>> ionicPotentialVectors_becke_debug;
        computeIonicPotentialVectorsFromBecke(becke, ionicPotentialVectors_becke_debug, charges, chargesPositions, loopOnAtoms, beckeParams[0], beckeParams[1], beckeParams[2]);


        // Multiply by eigenvectors values to get sigma vectors (i.e. the pseudo orbitals coefficients in the basis of the original orbitals).
        std::vector<std::vector<std::vector<std::vector<double>>>> sigmaVectors_becke_debug(nbCharges);
        if (loopOnAtoms)
        {
            for (size_t i = 0; i < nbCharges; ++i)
            {
                // In the looping case, each charge has multiple positions (one for each atom).
                // So we need to compute the ionic matrixes for each position of the charge.
                sigmaVectors_becke_debug[i].resize(nbChargePositions, std::vector<std::vector<double>>(2, std::vector<double>(eigenvalues_becke[0].size(), 0.0)));

                for (size_t j = 0; j < nbChargePositions; ++j)
                {
                    for (int spin = 0; spin < 2; ++spin)
                    {
                        for (size_t k = 0; k < eigenvalues_becke[spin].size(); ++k)
                        {
                            for (size_t l = 0; l < eigenvectors_becke[spin].size(); ++l)
                            {
                                sigmaVectors_becke_debug[i][j][spin][k] += eigenvectors_becke[spin][l][k] * ionicPotentialVectors_becke_debug[i][j][spin][l];
                            }
                        }
                    }
                }
            }
        }
        else
        {
            for (size_t i = 0; i < nbCharges; ++i)
            {
                // In the non-looping case, each charge has only one position.
                // So we compute the ionic matrixes only once for each charge and store it in the first position of the second dimension of the vector.
                sigmaVectors_becke_debug[i].resize(1, std::vector<std::vector<double>>(2, std::vector<double>(eigenvalues_becke[0].size(), 0.0)));

                for (int spin = 0; spin < 2; ++spin)
                {
                    for (size_t k = 0; k < eigenvalues_becke[spin].size(); ++k)
                    {
                        for (size_t l = 0; l < eigenvectors_becke[spin].size(); ++l)
                        {
                            sigmaVectors_becke_debug[i][0][spin][k] += eigenvectors_becke[spin][l][k] * ionicPotentialVectors_becke_debug[i][0][spin][l];
                        }
                    }
                }
            }
        }


        // Print results
        printResultsLinearResponseWithPointCharges(eigenvalues_becke, sigmaVectors_becke_debug, charges, chargesPositions, loopOnAtoms, atoms, outputStream, verbose);
        */
    }

    
    if (verbose != 0 && logFile)
    {
        logFile.close();
    }
}
