#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <Becke/Becke.h>
#include <Common/Atom.h>
#include <Common/Descriptors.h>
#include <JobControl/Job.h>
#include <JobControl/Jobs/ComputeCondensedLinearResponse.hpp>
#include <JobControl/Jobs/ComputeDescriptors.hpp>
#include <Orbitals/Orbitals.h>
#include <Utils/Utils.h>


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//----------------------------------------------------------------------------------------------------//

ComputeCondensedLinearResponse::ComputeCondensedLinearResponse(const std::string& inputFileName):
    Job(inputFileName)
{ }


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

void ComputeCondensedLinearResponse::run()
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


    // TODO: check number of analytic files


    // Read Becke grid parameters
    std::vector<int> beckeParams;
    readBecke(beckeParams);

    if (beckeParams.size() == 0)
    {
        beckeParams = { 3, 41, 5 };
    }


    // Build Becke grid
    std::cout << "Building Becke object... ";
    Becke becke;
    computeOrbitalsOrBecke<Becke>(becke, analyticFilesNames[0]);

    // Keep a const reference on atoms
    std::vector<Atom> atoms = becke.get_orbitals().get_struct().get_atoms();


    // Compute condensed linear response function (CLRF) matrix using Becke grid
    std::vector<std::vector<double>> clrfMatrix_becke;
    becke.chiAtomic(clrfMatrix_becke, beckeParams[0], beckeParams[1], beckeParams[2]);


    // Print CLRF matrix
    logStream << "Condensed linear response function (CLRF) matrix:" << std::endl;
    log(logStream, outputStream);
    logStream << std::scientific;
    logStream << std::setprecision(10);
    for (size_t i = 0; i < clrfMatrix_becke.size(); ++i)
    {
        for (size_t j = 0; j < clrfMatrix_becke[i].size(); ++j)
        {
            logStream << std::right << std::setw(17) << clrfMatrix_becke[i][j] << ' ';
        }
        logStream << std::endl;
    }
    logStream << std::defaultfloat << std::endl;
    log(logStream, outputStream);


    // Print χ_{AA} + \sum χ_{AB} for each atom
    logStream << "χ_{AA} + sum(χ_{AB}) for each atom (should be around 0):" << std::endl;
    log(logStream, outputStream);
    logStream << std::scientific;
    logStream << std::setprecision(10);
    for (size_t i = 0; i < clrfMatrix_becke.size(); ++i)
    {
        double sumXi = 0.0;

        for (size_t j = 0; j < clrfMatrix_becke[i].size(); ++j)
        {
            sumXi += (j <= i  ? clrfMatrix_becke[i][j] : clrfMatrix_becke[j][i]); // Because the matrix is triangular
        }

        logStream << "Atom " << i + 1 << " (" << atoms[i].get_symbol() << "): " << sumXi << std::endl;
    }
    logStream << std::defaultfloat << std::endl;
    log(logStream, outputStream);


    // Diagonalize CLRF matrix
    std::vector<double> eigenvalues_becke;
    std::vector<std::vector<double>> eigenvectors_becke;

    // Compute results
    findEigenValuesAndEigenVectorsOfSymmetricalMatrix(clrfMatrix_becke, eigenvalues_becke, eigenvectors_becke);
    sortEigenValuesAndEigenVectors(eigenvalues_becke, eigenvectors_becke);


    // Save results
    std::ofstream outputFile(outputPrefix + "_eigenvalues_becke.cdftt");
    if (!outputFile)
    {
        std::stringstream errorMessage;
        errorMessage << "Error in Job::run_computeCondensedLinearResponse(): could not open output file " << outputPrefix << "_eigenvalues_becke.cdftt for writing." << std::endl;

        print_error(errorMessage.str(), outputStream);
        
        std::exit(1);
    }

    logStream << "Sorted Eigenvalues:" << std::endl;
    log(logStream, outputStream);
    logStream << std::scientific;
    logStream << std::setprecision(10);
    outputFile << std::scientific;
    outputFile << std::setprecision(10);
    for (size_t k = 0; k < eigenvalues_becke.size(); ++k)
    {
        logStream << eigenvalues_becke[k] << ' ';
        outputFile << eigenvalues_becke[k] << std::endl;
    }
    logStream << std::endl << std::endl;
    log(logStream, outputStream);
    outputFile.close();

    outputFile.open(outputPrefix + "_eigenvectors_becke.cdftt");
    if (!outputFile)
    {
        std::stringstream errorMessage;
        errorMessage << "Error in Job::run_computeCondensedLinearResponse(): could not open output file " << outputPrefix << "_eigenvectors_becke.cdftt for writing." << std::endl;

        print_error(errorMessage.str(), outputStream);
        
        std::exit(1);
    }

    logStream << "Sorted Eigenvectors (columns): " << std::endl;
    log(logStream, outputStream);
    logStream << std::scientific;
    logStream << std::setprecision(10);
    outputFile << std::scientific;
    outputFile << std::setprecision(10);
    for (size_t i = 0; i < eigenvectors_becke.size(); ++i)
    {
        for (size_t j = 0; j < eigenvectors_becke[i].size(); ++j)
        {
            logStream << std::right << std::setw(17) << eigenvectors_becke[i][j] << ' ';
            outputFile << std::right << std::setw(17) << eigenvectors_becke[i][j] << ' ';
        }

        logStream << std::endl;
        outputFile << std::endl;
    }
    logStream << std::defaultfloat << std::endl;
    log(logStream, outputStream);
    outputFile.close();


    // Lowest electron density distorsion mode
    for (size_t i = 0; i < eigenvalues_becke.size(); ++i)
    {
        logStream << "Lowest electron density distorsion mode condensed on atom " << i + 1 << " (" << atoms[i].get_symbol() << "): " << std::abs(eigenvalues_becke[0] * eigenvectors_becke[i][0] * eigenvectors_becke[i][0]) << std::endl;
    }
    log(logStream, outputStream);


    /***********************/
    /* ENERGY - FMO method */
    /***********************/

    logStream << std::endl << std::endl << "====================== ENERGY - FMO method ======================" << std::endl << std::endl;
    log(logStream, outputStream);

    Orbitals orbitals = becke.get_orbitals();
    orbitals.init_homoLumoIndexes();
    orbitals.get_f();
    Descriptors descriptors_FMO(orbitals.get_descriptors());
    descriptors_FMO.set_mu_fk_data(orbitals.get_all_f(), orbitals.getHomoEnergy()[static_cast<int>(SpinType::ALPHA)], orbitals.getLumoEnergy()[static_cast<int>(SpinType::ALPHA)]);
    descriptors_FMO.compute_all();

    // Print deltafk
    const std::vector<double>& deltafk_FMO = descriptors_FMO.get_deltafk();

    logStream << "Δf (column)" << std::endl;
    logStream << std::scientific;
    logStream << std::setprecision(10);
    for (size_t i = 0; i < deltafk_FMO.size(); ++i)
    {
        logStream << std::right << std::setw(17) << deltafk_FMO[i] << std::endl;
    }
    logStream << std::endl;
    log(logStream, outputStream);

    // Get dk coefficient
    std::vector<double> dk_coeffs(deltafk_FMO.size(), 0.0);
    for (size_t i = 0; i < deltafk_FMO.size(); ++i)
    {
        for (size_t j = 0; j < deltafk_FMO.size(); ++j)
        {
            dk_coeffs[i] += deltafk_FMO[j] * eigenvectors_becke[j][i];
        }
    }

    logStream << "Δf = ";
    for (size_t i = 0; i < deltafk_FMO.size(); ++i)
    {
        logStream << (dk_coeffs[i] >= 0 ? "+ " : "- ") << std::abs(dk_coeffs[i]) << " χ_" << i + 1 << ' ';
    }
    logStream << std::endl;
    log(logStream, outputStream);

    double energy_FMO = 0.0;
    for (size_t i = 0; i < dk_coeffs.size(); ++i)
    {
        if (std::abs(eigenvalues_becke[i]) > 1e-5)
        {
            energy_FMO += dk_coeffs[i] * dk_coeffs[i] / std::abs(eigenvalues_becke[i]);
        }
    }
    energy_FMO *= (-0.5);

    logStream << "Energy = " << std::scientific << std::setprecision(10) << energy_FMO << " H" << std::endl;
    log(logStream, outputStream);


    /***********************************************/
    /* ENERGY - FD method with Q0, Q+ and Q- files */
    /***********************************************/

    logStream << std::endl << std::endl << "====================== ENERGY - FD method with Q0, Q- and Q+ files ======================" << std::endl << std::endl;
    log(logStream, outputStream);

    Descriptors descriptors_FD = ComputeDescriptors::computeWithFiniteDifferenceBecke(analyticFilesNames[1], analyticFilesNames[2], analyticFilesNames[3], beckeParams[0], beckeParams[1], beckeParams[2]);

    std::cout << descriptors_FD << std::endl;

    // Print deltafk
    const std::vector<double>& deltafk_FD_threeFiles = descriptors_FD.get_deltafk();

    logStream << "Δf (column)" << std::endl;
    logStream << std::scientific;
    logStream << std::setprecision(10);
    for (size_t i = 0; i < deltafk_FD_threeFiles.size(); ++i)
    {
        logStream << std::right << std::setw(17) << deltafk_FD_threeFiles[i] << std::endl;
    }
    logStream << std::endl;
    log(logStream, outputStream);

    // Get dk coefficient
    std::vector<double> dk_coeffs_FD_threeFiles(deltafk_FD_threeFiles.size(), 0.0);
    for (size_t i = 0; i < deltafk_FD_threeFiles.size(); ++i)
    {
        for (size_t j = 0; j < deltafk_FD_threeFiles.size(); ++j)
        {
            dk_coeffs_FD_threeFiles[i] += deltafk_FD_threeFiles[j] * eigenvectors_becke[j][i];
        }
    }

    logStream << "Δf = ";
    for (size_t i = 0; i < deltafk_FD_threeFiles.size(); ++i)
    {
        logStream << (dk_coeffs_FD_threeFiles[i] >= 0 ? "+ " : "- ") << std::abs(dk_coeffs_FD_threeFiles[i]) << " χ_" << i + 1 << ' ';
    }
    logStream << std::endl;
    log(logStream, outputStream);

    double energy_FD_threeFiles = 0.0;
    for (size_t i = 0; i < dk_coeffs_FD_threeFiles.size(); ++i)
    {
        if (std::abs(eigenvalues_becke[i]) > 1e-5)
        {
            energy_FD_threeFiles += dk_coeffs_FD_threeFiles[i] * dk_coeffs_FD_threeFiles[i] / std::abs(eigenvalues_becke[i]);
        }
    }
    energy_FD_threeFiles *= (-0.5);

    logStream << "Energy = " << std::scientific << std::setprecision(10) << energy_FD_threeFiles << " H" << std::endl;
    log(logStream, outputStream);

    /****************************************/
    /* ENERGY - FD method with Q0 file only */
    /****************************************/

    logStream << std::endl << std::endl << "====================== ENERGY - FD method with Q0 file only ======================" << std::endl << std::endl;
    log(logStream, outputStream);

    // Here we assume that Δf = ρ_LUMO - ρ_HOMO
    std::vector<double> deltafk_FD_oneFile = becke.getRhoLumoMinusRhoHomo(beckeParams[0], beckeParams[1], beckeParams[2]);

    logStream << "Δf (column)" << std::endl;
    logStream << std::scientific;
    logStream << std::setprecision(10);
    for (size_t i = 0; i < deltafk_FD_oneFile.size(); ++i)
    {
        logStream << std::right << std::setw(17) << deltafk_FD_oneFile[i] << std::endl;
    }
    logStream << std::endl;
    log(logStream, outputStream);

    // Get dk coefficient
    std::vector<double> dk_coeffs_FD_oneFile(deltafk_FD_oneFile.size(), 0.0);
    for (size_t i = 0; i < deltafk_FD_oneFile.size(); ++i)
    {
        for (size_t j = 0; j < deltafk_FD_oneFile.size(); ++j)
        {
            dk_coeffs_FD_oneFile[i] += deltafk_FD_oneFile[j] * eigenvectors_becke[j][i];
        }
    }

    logStream << "Δf = ";
    for (size_t i = 0; i < deltafk_FD_oneFile.size(); ++i)
    {
        logStream << (dk_coeffs_FD_oneFile[i] >= 0 ? "+ " : "- ") << std::abs(dk_coeffs_FD_oneFile[i]) << " χ_" << i + 1 << ' ';
    }
    logStream << std::endl;
    log(logStream, outputStream);

    double energy_FD_oneFile = 0.0;
    for (size_t i = 0; i < dk_coeffs_FD_oneFile.size(); ++i)
    {
        if (std::abs(eigenvalues_becke[i]) > 1e-5)
        {
            energy_FD_oneFile += dk_coeffs_FD_oneFile[i] * dk_coeffs_FD_oneFile[i] / std::abs(eigenvalues_becke[i]);
        }
    }
    energy_FD_oneFile *= (-0.5);

    logStream << "Energy = " << std::scientific << std::setprecision(10) << energy_FD_oneFile << " H" << std::endl;
    log(logStream, outputStream);

    if (verbose != 0 && logFile)
    {
        logFile.close();
    }
}