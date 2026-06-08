#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <Common/Descriptors.h>
#include <Orbitals/Orbitals.h>
#include <JobControl/Job.h>
#include <JobControl/Jobs/ComputeDescriptors.hpp>
#include <Utils/Enums.hpp>
#include <Utils/Utils.h>


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//----------------------------------------------------------------------------------------------------//

ComputeDescriptors::ComputeDescriptors(const std::string& inputFileName):
    Job(inputFileName)
{ }


//----------------------------------------------------------------------------------------------------//
// STATIC METHODS
//----------------------------------------------------------------------------------------------------//

Descriptors ComputeDescriptors::computeWithIonisationEnergyAndElectronAffinity(const std::string& gridFileName1, const std::string& gridFileName2, const std::string& gridFileName3, double ionisationEnergy, double electronAffinity, PartitionMethod partitionMethod)
{
    std::ifstream gridFile1(gridFileName1);
    if (!gridFile1)
    {
        print_error("Error in ComputeDescriptors::computeWithIonisationEnergyAndElectronAffinity(): could not read file " + gridFileName1 + ".");
        std::exit(1);
    }

    std::ifstream gridFile2(gridFileName2);
    if (!gridFile2)
    {
        print_error("Error in ComputeDescriptors::computeWithIonisationEnergyAndElectronAffinity(): could not read file " + gridFileName2 + ".");
        std::exit(1);
    }

    std::ifstream gridFile3(gridFileName3);
    if (!gridFile3)
    {
        print_error("Error in ComputeDescriptors::computeWithIonisationEnergyAndElectronAffinity(): could not read file " + gridFileName3 + ".");
        std::exit(1);
    }

    Descriptors D(gridFile1, gridFile2, gridFile3, ionisationEnergy, electronAffinity, partitionMethod);
    gridFile1.close();
    gridFile2.close();
    gridFile3.close();

    return D;
}

Descriptors ComputeDescriptors::computeWithEnergies(const std::string& gridFileName1, const std::string& gridFileName2, const std::string& gridFileName3, const std::vector<double>& energies, PartitionMethod partitionMethod)
{
    std::ifstream gridFile1(gridFileName1);
    if (!gridFile1)
    {
        print_error("Error in ComputeDescriptors::computeWithEnergies(): could not read file " + gridFileName1 + ".");
        std::exit(1);
    }

    std::ifstream gridFile2(gridFileName2);
    if (!gridFile2)
    {
        print_error("Error in ComputeDescriptors::computeWithEnergies(): could not read file " + gridFileName2 + ".");
        std::exit(1);
    }

    std::ifstream gridFile3(gridFileName3);
    if (!gridFile3)
    {
        print_error("Error in ComputeDescriptors::computeWithEnergies(): could not read file " + gridFileName3 + ".");
        std::exit(1);
    }

    Descriptors D(gridFile1, gridFile2, gridFile3, energies, partitionMethod);
    gridFile1.close();
    gridFile2.close();
    gridFile3.close();

    return D;
}

Descriptors ComputeDescriptors::computeWithFiniteDifferenceBecke(const std::string& ANAFileName1, const std::string& ANAFileName2, const std::string& ANAFileName3, int kmax, int lebedev_order, int radial_grid_factor)
{
    std::cout << "Building Becke object from " << ANAFileName1 << "... ";
    Becke B1;
    computeOrbitalsOrBecke<Becke>(B1, ANAFileName1);

    std::cout << "Building Becke object from " << ANAFileName2 << "... ";
    Becke B2;
    computeOrbitalsOrBecke<Becke>(B2, ANAFileName2);

    std::cout << "Building Becke object from " << ANAFileName3 << "... ";
    Becke B3;
    computeOrbitalsOrBecke<Becke>(B3, ANAFileName3);

    const std::vector<std::vector<double>>& partialChargesAndEnergy1 = B1.PartialChargesAndEnergy(kmax, lebedev_order, radial_grid_factor);
    const std::vector<std::vector<double>>& partialChargesAndEnergy2 = B2.PartialChargesAndEnergy(kmax, lebedev_order, radial_grid_factor);
    const std::vector<std::vector<double>>& partialChargesAndEnergy3 = B3.PartialChargesAndEnergy(kmax, lebedev_order, radial_grid_factor);

    std::vector<double> energies;
    energies.push_back(partialChargesAndEnergy1[0][0]);
    energies.push_back(partialChargesAndEnergy2[0][0]);
    energies.push_back(partialChargesAndEnergy3[0][0]);

    return Descriptors(B1.get_molecule(), partialChargesAndEnergy1[1], partialChargesAndEnergy2[1], partialChargesAndEnergy3[1], energies);
}


//----------------------------------------------------------------------------------------------------//
// PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

void ComputeDescriptors::run()
{
    // Read partition method
    PartitionMethod partitionMethod;
    readPartitionMethod(partitionMethod);


    // Check partition method validity
    if (partitionMethod == PartitionMethod::BBS || partitionMethod == PartitionMethod::B2S)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: partitionMethod \"" << to_string(partitionMethod) << "\" invalid for this job." << std::endl;
        errorMessage << "Please check documentation and \"PartitionMethod\" parameter value in " << _inputFileName << '.';

        print_error(errorMessage.str());
        
        std::exit(1);
    }

    // Print partition method
    std::cout << "Volume partition method: " << to_string(partitionMethod) << std::endl << std::endl;

    if (partitionMethod == PartitionMethod::FD || partitionMethod == PartitionMethod::FMO)
    {
        // Read analytic files names
        std::vector<std::string> analyticFilesNames;
        readAnalyticFilesNames(analyticFilesNames);

        if (partitionMethod == PartitionMethod::FD)
        {
            // Check number of analytic files names
            if (analyticFilesNames.size() != 3)
            {
                std::stringstream errorMessage;
                errorMessage << "Error: incorrect number of analytic files names (three files expected)." << std::endl;
                errorMessage << "Please check the documentation and the number of files specified in the \"AnalyticFiles\" parameter in " << _inputFileName << '.';

                print_error(errorMessage.str());

                std::exit(1);
            }
            
            // Compute descriptors
            Descriptors descriptors = computeWithFiniteDifferenceBecke(analyticFilesNames[0], analyticFilesNames[1], analyticFilesNames[2]);
            std::cout << descriptors << std::endl;
        }
        else
        {
            // Check number of analytic files names
            if (analyticFilesNames.size() > 1)
            {
                std::stringstream errorMessage;
                errorMessage << "Error: too many analytic files names (one file expected)." << std::endl;
                errorMessage << "Please check the documentation and the number of files specified in the \"AnalyticFiles\" parameter in " << _inputFileName << '.';

                print_error(errorMessage.str());

                std::exit(1);
            }


            // Compute descriptors
            Orbitals o;
            computeOrbitalsOrBecke<Orbitals>(o, analyticFilesNames[0]);
            o.printDescriptors();
        }
    }
    else
    {
        // Read grid files names
        std::vector<std::string> gridFilesNames;
        readGridFilesNames(gridFilesNames);

        // Check number of grid files names
        if (gridFilesNames.size() != 3)
        {
            std::stringstream errorMessage;
            errorMessage << "Error: incorrect number of grid files names (three files expected)." << std::endl;
            errorMessage << "Please check the documentation and the number of files specified in the \"GridFiles\" parameter in " << _inputFileName << '.';

            print_error(errorMessage.str());

            std::exit(1);
        }


        // Read energies
        std::vector<double> energies;
        readEnergies(energies);

        if (energies.size() == 2)
        {
            std::cout << "Reading Ionization potential I = " << energies[0] << " and electron Affinity A = " << energies[1] << std::endl;
            
            // Compute descriptors
            Descriptors D = computeWithIonisationEnergyAndElectronAffinity(gridFilesNames[0], gridFilesNames[1], gridFilesNames[2], energies[0], energies[1], partitionMethod);
            std::cout << D;
        }
        else if (energies.size() == 3)
        {
            std::cout << " Reading Total Energies: E1 = " << energies[0] << ", E2 = " << energies[1] << " and E3 = " << energies[2] << std::endl;
            Descriptors D = computeWithEnergies(gridFilesNames[0], gridFilesNames[1], gridFilesNames[2], energies, partitionMethod);
            std::cout << D;
        }
        else
        {
            std::stringstream errorMessage;
            errorMessage << "Error: incorrect number of energies (two or three expected)." << std::endl;
            errorMessage << "Please check the documentation and the number of energies specified in the \"Energies\" parameter in " << _inputFileName << '.';

            print_error(errorMessage.str());

            std::exit(1);
        }
    }
}
