#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifdef ENABLE_OMP
#include <omp.h>
#endif

#include <Becke/Becke.h>
#include <Common/Descriptors.h>
#include <Common/Element.h>
#include <Common/Constants.h>
#include <Cube/Grid.h>
#include <Cube/GridCP.h>
#include <JobControl/Job.h>
#include <JobControl/Jobs/ComputeDescriptors.hpp>
#include <Orbitals/Orbitals.h>
#include <Orbitals/ExcitedState.hpp>
#include <Orbitals/SlaterDeterminant.hpp>
#include <Utils/Utils.h>


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTORS AND DESTRUCTOR
//----------------------------------------------------------------------------------------------------//

Job::Job():
    _inputFileName(""),
    _inputFile()
{ }

Job::Job(const std::string& inputFileName):
    _inputFileName(inputFileName)
{
    openInputFile();
}

Job::~Job()
{
    _inputFile.close();
}


//----------------------------------------------------------------------------------------------------//
// SPECIFIC PARAMETERS READING FROM INPUT FILE
//----------------------------------------------------------------------------------------------------//

bool Job::readAnalyticFilesNames(std::vector<std::string>& analyticFilesNames)
{
    bool read = readListType<std::string>(_inputFile, "AnalyticFiles", analyticFilesNames);

    if (!read)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: could not find analytic files names." << std::endl;
        errorMessage << "Please check if the \"AnalyticFiles\" parameter is defined and set in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }

    // Checking if the files have an extension
    std::regex fileExtensionRegex("\\.([a-zA-Z0-9]+)$");
    std::smatch fileExtensionMatch;

    for (const std::string &fileName : analyticFilesNames)
    {
        if (!std::regex_search(fileName, fileExtensionMatch, fileExtensionRegex))
        {
            std::stringstream errorMessage;
            errorMessage << "Error: cannot determine file format from file name \"" << fileName << "\" (no file extension)." << std::endl;
            errorMessage << "Please check the documentation and the \"AnalyticFiles\" parameter value in the provided input file (" << _inputFileName << ").";

            print_error(errorMessage.str());

            std::exit(1);
        }
    }

    return read;
}

bool Job::readBecke(std::vector<int>& beckeParameters)
{
    bool read = readListType<int>(_inputFile, "Becke", beckeParameters);

    if (read && !beckeParameters.empty() && beckeParameters.size() != 3)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect number of values for the \"Becke\" parameter (three values expected)." << std::endl;
        errorMessage << "Please check documentation and the \"Becke\" parameter values in " << _inputFileName << '.';

        print_error(errorMessage.str());

        std::exit(1);
    }
    else if (!read)
    {
        std::string defaultBeckeStr;
        read = readOneString(_inputFile, "Becke", defaultBeckeStr);

        if (to_lower(defaultBeckeStr) == "default")
        {
            beckeParameters = { 3, 41, 5 };
        }
        else
        {
            std::stringstream errorMessage;
            errorMessage << "Error: incorrect value for the \"Becke\" parameter (either three integers or \"Default\" string expected)." << std::endl;
            errorMessage << "Please check documentation and the \"Becke\" parameter values in " << _inputFileName << '.';

            print_error(errorMessage.str());

            std::exit(1);
        }
    }

    return read;
}

bool Job::readCharges(std::vector<double>& charges)
{
    bool read = readListType<double>(_inputFile, "Charges", charges);

    if (!read)
    {
        std::cout << "Note: the \"Charges\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (Charges = -1)." << std::endl << std::endl;

        charges = { - 1.0 };
    }

    return read;
}

bool Job::readChargesPositionsBijections(bool& chargesPositionsBijections)
{
    std::string strChargesPositionsBijections;
    bool read = readOneString(_inputFile, "ChargesPositionsBijections", strChargesPositionsBijections);

    chargesPositionsBijections = false;
    if (!read)
    {
        std::cout << "Note: the \"ChargesPositionsBijections\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (ChargesPositionsBijections = False)." << std::endl << std::endl;
    }
    else if (to_lower(strChargesPositionsBijections) == "true")
    {
        chargesPositionsBijections = true;
    }
    else if (to_lower(strChargesPositionsBijections) != "false")
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect value for the \"ChargesPositionsBijections\" parameter (" << strChargesPositionsBijections << ")." << std::endl;
        errorMessage << "Please check the documentation and the \"ChargesPositionsBijections\" parameter value in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }

    return read;
}

bool Job::readCutoff(double& cutoff)
{
    bool read = readOneType<double>(_inputFile, "Cutoff", cutoff);

    if (!read)
    {
        std::cout << "Note: the \"Cutoff\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (Cutoff = 0.0)." << std::endl << std::endl;

        cutoff = 0.0;
    }

    return read;
}

bool Job::readELFMethod(ELFMethod& elfMethod)
{
    std::string strELFMethod;
    bool read = readOneString(_inputFile, "ELFMethod", strELFMethod);

    if (!read)
    {
        std::cout << "Note: the \"ELFMethod\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (ELFMethod = Savin)." << std::endl << std::endl;
        elfMethod = ELFMethod::SAVIN;
    }
    else
    {
        elfMethod = elfMethod_from_string(strELFMethod);
    }

    // Handle unknown ELF method
    if (elfMethod == ELFMethod::UNKNOWN)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: ELF method \"" << strELFMethod << "\" unknown." << std::endl;
        errorMessage << "Please check the documentation and the \"ELFMethod\" parameter value in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }

    return read;
}

bool Job::readEnergies(std::vector<double>& energies)
{
    bool read = readListType<double>(_inputFile, "Energies", energies);

    if (!read)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: could not find energies." << std::endl;
        errorMessage << "Please check if the \"Energies\" parameter is defined and set in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }

    return read;
}

bool Job::readEnergyPointChargeMethods(std::vector<EnergyPointChargeMethod>& energyPointChargeMethods)
{
    std::vector<std::string> strEnergyPointChargeMethods;
    bool read = readListType<std::string>(_inputFile, "EnergyPointChargeMethods", strEnergyPointChargeMethods);

    if (!read)
    {
        std::cout << "Note: the \"EnergyPointChargeMethods\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (EnergyPointChargeMethods = Variational)." << std::endl << std::endl;

        energyPointChargeMethods = { EnergyPointChargeMethod::VARIATIONAL };
    }

    for (const std::string& strMethod : strEnergyPointChargeMethods)
    {
        EnergyPointChargeMethod method = energyPointChargeMethod_from_string(strMethod);

        // Handle unknown method: exit program with error message.
        if (method == EnergyPointChargeMethod::UNKNOWN)
        {
            std::stringstream errorMessage;
            errorMessage << "Error: Energy point charge method \"" << strMethod << "\" unknown." << std::endl;
            errorMessage << "Please check the documentation and the \"EnergyPointChargeMethods\" parameter values in the provided input file (" << _inputFileName << ").";

            print_error(errorMessage.str());

            std::exit(1);
        }

        energyPointChargeMethods.push_back(method);
    }

    return read;
}

bool Job::readExcitedStatesNumbers(std::vector<int>& excitedStatesNumbers)
{
    bool read = readListType<int>(_inputFile, "ExcitedStatesNumbers", excitedStatesNumbers);

    if (!read)
    {
        std::stringstream errorMessage;
        errorMessage << "Note: the \"ExcitedStatesNumbers\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        errorMessage << "The program will keep all excited states." << std::endl << std::endl;
    }

    return read;
}

bool Job::readExcludedOrbitals(std::vector<int>& excludedOrbitals)
{
    bool read = readListType<int>(_inputFile, "ExcludedOrbitals", excludedOrbitals);

    if (!read)
    {
        std::stringstream errorMessage;
        errorMessage << "Note: the \"ExcludedOrbitalsNumbers\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        errorMessage << "The program will keep all orbitals for the computation." << std::endl << std::endl;
    }

    return read;
}

bool Job::readGridFilesNames(std::vector<std::string>& gridFilesNames)
{
    bool read = readListType<std::string>(_inputFile, "Grids", gridFilesNames);

    if (!read)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: could not find grid files names." << std::endl;
        errorMessage << "Please check if the \"Grids\" parameter is defined and set in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }

    // Checking if the files have an extension
    std::regex fileExtensionRegex("\\.([a-zA-Z0-9]+)$");
    std::smatch fileExtensionMatch;
    
    for (const std::string& fileName : gridFilesNames)
    {
        if (!std::regex_search(fileName, fileExtensionMatch, fileExtensionRegex))
        {
            std::stringstream errorMessage;
            errorMessage << "Error: cannot determine file format from file name \"" << fileName << "\" (no file extension)." << std::endl;
            errorMessage << "Please check the documentation and the \"Grids\" parameter value in the provided input file (" << _inputFileName << ").";

            print_error(errorMessage.str());

            std::exit(1);
        }
    }

    return read;
}

bool Job::readGroundStateEnergy(double& energy)
{
    return readOneType<double>(_inputFile, "GroundStateEnergy", energy);
}

bool Job::readMaxNumberOfExcitedStates(int& maxNumberOfExcitedStates)
{
    bool read = readOneType<int>(_inputFile, "MaxNumberOfExcitedStates", maxNumberOfExcitedStates);

    if (!read)
    {
        std::cout << "Note: the \"MaxNumberOfExcitedStates\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will consider all excited states available in the excited states file (MaxNumberOfExcitedStates = -1)." << std::endl << std::endl;

        maxNumberOfExcitedStates = -1;
    }

    return read;
}

bool Job::readNuclearCutoff(double& nuclearCutoff)
{
    bool read = readOneType<double>(_inputFile, "NuclearCutoff", nuclearCutoff);

    if (!read)
    {
        std::cout << "Note: the \"NuclearCutoff\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (NuclearCutoff = 0.1 Å)." << std::endl << std::endl;

        nuclearCutoff = 0.1 * Constants::ANGSTROM_TO_BOHR_RADIUS;
    }
    else
    {
        nuclearCutoff *= Constants::ANGSTROM_TO_BOHR_RADIUS;
    }

    return read;
}

bool Job::readOrbitalsNumbers(std::vector<int>& orbitalsNumbers)
{
    bool read = readListType<int>(_inputFile, "OrbitalsNumbers", orbitalsNumbers);

    if (!read)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: could not find orbitals numbers." << std::endl;
        errorMessage << "Please check if the \"OrbitalsNumbers\" parameter is defined and set in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }

    return read;
}

bool Job::readOrbitalsSpins(std::vector<SpinType>& orbitalsSpins)
{
    std::vector<std::string> strOrbitalsSpins;
    bool read = readListType<std::string>(_inputFile, "OrbitalsSpins", strOrbitalsSpins);

    if (!read)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: could not find orbitals spins." << std::endl;
        errorMessage << "Please check if the \"OrbitalsSpins\" parameter is defined and set in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }

    for (const std::string& strSpin : strOrbitalsSpins)
    {
        SpinType spin = spinType_from_string(strSpin);

        // Handle unknown spin type: exit program with error message.
        if (spin == SpinType::UNKNOWN)
        {
            std::stringstream errorMessage;
            errorMessage << "Error: Orbital spin type \"" << strSpin << "\" unknown." << std::endl;
            errorMessage << "Please check the documentation and the \"OrbitalsSpins\" parameter value in the provided input file (" << _inputFileName << ").";

            print_error(errorMessage.str());

            std::exit(1);
        }

        // Handle forbidden ALPHA_BETA spin type in the orbitals spin list.
        if (spin == SpinType::ALPHA_BETA)
        {
            std::stringstream errorMessage;
            errorMessage << "Error: Spin type \"Alpha-Beta\" cannot be used in a custom spin list." << std::endl;
            errorMessage << "Please check the documentation and the \"OrbitalsSpins\" parameter values in the provided input file (" << _inputFileName << ").";

            print_error(errorMessage.str());

            std::exit(1);
        }

        orbitalsSpins.push_back(spin);
    }

    return read;
}

bool Job::readOutputPrefix(std::string& outputPrefix)
{
    bool read = readOneString(_inputFile, "OutputPrefix", outputPrefix);

    if (!read)
    {
        std::cout << "Note: the \"OutputPrefix\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (OutputPrefix = \"\")." << std::endl << std::endl;

        outputPrefix = "";
    }

    if (!outputPrefix.empty() && outputPrefix.back() != '_')
    {
        outputPrefix += '_';
    }

    return read;
}

bool Job::readOrbitalType(OrbitalType& orbitalType)
{
    std::string strOrbitalType;
    bool read = readOneString(_inputFile, "OrbitalType", strOrbitalType);

    if (!read)
    {
        std::cout << "Note: the \"OrbitalType\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (OrbitalType = All)." << std::endl << std::endl;
        orbitalType = OrbitalType::ALL;
    }
    else
    {
        orbitalType = orbitalType_from_string(strOrbitalType);
    }

    // Handle unknown orbital type: exit program with error message.
    if (orbitalType == OrbitalType::UNKNOWN)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: Orbital type \"" << strOrbitalType << "\" unknown." << std::endl;
        errorMessage << "Please check the documentation and the \"OrbitalType\" parameter value in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }

    return read;
}

bool Job::readPartitionMethod(PartitionMethod& partitionMethod)
{
    std::string strMethod;
    bool read = readOneString(_inputFile, "PartitionMethod", strMethod);

    if (!read)
    {
        std::cout << "Note: the \"PartitionMethod\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (PartitionMethod = On-Grid)." << std::endl << std::endl;
        partitionMethod = PartitionMethod::AIM_ON_GRID;
    }
    else
    {
        partitionMethod = partitionMethod_from_string(strMethod);
    }

    // Handle unknown partition method: exit program with error message.
    if (partitionMethod == PartitionMethod::UNKNOWN)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: Volume partition method \"" << strMethod << "\" unknown." << std::endl;
        errorMessage << "Please check the documentation and the \"PartitionMethod\" parameter value in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }

    return read;
}

bool Job::readPositions(std::vector<std::array<double, 3>>& positions)
{
    bool ok = readListTypeArray<double, 3>(_inputFile, "Positions", positions);

    for (size_t i = 0; i < positions.size(); ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            positions[i][j] *= Constants::ANGSTROM_TO_BOHR_RADIUS;
        }
    }

    return ok;
}

bool Job::readPrecision(int& precision)
{
    bool read = readOneType<int>(_inputFile, "Precision", precision);

    if (!read)
    {
        std::cout << "Note: the \"Precision\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (Precision = 10)." << std::endl << std::endl;

        precision = 10;
    }

    return read;
}

bool Job::readRDMMethod(RDMMethod& rdmMethod)
{
    std::string strRDMMethod;
    bool read = readOneString(_inputFile, "RDMMethod", strRDMMethod);

    if (!read)
    {
        std::cout << "Note: the \"RDMMethod\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (RDMMethod = Gamma)." << std::endl << std::endl;
        rdmMethod = RDMMethod::GAMMA;
    }
    else
    {
        rdmMethod = rdmMethod_from_string(strRDMMethod);
    }

    // Handle unknown RDM method: exit program with error message.
    if (rdmMethod == RDMMethod::UNKNOWN)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: RDM method \"" << strRDMMethod << "\" unknown." << std::endl;
        errorMessage << "Please check the documentation and the \"RDMMethod\" parameter value in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }

    return read;
}

bool Job::readRunType(RunType& runType)
{
    std::string strRunType;
    bool read = readOneString(_inputFile, "RunType", strRunType);

    if (!read)
    {
        std::cout << "Note: the \"RunType\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (RunType = Help)." << std::endl << std::endl;
        runType = RunType::HELP;
    }
    else
    {
        runType = runType_from_string(strRunType);
    }

    // Handle unknown run type: exit program with error message.
    if (runType == RunType::UNKNOWN)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: Run type \"" << strRunType << "\" unknown." << std::endl;
        errorMessage << "Please check the \"RunType\" parameter value in " << _inputFileName << '.';

        print_error(errorMessage.str());

        std::exit(1);
    }

    return read;
}

bool Job::readSavePseudoOrbitals(bool& savePseudoOrbitals)
{
    std::string strSavePseudoOrbitals;
    bool read = readOneString(_inputFile, "SavePseudoOrbitals", strSavePseudoOrbitals);
    
    savePseudoOrbitals = false;
    if (!read)
    {
        std::cout << "Note: the \"SavePseudoOrbitals\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (SavePseudoOrbitals = False)." << std::endl << std::endl;
    }
    else if (to_lower(strSavePseudoOrbitals) == "true")
    {
        savePseudoOrbitals = true;
    }
    else if (to_lower(strSavePseudoOrbitals) != "false")
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect value for the \"SavePseudoOrbitals\" parameter (" << strSavePseudoOrbitals << ")." << std::endl;
        errorMessage << "Please check the documentation and the \"SavePseudoOrbitals\" parameter value in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }

    return read;
}

bool Job::readSaveReducedDensityMatrix(bool& saveReducedDensityMatrix)
{
    std::string strSaveReducedDensityMatrix;
    bool read = readOneString(_inputFile, "SaveReducedDensityMatrix", strSaveReducedDensityMatrix);
    
    saveReducedDensityMatrix = false;
    if (!read)
    {
        std::cout << "Note: the \"SaveReducedDensityMatrix\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (SaveReducedDensityMatrix = False)." << std::endl << std::endl;
    }
    else if (to_lower(strSaveReducedDensityMatrix) == "true")
    {
        saveReducedDensityMatrix = true;
    }
    else if (to_lower(strSaveReducedDensityMatrix) != "false")
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect value for the \"SaveReducedDensityMatrix\" parameter (" << strSaveReducedDensityMatrix << ")." << std::endl;
        errorMessage << "Please check the documentation and the \"SaveReducedDensityMatrix\" parameter value in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }

    return read;
}

bool Job::readShowProgress(bool& showProgress)
{
    std::string strShowProgress;
    bool read = readOneString(_inputFile, "ShowProgress", strShowProgress);

    showProgress = false;
    if (!read)
    {
        std::cout << "Note: the \"ShowProgress\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (ShowProgress = False)." << std::endl << std::endl;
    }
    else if (to_lower(strShowProgress) == "true")
    {
        showProgress = true;
    }
    else if (to_lower(strShowProgress) != "false")
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect value for the \"ShowProgress\" parameter (" << strShowProgress << ")." << std::endl;
        errorMessage << "Please check the documentation and the \"ShowProgress\" parameter value in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }

    return read;
}

bool Job::readSingleCharge(bool& singleCharge)
{
    std::string strSingleCharge;
    bool read = readOneString(_inputFile, "SingleCharge", strSingleCharge);

    singleCharge = true;
    if (!read)
    {
        std::cout << "Note: the \"SingleCharge\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (SingleCharge = True)." << std::endl << std::endl;
    }
    else if (to_lower(strSingleCharge) == "false")
    {
        singleCharge = false;
    }
    else if (to_lower(strSingleCharge) != "true")
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect value for the \"SingleCharge\" parameter (" << strSingleCharge << ")." << std::endl;
        errorMessage << "Please check the documentation and the \"SingleCharge\" parameter value in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }

    return read;
}

bool Job::readSize(GridSize& gridSize, CustomSizeData& customSizeData)
{
    std::string strGridSize;
    bool read = readOneString(_inputFile, "Size", strGridSize);

    if (!read)
    {
        std::cout << "Note: the \"Size\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (Size = Medium)." << std::endl << std::endl;

        gridSize = GridSize::MEDIUM;
    }
    else
    {
        gridSize = gridSize_from_string(strGridSize);
    }

    // Handle unknown grid size: exit program with error message.
    if (gridSize == GridSize::UNKNOWN)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: Grid size \"" << strGridSize << "\" unknown." << std::endl;
        errorMessage << "Please check the documentation and the \"Size\" parameter value in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }

    if (gridSize == GridSize::CUSTOM)
    {
        std::vector<double> customSizeDataValues;
        if (!readListType<double>(_inputFile, "CustomSizeData", customSizeDataValues))
        {
            std::stringstream errorMessage;
            errorMessage << "Error: could not find custom size data." << std::endl;
            errorMessage << "Please check if the \"CustomSizeData\" parameter is defined and set in the provided input file (" << _inputFileName << ").";

            print_error(errorMessage.str());

            std::exit(1);
        }

        if (customSizeDataValues.size() != 15)
        {
            std::stringstream errorMessage;
            errorMessage << "Error: incorrect number of values for the \"CustomSizeData\" parameter (fifteen values expected)." << std::endl;
            errorMessage << "Please check documentation and the \"CustomSizeData\" parameter values in " << _inputFileName << '.';

            print_error(errorMessage.str());

            std::exit(1);
        }

        customSizeData = CustomSizeData(customSizeDataValues);
    }

    return read;
}

bool Job::readSpinList(std::vector<SpinType>& spinList)
{
    std::vector<std::string> strSpinList;
    bool read = readListType<std::string>(_inputFile, "SpinList", strSpinList);

    if (!read)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: could not find spin list." << std::endl;
        errorMessage << "Please check if the \"SpinList\" parameter is defined and set in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }

    for (const std::string& strSpin : strSpinList)
    {
        SpinType spin = spinType_from_string(strSpin);

        // Handle unknown spin type: exit program with error message.
        if (spin == SpinType::UNKNOWN)
        {
            std::stringstream errorMessage;
            errorMessage << "Error: Orbital spin type \"" << strSpin << "\" unknown." << std::endl;
            errorMessage << "Please check the documentation and the \"SpinList\" parameter values in the provided input file (" << _inputFileName << ").";

            print_error(errorMessage.str());

            std::exit(1);
        }

        // Handle forbidden ALPHA_BETA spin type in spin list.
        if (spin == SpinType::ALPHA_BETA)
        {
            std::stringstream errorMessage;
            errorMessage << "Error: Spin type \"Alpha-Beta\" cannot be used in a custom spin list." << std::endl;
            errorMessage << "Please check the documentation and the \"SpinList\" parameter values in the provided input file (" << _inputFileName << ").";

            print_error(errorMessage.str());

            std::exit(1);
        }

        spinList.push_back(spin);
    }

    return read;
}

bool Job::readSpinType(SpinType& spinType)
{
    std::string strSpinType;
    bool read = readOneString(_inputFile, "SpinType", strSpinType);

    if (!read)
    {
        std::cout << "Note: the \"SpinType\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (SpinType = Alpha-Beta)." << std::endl << std::endl;

        spinType = SpinType::ALPHA_BETA;
    }
    else
    {
        spinType = spinType_from_string(strSpinType);
    }

    // Handle unknown spin type: exit program with error message.
    if (spinType == SpinType::UNKNOWN)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: Spin type \"" << strSpinType << "\" unknown." << std::endl;
        errorMessage << "Please check the documentation and the \"SpinType\" parameter value in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }

    return read;
}

bool Job::readTransitionDensities(std::vector<std::array<int, 2>>& transitionDensities)
{
    bool read = readListTypeArray<int, 2>(_inputFile, "TransitionDensities", transitionDensities);

    if (!read)
    {
        std::cout << "Note: the \"TransitionDensities\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will consider all transition densities." << std::endl << std::endl;
    }

    return read;
}

bool Job::readTransitionsFileName(std::string& transitionsFileName)
{
    bool read = readOneString(_inputFile, "TransitionsFile", transitionsFileName);

    if (!read)
    {
        std::cout << "Note: the \"TransitionsFile\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will try to read the transitions from the analytic file." << std::endl << std::endl;

        transitionsFileName = "";
    }

    return read;
}

bool Job::readVerbose(int& verbose)
{
    bool read = readOneType<int>(_inputFile, "Verbose", verbose);

    if (!read)
    {
        verbose = 0;
    }

    return read;
}


//----------------------------------------------------------------------------------------------------//
// OTHER PRIVATE METHODS
//----------------------------------------------------------------------------------------------------//

Domain Job::buildDomainForCube(Orbitals& orb, const GridSize gridSize, const CustomSizeData& customSizeData, const int& Nval)
{
    Domain d;
    double size = d.sizeUpMol(orb.get_struct(), 1.5);
    double xmax = size;
    double ymax = size;
    double zmax = size;
    int N1 = 0;
    int N2 = 0;
    int N3 = 0;

    if (gridSize == GridSize::COARSE)
    {
        N1 = static_cast<int>(std::floor(size * 6));
        N2 = static_cast<int>(std::floor(size * 6));
        N3 = static_cast<int>(std::floor(size * 6));
    }
    else if (gridSize == GridSize::MEDIUM)
    {
        N1 = static_cast<int>(std::floor(size * 12));
        N2 = static_cast<int>(std::floor(size * 12));
        N3 = static_cast<int>(std::floor(size * 12));
    }
    else if (gridSize == GridSize::FINE)
    {
        N1 = static_cast<int>(std::floor(size * 24));
        N2 = static_cast<int>(std::floor(size * 24));
        N3 = static_cast<int>(std::floor(size * 24));
    }

    if (gridSize == GridSize::COARSE || gridSize == GridSize::MEDIUM || gridSize == GridSize::FINE)
    {
        std::array<double, 3> tx = {2.0 * xmax / (N1 - 1), 0.0, 0.0};
        std::array<double, 3> ty = {0.0, 2.0 * ymax / (N2 - 1), 0.0};
        std::array<double, 3> tz = {0.0, 0.0, 2.0 * zmax / (N3 - 1)};

        std::array<std::array<double, 3>, 3>  t = {tx, ty, tz};

        d.set_all(Nval, N1, N2, N3, xmax, ymax, zmax, t);
    }
    else
    {
        std::array<double, 3> dx = {customSizeData[6], customSizeData[7], customSizeData[8]};
        std::array<double, 3> dy = {customSizeData[9], customSizeData[10], customSizeData[11]};
        std::array<double, 3> dz = {customSizeData[12], customSizeData[13], customSizeData[14]};

        std::array<std::array<double, 3>, 3>  t = {dx, dy, dz};

        d.set_all(Nval, int(customSizeData[0]), int(customSizeData[1]), int(customSizeData[2]), customSizeData[3], customSizeData[4], customSizeData[5], t);
    }
    return d;
}

Orbitals Job::computePseudoOrbitalsFromLrfMatrix(const Orbitals& orbitals, const std::vector<std::vector<std::vector<double>>>& lrfMatrix, std::vector<std::vector<double>>& eigenvalues, std::vector<std::vector<std::vector<double>>>& eigenvectors, const std::string& outputPrefix, bool savePseudoOrbitals, std::ostream& outputStream, int verbose, bool showProgress)
{
    std::stringstream logStream;


    // Diagonalize LRF matrixes for both spin
    eigenvalues.resize(2);
    eigenvectors.resize(2);


    // Compute and save results for alpha spin
    findEigenValuesAndEigenVectorsOfSymmetricalMatrix(lrfMatrix[0], eigenvalues[0], eigenvectors[0]);

    std::ofstream outputFile(outputPrefix + "_eigenvalues_alpha.cdftt");
    if (!outputFile)
    {
        std::stringstream errorMessage;
        errorMessage << "Error in Job::computePseudoOrbitalsLinearResponseWithPointCharges(): could not open output file " << outputPrefix << "_eigenvalues_alpha.cdftt for writing." << std::endl;

        print_error(errorMessage.str());
        
        std::exit(1);
    }

    logStream << "Sorted Eigenvalues (alpha spin):" << std::endl;
    log(logStream, outputStream);
    logStream << std::scientific;
    logStream << std::setprecision(10);
    outputFile << std::scientific;
    outputFile << std::setprecision(10);
    for (size_t k = 0; k < eigenvalues[0].size(); ++k)
    {
        logStream << eigenvalues[0][k] << ' ';
        outputFile << eigenvalues[0][k] << std::endl;
    }
    logStream << std::endl << std::endl;
    log(logStream, outputStream);
    outputFile.close();

    outputFile.open(outputPrefix + "_eigenvectors_alpha.cdftt");
    if (!outputFile)
    {
        std::stringstream errorMessage;
        errorMessage << "Error in Job::computePseudoOrbitalsLinearResponseWithPointCharges(): could not open output file " << outputPrefix << "_eigenvectors_alpha.cdftt for writing." << std::endl;

        print_error(errorMessage.str());
        
        std::exit(1);
    }

    logStream << "Sorted Eigenvectors (columns, alpha spin): " << std::endl;
    log(logStream, outputStream);
    logStream << std::scientific;
    logStream << std::setprecision(10);
    outputFile << std::scientific;
    outputFile << std::setprecision(10);
    for (size_t i = 0; i < eigenvectors[0].size(); ++i)
    {
        for (size_t j = 0; j < eigenvectors[0][i].size(); ++j)
        {
            logStream << std::right << std::setw(17) << eigenvectors[0][i][j] << '\t';
            outputFile << std::right << std::setw(17) << eigenvectors[0][i][j] << ' ';
        }

        logStream << std::endl;
        outputFile << std::endl;
    }
    logStream << std::defaultfloat << std::endl;
    log(logStream, outputStream);
    outputFile.close();


    // Compute and save results for beta spin
    findEigenValuesAndEigenVectorsOfSymmetricalMatrix(lrfMatrix[1], eigenvalues[1], eigenvectors[1]);

    outputFile.open(outputPrefix + "_eigenvalues_beta.cdftt");
    if (!outputFile)
    {
        std::stringstream errorMessage;
        errorMessage << "Error in Job::computePseudoOrbitalsLinearResponseWithPointCharges(): could not open output file " << outputPrefix << "_eigenvalues_beta.cdftt for writing." << std::endl;

        print_error(errorMessage.str());
        
        std::exit(1);
    }

    logStream << "Sorted Eigenvalues (beta spin):" << std::endl;
    log(logStream, outputStream);
    logStream << std::scientific;
    logStream << std::setprecision(10);
    outputFile << std::scientific;
    outputFile << std::setprecision(10);
    for (size_t k = 0; k < eigenvalues[1].size(); ++k)
    {
        logStream << eigenvalues[1][k] << ' ';
        outputFile << eigenvalues[1][k] << std::endl;
    }
    logStream << std::endl << std::endl;
    log(logStream, outputStream);
    outputFile.close();

    outputFile.open(outputPrefix + "_eigenvectors_beta.cdftt");
    if (!outputFile)
    {
        std::stringstream errorMessage;
        errorMessage << "Error in Job::computePseudoOrbitalsLinearResponseWithPointCharges(): could not open output file " << outputPrefix << "_eigenvectors_beta.cdftt for writing." << std::endl;

        print_error(errorMessage.str());
        
        std::exit(1);
    }

    logStream << "Sorted Eigenvectors (columns, beta spin): " << std::endl;
    log(logStream, outputStream);
    logStream << std::scientific;
    logStream << std::setprecision(10);
    outputFile << std::scientific;
    outputFile << std::setprecision(10);
    for (size_t i = 0; i < eigenvectors[1].size(); ++i)
    {
        for (size_t j = 0; j < eigenvectors[1][i].size(); ++j)
        {
            logStream << std::right << std::setw(17) << eigenvectors[1][i][j] << '\t';
            outputFile << std::right << std::setw(17) << eigenvectors[1][i][j] << ' ';
        }

        logStream << std::endl;
        outputFile << std::endl;
    }
    logStream << std::defaultfloat << std::endl;
    log(logStream, outputStream);
    outputFile.close();


    // Expand LRF eigenvector in AO basis (eigenvectors are in vertical format, i.e. columns are eigenvectors)
    std::vector<std::vector<std::vector<double>>> lrfEigenvectorsInAoBasis(2); // First index for alpha spin, second index for beta spin
    const std::vector<std::vector<std::vector<double>>>& coefficients = orbitals.get_coefficients();

    for (int spin = 0; spin < 2; ++spin)
    {
        lrfEigenvectorsInAoBasis[spin] = std::vector<std::vector<double>>(orbitals.get_numberOfMo(), std::vector<double>(orbitals.get_numberOfAo(), 0.0));

        for (size_t i = 0; i < eigenvectors[spin].size(); ++i) // phi
        {
            for (size_t j = 0; j < eigenvectors[spin][i].size(); ++j) // sigma
            {
                for (size_t k = 0; k < coefficients[spin].size(); ++k) // xi
                {
                    lrfEigenvectorsInAoBasis[spin][j][k] += eigenvectors[spin][i][j] * coefficients[spin][i][k];
                }
            }
        }
    }


    // Copy orbitals to pseudoOrbitals to keep the same structure and only change energies and coefficients
    Orbitals pseudoOrbitals(orbitals);

    // Replace MO energies by LRF eigenvalues
    pseudoOrbitals.set_orbitalEnergy(eigenvalues);

    // Replace coefficients by LRF eigenvectors in AO basis
    pseudoOrbitals.set_coefficients(lrfEigenvectorsInAoBasis);

    // Save pseudoOrbitals in cube format if desired by the user
    if (savePseudoOrbitals)
    {
        // Read cube grid parameters and build domain
        GridSize gridSize;
        CustomSizeData customSizeData;
        readSize(gridSize, customSizeData);
        Domain domain = buildDomainForCube(pseudoOrbitals, gridSize, customSizeData, pseudoOrbitals.get_numberOfMo());

        std::vector<int> pseudoOrbitalsIndexes;
        std::vector<SpinType> pseudoOrbitalsSpinTypes_alpha;
        std::vector<SpinType> pseudoOrbitalsSpinTypes_beta;
        for (int i = 0; i < pseudoOrbitals.get_numberOfMo(); ++i)
        {
            pseudoOrbitalsIndexes.push_back(i);
            pseudoOrbitalsSpinTypes_alpha.push_back(SpinType::ALPHA);
            pseudoOrbitalsSpinTypes_beta.push_back(SpinType::BETA);
        }

        // Save pseudo orbitals
        std::cout << "Building pseudo orbitals grid for alpha spin..." << std::endl;
        createCube(pseudoOrbitals, domain, outputPrefix + "_lrf_pseudoOrbitals_alpha.cube", CubeType::ORBITALS, showProgress, ELFMethod::UNKNOWN, pseudoOrbitalsIndexes, pseudoOrbitalsSpinTypes_alpha);
        logStream << "Pseudo orbitals with alpha spin saved to " << outputPrefix << "_lrf_pseudoOrbitals_alpha.cube." << std::endl;
        log(logStream, outputStream);
        if (showProgress)
        {
            std::cout << std::endl;
        }

        std::cout << "Saving pseudo orbitals in cube format for beta spin..." << std::endl;
        createCube(pseudoOrbitals, domain, outputPrefix + "_lrf_pseudoOrbitals_beta.cube", CubeType::ORBITALS, showProgress, ELFMethod::UNKNOWN, pseudoOrbitalsIndexes, pseudoOrbitalsSpinTypes_beta);
        logStream << "Pseudo orbitals with beta spin saved to " << outputPrefix << "_lrf_pseudoOrbitals_beta.cube." << std::endl;
        log(logStream, outputStream);
        if (showProgress)
        {
            std::cout << std::endl;
        }
    }

    return pseudoOrbitals;
}

void Job::createCube(Orbitals& orbitals, const Domain& domain, const std::string& outputCubeFileName, CubeType cubeType, bool showProgress, const ELFMethod elfMethod, std::vector<int> nums, std::vector<SpinType> typesSpin)
{
    Grid g;
    if (cubeType == CubeType::DENSITY)
    {
        std::cout << "Creating density grid, please wait..." << std::endl;
        g = orbitals.makeGrid(domain, showProgress);
    }
    else if (cubeType == CubeType::ORBITALS)
    {
        g = orbitals.makeOrbGrid(domain, nums, typesSpin, showProgress);
    }
    else if (cubeType == CubeType::ELF)
    {
        g = orbitals.makeELFgrid(domain, elfMethod);
    }
    else
    {
        std::stringstream errorMessage;
        errorMessage << "Error in Job::createCube(): invalid cube type." << std::endl;
        errorMessage << "Please check the documentation and the \"CubeType\" parameter value in the provided input file.";

        print_error(errorMessage.str());

        std::exit(1);
    }

    std::cout << "Writing cube file, please wait..." << std::endl;
    std::ofstream outputFile(outputCubeFileName);
    g.save(outputFile, showProgress);
    outputFile.close();

    std::cout << "Density cube saved to " << outputCubeFileName << '.' << std::endl;
}

void Job::openInputFile()
{
    _inputFile.open(_inputFileName);
    if (_inputFile.fail())
    {
        std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        std::cerr << "Sorry, I cannot open the input file : " << _inputFileName << std::endl;
        std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;

        std::exit(1);
    }
}

void Job::printCriticalPoints()
{
    std::cerr << "Function Job::printCriticalPoints() not implemented yet." << std::endl;
}

void Job::setAllOrbitals(std::vector<int>& orbnums, std::vector<SpinType>& orbspin, Orbitals& o, SpinType spinType, int N)
{
    orbnums.resize(N);
    for(int i = 0; i < N; ++i)
    {
        orbnums[i] = i;
    }

    if (spinType == SpinType::ALPHA)
    {
        orbspin.resize(N, SpinType::ALPHA);
    }
    else if (!o.get_alphaAndBeta() || spinType == SpinType::ALPHA_BETA)
    {
        //orbspin.resize(N, 0);
        orbspin.resize(2 * N, SpinType::ALPHA);
        orbnums.resize(2 * N, 1);
        for(int i = N; i < 2 * N; ++i)
        {
            orbnums[i] = i - N;
            orbspin[i] = SpinType::BETA;
        }
    }
    else
    {
        orbspin.resize(N, SpinType::BETA);
    }
}

void Job::setCustomOrbitals(std::vector<int>& orbnums, std::vector<SpinType>& orbspin, const std::vector<SpinType>& spinList)
{
    for (size_t i = 0; i < spinList.size(); ++i)
    {
        if (spinList[i] == SpinType::ALPHA)
        {
            orbspin.push_back(SpinType::ALPHA);
        }
        else if (spinList[i] == SpinType::BETA)
        {
            orbspin.push_back(SpinType::BETA);
        }
    }

    for (size_t i = 0; i < orbnums.size(); ++i)
    {
        orbnums[i] -= 1;
    }

    if (orbspin.size() < orbnums.size())
    {
        SpinType last = orbspin.back();
        for (size_t i = orbspin.size(); i < orbnums.size(); ++i)
        {
            orbspin.push_back(last);
        }
    }
}

void Job::setOccupiedOrbitals(std::vector<int>& orbnums, std::vector<SpinType>& orbspin, Orbitals& o, SpinType spinType, int N)
{
    std::vector<std::vector<double>> occ=o.get_occupationNumber();
    if (!o.get_alphaAndBeta() and spinType == SpinType::ALPHA) 
    {
        int k = 0;
        for(int i = 0; i < N; ++i)
        {
            if (occ[0][i] > 1e-10)
            {
                orbnums[k] = i;
                k++;
            }
        }
        orbnums.resize(k);
        orbspin.resize(k, SpinType::ALPHA);
    }
    else if (!o.get_alphaAndBeta() and spinType == SpinType::ALPHA_BETA)
    {
        int k = 0;
        orbnums.resize(2 * N, 0);
        for(int i = 0; i < N; ++i)
        {
            if (occ[0][i] > 1e-10)
            {
                orbnums[k] = i;
                k++;
            }
        }
        orbspin.resize(k, SpinType::ALPHA);
        int j = k;
        for(int i = 0; i < N; ++i)
        {
            if (occ[1][i] > 1e-10)
            {
                orbnums[j] = i;
                j++;
            }
        }
        orbnums.resize(j);
        orbspin.resize(j, SpinType::BETA);    
    }
    else
    {
        int j = 0;
        for(int i = 0; i < N; ++i)
        {
            if (occ[1][i] > 1e-10)
            {
                orbnums[j] = i;
                j++;
            }
        }
        orbnums.resize(j);
        orbspin.resize(j, SpinType::BETA);
    }
}

void Job::setOrbitals(Orbitals& o, const int numberOfOrbitals, std::vector<int>& orbitalsNumbers, std::vector<SpinType>& orbitalsSpins, const OrbitalType orbitalType, SpinType spinType, const std::vector<SpinType>& spinList)
{
    switch (orbitalType)
    {
        case OrbitalType::ALL:
        {
            setAllOrbitals(orbitalsNumbers, orbitalsSpins, o, spinType, numberOfOrbitals);
            break;
        }
        case OrbitalType::CUSTOM:
        {
            setCustomOrbitals(orbitalsNumbers, orbitalsSpins, spinList);
            break;
        }
        case OrbitalType::HOMO:
        {
            setHomo(orbitalsNumbers, orbitalsSpins, o, spinType);
            break;
        }
        case OrbitalType::HOMO_LUMO:
        {
            setHomoLumo(orbitalsNumbers, orbitalsSpins, o, spinType);
            break;
        }
        case OrbitalType::LUMO:
        {
            setLumo(orbitalsNumbers, orbitalsSpins, o, spinType);
            break;
        }
        case OrbitalType::OCCUPIED:
        {
            setOccupiedOrbitals(orbitalsNumbers, orbitalsSpins, o, spinType, numberOfOrbitals);
            break;
        }
        case OrbitalType::VIRTUAL:
        {
            setVirtualOrbitals(orbitalsNumbers, orbitalsSpins, o, spinType, numberOfOrbitals);
            break;
        }
        default:
        {
            // should not happen: verified before function called
            std::cerr << "Error: Unknown orbital type encountered in setOrbitals(). This should not happen as the orbital type is verified before this function is called." << std::endl;
        }
    }
}

void Job::setHomo(std::vector<int>& orbnums, std::vector<SpinType>& orbspin, Orbitals& o, SpinType spinType)
{
    int i = 0;
    std::vector<std::vector<double>> occ = o.get_occupationNumber();
    if (!o.get_alphaAndBeta() and spinType == SpinType::ALPHA)
    {
        while (occ[0][i] > 1e-10)
        {
            i++;
        }
        orbnums = {i - 1};
        orbspin = {SpinType::ALPHA};
    }
    else if (!o.get_alphaAndBeta() and spinType == SpinType::ALPHA_BETA)
    {
        int j = 0;
        while (occ[0][i] > 1e-10)
        {
            i++;
        }
        while (occ[1][j] > 1e-10)
        {
            j++;
        }
        orbnums = {i - 1, j - 1};
        orbspin = {SpinType::ALPHA, SpinType::BETA};
    }
    else
    {
        int j = 0;
        while (occ[1][j] > 1e-10)
        {
            j++;
        }
        orbnums = {j - 1};
        orbspin = {SpinType::BETA};
    }
}

void Job::setHomoLumo(std::vector<int>& orbnums, std::vector<SpinType>& orbspin, Orbitals& o, SpinType spinType)
{
    int i = 0;
    std::vector<std::vector<double>> occ = o.get_occupationNumber();
    if (!o.get_alphaAndBeta() and spinType == SpinType::BETA)
    {
        while (occ[0][i] > 1e-10)
        {
            i++;
        }
        orbnums = {i - 1, i};
        orbspin = {SpinType::ALPHA, SpinType::ALPHA};
    }
    else if (!o.get_alphaAndBeta() and spinType == SpinType::ALPHA_BETA)
    {
        int j = 0;
        while (occ[0][i] > 1e-10)
        {
            i++;
        }
        while (occ[1][j] > 1e-10)
        {
            j++;
        }
        orbnums = {i - 1, j - 1, i, j};
        orbspin = {SpinType::ALPHA, SpinType::BETA, SpinType::ALPHA, SpinType::BETA};
    }
    else
    {
        int j = 0;
        while (occ[1][j] > 1e-10)
        {
            j++;
        }
        orbnums = {j - 1, j};
        orbspin = {SpinType::BETA, SpinType::BETA};
    }
}

void Job::setLumo(std::vector<int>& orbnums, std::vector<SpinType>& orbspin, Orbitals& o, SpinType spinType)
{
    int i = 0;
    std::vector<std::vector<double>> occ = o.get_occupationNumber();
    if (!o.get_alphaAndBeta() and spinType == SpinType::ALPHA)
    {
        while (occ[0][i] > 1e-10)
        {
            i++;
        }
        orbnums = {i};
        orbspin = {SpinType::ALPHA};
    }
    else if (!o.get_alphaAndBeta() and spinType == SpinType::ALPHA_BETA)
    {
        int j = 0;
        while (occ[0][i] > 1e-10)
        {
            i++;
        }
        while (occ[1][j] > 1e-10)
        {
            j++;
        }
        orbnums = {i, j};
        orbspin = {SpinType::ALPHA, SpinType::BETA};
    }
    else
    {
        int j = 0;
        while (occ[1][j] > 1e-10)
        {
            j++;
        }
        orbnums = {j};
        orbspin = {SpinType::BETA};
    }
}

void Job::setVirtualOrbitals(std::vector<int>& orbnums, std::vector<SpinType>& orbspin, Orbitals& o, SpinType spinType, int N)
{
    std::vector<std::vector<double>> occ=o.get_occupationNumber();
    if (!o.get_alphaAndBeta() and spinType == SpinType::ALPHA)
    {
        int k = 0;
        for(int i = 0; i < N; ++i)
        {
            if (occ[0][i] < 1e-10)
            {
                orbnums[k] = i;
                k++;
            }
        }
        orbnums.resize(k);
        orbspin.resize(k, SpinType::ALPHA);
    }
    else if (!o.get_alphaAndBeta() and spinType == SpinType::ALPHA_BETA)
    {
        int k = 0;
        orbnums.resize(2 * N, 0);
        for(int i = 0; i < N; ++i)
        {
            if (occ[0][i] < 1e-10)
            {
                orbnums[k] = i;
                k++;
            }
        }
        orbspin.resize(k, SpinType::ALPHA);
        int j = k;
        for(int i = 0; i < N; ++i)
        {
            if (occ[1][i] < 1e-10)
            {
                orbnums[j] = i;
                j++;
            }
        }
        orbnums.resize(j);
        orbspin.resize(j, SpinType::BETA);    
    }
    else
    {
        int j = 0;
        for(int i = 0; i < N; ++i)
        {
            if (occ[1][i] < 1e-10)
            {
                orbnums[j] = i;
                j++;
            }
        }
        orbnums.resize(j);
        orbspin.resize(j, SpinType::BETA);
    }
}
