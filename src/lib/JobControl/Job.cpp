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

Job::Job(std::string inputFileName):
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

    return read;
}

bool Job::readCharges(std::vector<double>& charges)
{
    bool read = readListType<double>(_inputFile, "Charges", charges);

    if (!read)
    {
        std::cout << "Note: the \"Charges\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (Charges=1)." << std::endl << std::endl;

        charges = { 1.0 };
    }

    return read;
}

bool Job::readCutoff(double& cutoff)
{
    bool read = readOneType<double>(_inputFile, "Cutoff", cutoff);

    if (!read)
    {
        std::cout << "Note: the \"Cutoff\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (Cutoff=0)." << std::endl << std::endl;

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
        std::cout << "The program will use the default value (ELFMethod=Savin)." << std::endl << std::endl;
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
        std::cout << "The program will consider all excited states available in the excited states file (MaxNumberOfExcitedStates=-1)." << std::endl << std::endl;

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
        std::cout << "The program will use the default value (NuclearCutoff=0.0529177210544)." << std::endl << std::endl;

        nuclearCutoff = 0.1;
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
        std::cout << "The program will use the default value (OutputPrefix=\"output\")." << std::endl << std::endl;

        outputPrefix = "output";
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
        std::cout << "The program will use the default value (OrbitalType=All)." << std::endl << std::endl;
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
        std::cout << "The program will use the default value (PartitionMethod=On-Grid)." << std::endl << std::endl;
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
    std::vector<double> positionValues;
    bool read = readListType<double>(_inputFile, "Positions", positionValues);

    if (read && !positionValues.empty())
    {
        if (positionValues.size() % 3 != 0)
        {
            std::stringstream errorMessage;
            errorMessage << "Error: incorrect number of values for the \"Positions\" parameter (multiple of three values expected)." << std::endl;
            errorMessage << "Please check documentation and the \"Positions\" parameter values in " << _inputFileName << '.';

            print_error(errorMessage.str());

            std::exit(1);
        }

        for (size_t i = 0; i < positionValues.size(); i += 3)
        {
            positions.push_back({ positionValues[i] * Constants::ANGSTROM_TO_BOHR_RADIUS,
                                  positionValues[i + 1] * Constants::ANGSTROM_TO_BOHR_RADIUS,
                                  positionValues[i + 2] * Constants::ANGSTROM_TO_BOHR_RADIUS });
        }
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
        std::cout << "The program will use the default value (RunType=Help)." << std::endl << std::endl;
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
    bool read = readOneString(_inputFile, "SavePseudoOrbitals", strSavePseudoOrbitals) && to_lower(strSavePseudoOrbitals) == "true";
    
    savePseudoOrbitals = false;
    if (!read)
    {
        std::cout << "Note: the \"SavePseudoOrbitals\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (SavePseudoOrbitals=False)." << std::endl << std::endl;
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

bool Job::readShowProgress(bool& showProgress)
{
    std::string strShowProgress;
    bool read = readOneString(_inputFile, "ShowProgress", strShowProgress);

    showProgress = false;
    if (!read)
    {
        std::cout << "Note: the \"ShowProgress\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (ShowProgress=False)." << std::endl << std::endl;
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

bool Job::readSize(GridSize& gridSize, CustomSizeData& customSizeData)
{
    std::string strGridSize;
    bool read = readOneString(_inputFile, "Size", strGridSize);

    if (!read)
    {
        std::cout << "Note: the \"Size\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (Size=Medium)." << std::endl << std::endl;

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
        std::cout << "The program will use the default value (SpinType=Alpha-Beta)." << std::endl << std::endl;

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
// JOBS
//----------------------------------------------------------------------------------------------------//

void Job::run_computeCondensedLinearResponse()
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

    Descriptors descriptors_FD = computeDescriptorsFD(analyticFilesNames[1], analyticFilesNames[2], analyticFilesNames[3], beckeParams[0], beckeParams[1], beckeParams[2]);

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

void Job::run_computeDescriptors()
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
            Descriptors descriptors = computeDescriptorsFD(analyticFilesNames[0], analyticFilesNames[1], analyticFilesNames[2]);
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
            Descriptors D = computeDescriptors(gridFilesNames[0], gridFilesNames[1], gridFilesNames[2], energies[0], energies[1], partitionMethod);
            std::cout << D;
        }
        else if (energies.size() == 3)
        {
            std::cout << " Reading Total Energies: E1 = " << energies[0] << ", E2 = " << energies[1] << " and E3 = " << energies[2] << std::endl;
            Descriptors D = computeDescriptors(gridFilesNames[0], gridFilesNames[1], gridFilesNames[2], energies, partitionMethod);
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

void Job::run_computeEnergyWithPointCharges()
{
    
}

void Job::run_computeGridDifference()
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

    
    // Compute grid difference
    std::cout << "Computing difference between " << gridFilesNames[0] << " and " << gridFilesNames[1] << '.' << std::endl;
    
    computeGridDifference(gridFilesNames[0], gridFilesNames[1], gridFilesNames[2]);
    
    std::cout << "Difference grid saved to file " << gridFilesNames[2] << '.' << std::endl;
}

void Job::run_computeIntegrals()
{
    // Read grid files names
    std::vector<std::string> gridFilesNames;
    readGridFilesNames(gridFilesNames);


    // Read partition method
    PartitionMethod partitionMethod;
    readPartitionMethod(partitionMethod);

    // Check partition method validity
    if (partitionMethod == PartitionMethod::BECKE || partitionMethod == PartitionMethod::FD
                                                  || partitionMethod == PartitionMethod::FMO)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: partitionMethod \"" << to_string(partitionMethod) << "\" invalid for this job." << std::endl;
        errorMessage << "Please check documentation and \"PartitionMethod\" parameter value in " << _inputFileName << '.';

        print_error(errorMessage.str());

        std::exit(1);
    }

    // Print partition method
    std::cout << "Volume partition method: " << to_string(partitionMethod) << std::endl << std::endl;

    
    // Read cutoff
    double cutoff = 0.0;
    readCutoff(cutoff);


    // Build basins
    std::cout << "Reading file " << gridFilesNames[0] << " to build basins." << std::endl << std::endl;

    GridCP gcp;
    if (partitionMethod == PartitionMethod::BBS)
    {
        buildBasinsBySign(gcp, gridFilesNames[0], cutoff, false);
    }
    else if (partitionMethod == PartitionMethod::B2S)
    {
        buildBasinsBySign(gcp, gridFilesNames[0], cutoff, true);
    }
    else
    {
        buildBasins(gcp, gridFilesNames[0], partitionMethod);
    }

    
    // Compute local integrals
    computeLocalIntegrals(gcp, gridFilesNames);
}

void Job::run_computePartialCharges()
{
    // Read grid files names
    std::vector<std::string> gridFilesNames;
    readGridFilesNames(gridFilesNames);


    // Read partition method
    PartitionMethod partitionMethod;
    readPartitionMethod(partitionMethod);


    // Check partition method validity
    if (partitionMethod == PartitionMethod::BBS || partitionMethod == PartitionMethod::B2S)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: partitionMethod \""  << to_string(partitionMethod) << "\" invalid for this job." << std::endl;
        errorMessage << "Please check documentation and the \"PartitionMethod\" parameter value in " << _inputFileName << '.';

        print_error(errorMessage.str());

        std::exit(1);
    }


    // Print partition method
    std::cout << "Volume partition method: " << to_string(partitionMethod) << std::endl <<std::endl;

    
    // Compute partial charges
    std::vector<double> charges = computePartialCharges(gridFilesNames[0], partitionMethod);
}

void Job::run_convertOrbitals()
{
    // Read analytic files names
    std::vector<std::string> analyticFilesNames;
    readAnalyticFilesNames(analyticFilesNames);


    // Check number of analytic files names
    if (analyticFilesNames.size() != 2)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect number of analytic files names (two files expected)." << std::endl;
        errorMessage << "Please check the documentation and the number of files specified in the \"AnalyticFiles\" parameter in " << _inputFileName << '.';

        print_error(errorMessage.str());

        std::exit(1);
    }


    // Checking if the files have a different format
    std::regex fileExtensionRegex("\\.([a-zA-Z0-9]+)$");
    std::smatch matchInputFileExtension;
    std::smatch matchOutputFileExtension;
    std::string inputFileExtension;
    std::string outputFileExtension;
    std::regex_search(analyticFilesNames[0], matchInputFileExtension, fileExtensionRegex);
    std::regex_search(analyticFilesNames[1], matchOutputFileExtension, fileExtensionRegex);

    if (matchInputFileExtension[0] != matchOutputFileExtension[0])
    {
        // Loading orbitals and saving in the new format
        Orbitals o;
        computeOrbitalsOrBecke<Orbitals>(o, analyticFilesNames[0]);
        o.Save(analyticFilesNames[1]);
    }
    else
    {
        std::cout << "Input and output files have the same format (" << inputFileExtension << "). Nothing to be done." << std::endl;
    }
}

void Job::run_lambdaDiagnostic()
{
    //Read analytic file name
    std::vector<std::string> analyticFilesNames;
    readAnalyticFilesNames(analyticFilesNames);


    // Check number of analytic files names
    if (analyticFilesNames.size() != 1)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect number of analytic files names (one file expected)." << std::endl;
        errorMessage << "Please check the documentation and the number of files specified in the \"AnalyticFiles\" parameter in " << _inputFileName << '.';

        print_error(errorMessage.str());

        std::exit(1);
    }


    // Read grid size
    GridSize gridSize;
    CustomSizeData customSizeData;
    readSize(gridSize, customSizeData);


    // Read transitions file
    std::string transitionsFileName;
    readTransitionsFileName(transitionsFileName);


    // Loading orbitals
    Orbitals orbitals;
    computeOrbitalsOrBecke<Orbitals>(orbitals, analyticFilesNames[0]);


    // Setting orbitals
    std::vector<SpinType> orbitalsSpins;
    std::vector<int> orbitalsNumbers;
    setOrbitals(orbitals, orbitals.get_numberOfMo(), orbitalsNumbers, orbitalsSpins, OrbitalType::ALL, SpinType::ALPHA_BETA);


    // Building domain and grid
    std::cout << "Building domain and grid, please wait..." << std::endl;

    Domain domain = buildDomainForCube(orbitals, gridSize, customSizeData, orbitalsNumbers.size());
    Grid orbitalsGrid = orbitals.makeOrbGrid(domain, orbitalsNumbers, orbitalsSpins);


    // Reading transitions file
    std::vector<ExcitedState> excitedStates;
    ExcitedState::readTransitions(transitionsFileName, excitedStates, orbitals.get_energy());
    std::cout << "Number of excited states read: " << excitedStates.size() << std::endl << std::endl;


    // Printing lambda diagnostic for each excited state
    for (const ExcitedState& excitedState : excitedStates)
    {
        std::cout << excitedState << std::endl;
        excitedState.printLambdaDiagnostic(orbitalsGrid);
        std::cout << "-----------------" << std::endl << std::endl;
    }
}

void Job::run_makeDensityCube()
{
    // Read analytic file
    std::vector<std::string> analyticFilesNames;
    readAnalyticFilesNames(analyticFilesNames);


    // Check number of analytic files names
    if (analyticFilesNames.size() != 1)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect number of analytic files names (one file expected)." << std::endl;
        errorMessage << "Please check the documentation and the number of files specified in the \"AnalyticFiles\" parameter in " << _inputFileName << '.';

        print_error(errorMessage.str());

        std::exit(1);
    }

    
    // Loading orbitals
    Orbitals o;
    computeOrbitalsOrBecke<Orbitals>(o, analyticFilesNames[0]);
    

    // Read size
    GridSize gridSize;
    CustomSizeData customSizeData;
    readSize(gridSize, customSizeData);


    // Read grid file name
    std::vector<std::string> gridFilesName;
    readGridFilesNames(gridFilesName);


    // Check number of grid files names
    if (gridFilesName.size() != 1)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect number of grid files names (one file expected)." << std::endl;
        errorMessage << "Please check the documentation and the number of files specified in the \"GridFiles\" parameter in " << _inputFileName << '.';

        print_error(errorMessage.str());

        std::exit(1);
    }


    // Building domain
    std::cout << "Building domain, please wait..." << std::endl;

    Domain d = buildDomainForCube(o, gridSize, customSizeData, 1);


    // Creating density cube
    std::cout << "Creating density cube, please wait..." << std::endl;

    createCube(o, d, gridFilesName[0], 0);

    std::cout << "Density cube saved to file " << gridFilesName[0] << '.' << std::endl;
}

void Job::run_makeOrbitalsCube()
{
    // Read progress bar display option
    bool showProgress;
    readShowProgress(showProgress);


    // Read analytic file names
    std::vector<std::string> analyticFilesNames;
    readAnalyticFilesNames(analyticFilesNames);


    // Check number of analytic files names
    if (analyticFilesNames.size() != 1)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect number of analytic files names (one file expected)." << std::endl;
        errorMessage << "Please check the documentation and the number of files specified in the \"AnalyticFiles\" parameter in " << _inputFileName << '.';

        print_error(errorMessage.str());

        std::exit(1);
    }


    // Read grid files name
    std::vector<std::string> gridFileName;
    readGridFilesNames(gridFileName);


    // Check number of grid files names
    if (gridFileName.size() != 1)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect number of grid files names (one file expected)." << std::endl;
        errorMessage << "Please check the documentation and the number of files specified in the \"GridFiles\" parameter in " << _inputFileName << '.';

        print_error(errorMessage.str());

        std::exit(1);
    }


    // Read size
    GridSize gridSize;
    CustomSizeData customSizeData;
    readSize(gridSize, customSizeData);

    
    // Read spin type
    SpinType spinType;
    readSpinType(spinType);


    // Read orbitals type
    OrbitalType orbitalType;
    readOrbitalType(orbitalType);


    // Loading orbitals
    Orbitals o;
    computeOrbitalsOrBecke<Orbitals>(o, analyticFilesNames[0]);


    // Setting orbitals
    int numberOfOrbitals = o.get_numberOfMo();
    std::vector<SpinType> orbitalsSpins;
    std::vector<int> orbitalsNumbers;
    std::vector<SpinType> spinList;
    
    if (orbitalType == OrbitalType::CUSTOM)
    {
        // Read custom orbitals list
        readOrbitalsNumbers(orbitalsNumbers);
        readSpinList(spinList);
        

        // Check if the sizes of both lists match
        if (spinList.size() != orbitalsNumbers.size())
        {
            std::stringstream errorMessage;
            errorMessage << "Error: sizes of orbitals numbers list and spins list do not match." << std::endl;
            errorMessage << "Please check the documentation and the number of items specified in the \"OrbitalsList\" and \"SpinList\" parameters in " << _inputFileName << '.';

            print_error(errorMessage.str());

            std::exit(1);
        }
    }

    setOrbitals(o, numberOfOrbitals, orbitalsNumbers, orbitalsSpins, orbitalType, spinType, spinList);


    // Building domain
    std::cout << "Building domain, please wait..." << std::endl;

    Domain d = buildDomainForCube(o, gridSize, customSizeData, orbitalsNumbers.size());


    // Creating orbitals cube
    std::cout << "Creating orbitals cube, please wait..." << std::endl;

    createCube(o, d, gridFileName[0], 1, showProgress, ELFMethod::UNKNOWN, orbitalsNumbers, orbitalsSpins);

    std::cout << "Data saved to file: " << gridFileName[0] << std::endl;
}

void Job::run_makeELFCube()
{
    // Read progress bar display option
    bool showProgress;
    readShowProgress(showProgress);
    

    // Read analytic file
    std::vector<std::string> analyticFilesNames;
    readAnalyticFilesNames(analyticFilesNames);


    // Check number of analytic files names
    if (analyticFilesNames.size() != 1)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect number of analytic files names (one file expected)." << std::endl;
        errorMessage << "Please check the documentation and the number of files specified in the \"AnalyticFiles\" parameter in " << _inputFileName << '.';

        print_error(errorMessage.str());

        std::exit(1);
    }


    // Loading orbitals
    Orbitals o;
    computeOrbitalsOrBecke<Orbitals>(o, analyticFilesNames[0]);


    // Read size
    GridSize gridSize;
    CustomSizeData customSizeData;
    readSize(gridSize, customSizeData);


    // Read grid file name
    std::vector<std::string> gridFilesName;
    readGridFilesNames(gridFilesName);


    // Check number of grid files names
    if (gridFilesName.size() != 1)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: incorrect number of grid files names (one file expected)." << std::endl;
        errorMessage << "Please check the documentation and the number of files specified in the \"GridFiles\" parameter in " << _inputFileName << '.';

        print_error(errorMessage.str());

        std::exit(1);
    }


    // Read ELF method
    ELFMethod elfMethod;
    readELFMethod(elfMethod);


    // Building domain
    std::cout << "Building domain, please wait..." << std::endl;

    Domain d = buildDomainForCube(o, gridSize, customSizeData, 1);

    // Creating ELF cube
    std::cout << "Creating ELF cube, please wait..." << std::endl;

    createCube(o, d, gridFilesName[0], 2, showProgress, elfMethod);

    std::cout << "ELF cube saved to file " << gridFilesName[0] << '.' << std::endl;
}


//----------------------------------------------------------------------------------------------------//
// OTHER PRIVATE METHODS
//----------------------------------------------------------------------------------------------------//

void Job::buildBasins(GridCP& gridCP, const std::string& gridFileName, PartitionMethod partitionMethod)
{
    std::ifstream gridFile(gridFileName);
    if (!gridFile)
    {
        print_error("Error: could not read file " + gridFileName + ".");
        std::exit(1);
    }

    Grid g(gridFile, PeriodicTable());
    gridFile.close();

    gridCP.buildBasins(g, partitionMethod);
}

void Job::buildBasinsBySign(GridCP& gridCP, const std::string& gridFileName,  double cutoff, bool two) 
{
    std::ifstream gridFile(gridFileName);
    if (!gridFile)
    {
        print_error("Error in Job::buildBasinsBySign(): could not read file " + gridFileName + ".");
        std::exit(1);
    }

    Grid g(gridFile, PeriodicTable());
    gridFile.close();

    if (two)
    {
            gridCP.build2BasinSign(g);
    }    
    else
    {
            gridCP.buildBasinsBySign(g, cutoff);
    }    
}

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

Descriptors Job::computeDescriptors(const std::string& gridFileName1, const std::string& gridFileName2, const std::string& gridFileName3, double ionizationEnergy, double electronAffinity, PartitionMethod partitionMethod)
{
    std::ifstream gridFile1(gridFileName1);
    if (!gridFile1)
    {
        print_error("Error in Job::computeDescriptors(): could not read file " + gridFileName1 + ".");
        std::exit(1);
    }

    std::ifstream gridFile2(gridFileName2);
    if (!gridFile2)
    {
        print_error("Error in Job::computeDescriptors(): could not read file " + gridFileName2 + ".");
        std::exit(1);
    }

    std::ifstream gridFile3(gridFileName3);
    if (!gridFile3)
    {
        print_error("Error in Job::computeDescriptors(): could not read file " + gridFileName3 + ".");
        std::exit(1);
    }

    Descriptors D(gridFile1, gridFile2, gridFile3, ionizationEnergy, electronAffinity, partitionMethod);
    gridFile1.close();
    gridFile2.close();
    gridFile3.close();

    return D;
}

Descriptors Job::computeDescriptors(const std::string& gridFileName1, const std::string& gridFileName2, const std::string& gridFileName3, const std::vector<double>& energies, PartitionMethod partitionMethod)
{
    std::ifstream gridFile1(gridFileName1);
    if (!gridFile1)
    {
        print_error("Error in Job::computeDescriptors(): could not read file " + gridFileName1 + ".");
        std::exit(1);
    }

    std::ifstream gridFile2(gridFileName2);
    if (!gridFile2)
    {
        print_error("Error in Job::computeDescriptors(): could not read file " + gridFileName2 + ".");
        std::exit(1);
    }

    std::ifstream gridFile3(gridFileName3);
    if (!gridFile3)
    {
        print_error("Error in Job::computeDescriptors(): could not read file " + gridFileName3 + ".");
        std::exit(1);
    }

    Descriptors D(gridFile1, gridFile2, gridFile3, energies, partitionMethod);
    gridFile1.close();
    gridFile2.close();
    gridFile3.close();

    return D;
}

Descriptors Job::computeDescriptorsFD(const std::string& ANAFileName1, const std::string& ANAFileName2, const std::string& ANAFileName3, int kmax, int lebedev_order, int radial_grid_factor)
{
    std::cout << "Building Structure object... ";
    Structure s = returnStruct(ANAFileName1);
    std::cout << std::endl;

    std::vector<double> E(0);

    std::vector<double> Q0 = computePartialChargesAndEnergy(E, ANAFileName1, kmax, lebedev_order, radial_grid_factor);
    std::vector<double> Qm = computePartialChargesAndEnergy(E, ANAFileName2, kmax, lebedev_order, radial_grid_factor);
    std::vector<double> Qp = computePartialChargesAndEnergy(E, ANAFileName3, kmax, lebedev_order, radial_grid_factor);

    return Descriptors(s, Q0, Qm, Qp, E);
}

void Job::computeGridDifference(const std::string& minuendGridFileName, const std::string& subtrahendGridFileName, const std::string& outputGridFileName)
{
    std::ifstream minuendFile(minuendGridFileName);
    if (!minuendFile)
    {
        print_error("Error in Job::computeGridDifference(): could not read file " + minuendGridFileName + ".");
        std::exit(1);
    }

    Grid minuendGrid(minuendFile, PeriodicTable());
    minuendFile.close();

    std::ifstream subtrahendFile(subtrahendGridFileName);
    if (!subtrahendFile)
    {
        print_error("Error in Job::computeGridDifference(): could not read file " + subtrahendGridFileName + ".");
        std::exit(1);
    }

    Grid subtrahendGrid(subtrahendFile, PeriodicTable());
    subtrahendFile.close();

    Grid diff = minuendGrid - subtrahendGrid;

    std::ofstream outputGridFile(outputGridFileName, std::ios::out);
    if (!outputGridFile)
    {
        print_error("Error in Job::computeGridDifference(): failed to write to file " + outputGridFileName + ".");
        std::exit(1);
    }

    diff.save(outputGridFile);
    outputGridFile.close();

    std::cout << "Grid has been saved to : " << outputGridFileName << std::endl;
}

void Job::computeLocalIntegrals(GridCP& gridCP, const std::vector<std::string>& gridFileNames)
{
    for(size_t i = 0; i < gridFileNames.size(); ++i)
    {
        std::ifstream gridFile(gridFileNames[i]);
        if (!gridFile)
        {
            print_error("Error in Job::computeLocalIntegrals(): could not read file " + gridFileNames[i] + ".");
            std::exit(1);
        }

        Grid g(gridFile, PeriodicTable());
        gridFile.close();

        gridCP.computeIntegrals(g);
        gridCP.printCriticalPoints();
    }
}

template<typename T, typename U>
U Job::computeOrbitalsOrBecke(const std::string& analyticFileName)
{
    std::ifstream analyticFile(analyticFileName);
    if (!analyticFile)
    {
        print_error("Error in Job::computeOrbitalsOrBecke(): could not read file " + analyticFileName + ".");
        std::exit(1);
    }

    T analyticFileParser(analyticFile);
    analyticFile.close();

    U analyticObject(analyticFileParser, Binomial(100), PeriodicTable());
    
    return analyticObject;
}

std::vector<double> Job::computePartialCharges(const std::string& gridFileName, PartitionMethod partitionMethod)
{
    std::vector<double> charges;

    std::ifstream gridFile(gridFileName);
    if (!gridFile)
    {
        print_error("Error in Job::computeLocalIntegrals(): could not read file " + gridFileName + ".");
        std::exit(1);
    }

    Grid grid(gridFile, PeriodicTable());
    gridFile.close();

    if (partitionMethod == PartitionMethod::BECKE)
    {
        Becke B(grid);
        B.partial_charge(grid);
        B.printCharges();

        charges = B.get_partial_charge();
    }
    else
    {
        GridCP gridcp;
        gridcp.buildBasins(grid, partitionMethod);
        gridcp.computeIntegrals(grid);
        gridcp.printCriticalPoints();

        charges = gridcp.computeAIMCharges(grid);
    }

    return charges;
}

std::vector<double> Job::computePartialChargesAndEnergy(std::vector<double>& energies, const std::string& analyticFileName, int kmax, int lebedev_order, int radial_grid_factor)
{
    std::cout << "Building Becke object... ";
    Becke B;
    computeOrbitalsOrBecke<Becke>(B, analyticFileName);

    std::vector<std::vector<double>> partialChargesAndEnergy = B.PartialChargesAndEnergy(kmax, lebedev_order, radial_grid_factor);

    energies.push_back(partialChargesAndEnergy[0][0]);

    return partialChargesAndEnergy[1];
}

template<typename T>
void Job::computeOrbitalsOrBecke(T& analyticObject, const std::string& analyticFileName)
{
    if (analyticFileName.find(".wfx") != std::string::npos)
    {
        std::cout << "Reading data from " << analyticFileName << "... Please wait." << std::endl;
        analyticObject = computeOrbitalsOrBecke<WFX, T>(analyticFileName);
    }
    else if (analyticFileName.find(".fchk") != std::string::npos)
    {
        std::cout << "Reading data from " << analyticFileName << "... Please wait." << std::endl;
        analyticObject = computeOrbitalsOrBecke<FCHK, T>(analyticFileName);
    }
    else if (analyticFileName.find(".molden") != std::string::npos)
    {
        std::cout << "Reading data from " << analyticFileName << "... Please wait." << std::endl;
        analyticObject = computeOrbitalsOrBecke<MOLDENGAB, T>(analyticFileName);
    }
    else if (analyticFileName.find(".gab") != std::string::npos)
    {
        std::cout << "Reading data from " << analyticFileName << "... Please wait." << std::endl;
        analyticObject = computeOrbitalsOrBecke<MOLDENGAB, T>(analyticFileName);
    }
    else if (analyticFileName.find(".log") != std::string::npos)
    {
        std::cout << "Reading data from " << analyticFileName << "... Please wait." << std::endl;
        analyticObject = computeOrbitalsOrBecke<LOG, T>(analyticFileName);
    }
    else
    {
        std::stringstream errorMessage;
        errorMessage << "Error: invalid file format for file \"" << analyticFileName << "\"." << std::endl;
        errorMessage << "Please check the documentation and the \"AnalyticFiles\" parameter value in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }
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
        createCube(pseudoOrbitals, domain, outputPrefix + "_lrf_pseudoOrbitals_alpha.cube", 1, showProgress, ELFMethod::UNKNOWN, pseudoOrbitalsIndexes, pseudoOrbitalsSpinTypes_alpha);
        logStream << "Pseudo orbitals with alpha spin saved to " << outputPrefix << "_lrf_pseudoOrbitals_alpha.cube." << std::endl;
        log(logStream, outputStream);
        if (showProgress)
        {
            std::cout << std::endl;
        }

        std::cout << "Saving pseudo orbitals in cube format for beta spin..." << std::endl;
        createCube(pseudoOrbitals, domain, outputPrefix + "_lrf_pseudoOrbitals_beta.cube", 1, showProgress, ELFMethod::UNKNOWN, pseudoOrbitalsIndexes, pseudoOrbitalsSpinTypes_beta);
        logStream << "Pseudo orbitals with beta spin saved to " << outputPrefix << "_lrf_pseudoOrbitals_beta.cube." << std::endl;
        log(logStream, outputStream);
        if (showProgress)
        {
            std::cout << std::endl;
        }
    }

    return pseudoOrbitals;
}





// TypeFlag specifies the type of grid you wnat to make.
// For now there are 3 types available. electronic density, ELF and orbitals. Others can be added in the else ifs. additional parameters shoud be added before the default values.
void Job::createCube(Orbitals& orbitals, const Domain& domain, const std::string& cubeFileName, int TypeFlag, bool showProgress, const ELFMethod elfMethod, std::vector<int> nums, std::vector<SpinType> typesSpin)
{
    Grid g;
    if (TypeFlag == 0) // Electronic density
    {
        g=orbitals.makeGrid(domain);

    }
    else if (TypeFlag == 1) // Orbitals
    {
        g = orbitals.makeOrbGrid(domain, nums, typesSpin, showProgress);
    }
    else // ELF
    {
        if (elfMethod == ELFMethod::BECKE)
        {
            g = orbitals.makeELFgrid(domain, 0);
        }
        else // SAVIN
        {
            g=orbitals.makeELFgrid(domain);
        }

    }

    std::cout << "Writing cube file... Please wait." << std::endl;
    std::ofstream out(cubeFileName);
    g.save(out, showProgress);
    out.close();
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

Structure Job::returnStruct(const std::string& analyticFileName)
{
    Becke B;
    computeOrbitalsOrBecke<Becke>(B, analyticFileName);

    return B.get_molecule();
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
            orbspin[i] = SpinType::BETA; // TO BE TESTED
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


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

void Job::run() 
{
    // Determine run type
    RunType runType;
    readRunType(runType);


    // Print current job;
    print_title("Current job: " + to_string(runType));


    // Execute job
    switch (runType)
    {
        case RunType::COMPUTE_CONDENSED_LINEAR_RESPONSE:
        {
            run_computeCondensedLinearResponse();
            break;
        }
        case RunType::COMPUTE_DESCRIPTORS:
        {
            run_computeDescriptors();
            break;
        }
        case RunType::COMPUTE_GRID_DIFFERENCE:
        {
            run_computeGridDifference();
            break;
        }
        case RunType::COMPUTE_INTEGRALS:
        {
            run_computeIntegrals();
            break;
        }
        case RunType::COMPUTE_PARTIAL_CHARGES:
        {
            run_computePartialCharges();
            break;
        }
        case RunType::CONVERT_ORBITALS:
        {
            run_convertOrbitals();
            break;
        }
        case RunType::LAMBDA_DIAGNOSTIC:
        {
            run_lambdaDiagnostic();
            break;
        }
        case RunType::MAKE_DENSITY_CUBE:
        {
            run_makeDensityCube();
            break;
        }
        case RunType::MAKE_ORBITALS_CUBE:
        {
            run_makeOrbitalsCube();
            break;
        }
        case RunType::MAKE_ELF_CUBE:
        {
            run_makeELFCube();
            break;
        }
        default:
        {
            // should not happen: verified before switch
            std::cerr << "Error: Unknown run type encountered in Job::run(). This should not happen as the run type is verified before." << std::endl;
        }
    }
}


