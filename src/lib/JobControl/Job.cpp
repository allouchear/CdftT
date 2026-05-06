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
    _inputFileName("input.txt")
{
    setJobList();
    openInputFile();
}

Job::Job(std::string inputFileName):
    _inputFileName(inputFileName)
{
    setJobList();
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

void Job::printListOfRunTypes()
{
    std::cout << "Available jobs (runType=) :" << std::endl << std::endl;

    for(size_t i=0;i<_jobsList.size();i++)
    {
        std::cout << "--------------------------------------------------" << std::endl;
        std::cout << _jobsList[i] << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;
        std::cout << _jobDescription[i] << std::endl << std::endl;
    }
}

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
    printResultsEnergyWithPointCharges(states, ionicMatrixes, chargeNucleiContributions, charges, chargesPositions, loopOnAtoms, atoms, outputPrefix, outputStream, verbose);
    

    /****************/
    /* REGULAR GRID */
    /****************/
    /*
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
        printResultsEnergyWithPointCharges(states, ionicMatrixes_regularGrid, chargeNucleiContributions, charges, chargesPositions, loopOnAtoms, atoms, outputPrefix + "_regularGrid", outputStream, verbose);
    }
    */

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
        printResultsEnergyWithPointCharges(states, ionicMatrixes_becke, chargeNucleiContributions, charges, chargesPositions, loopOnAtoms, atoms, outputPrefix + "_becke", outputStream, verbose);
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

void Job::run_computeLinearResponseWithPointCharges()
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
    printResultsLinearResponseWithPointCharges(eigenvalues, ionicPotentialVectors, charges, chargesPositions, loopOnAtoms, atoms, outputStream, verbose);
    

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
        printResultsLinearResponseWithPointCharges(eigenvalues_becke, ionicPotentialVectors_becke, charges, chargesPositions, loopOnAtoms, atoms, outputStream, verbose);



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

void Job::run_help()
{
    printListOfRunTypes();
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

void Job::setJobList()
{
    _jobsList = { "Help",
                  "ComputeDescriptors",
                  "ComputeEnergyWithPointCharges",
                  "ComputeGridDifference",
                  "ComputeIntegrals",
                  "ComputePartialCharges",
                  "ConvertOrbitals",
                  "LambdaDiagnostic",
                  "LinearResponse",
                  "MakeDensityCube",
                  "MakeELFCube",
                  "MakeOrbitalsCube" };
    
    _jobDescription = { "Details are given for the available jobs run by this program.\nExample input files for each job are also given. In this format, comment lines are specified by # at the start of the line",
                        "Computation of chemical descriptors from analytic or cube files using on-grid, near-grid, near-grid-refinement and Becke. Frontier Molecular Orbitals(FMO) and finite difference(FD) are methods also provided for the computation. FMO requires 1 analytic file (.log, .wfx, .molden,...). FD requires 3 analytic files. The other methods require cube files of nucleophilic, electrophilic and radical attacks for the molecule. Energies must also be given by the user:\nif two are given, they are assumed to be the ionization energy and the electron affinity. If 3 are given they are assumed to be the total energies of each file. \n\n Example format for input file :\n\n#RunType=Help\n#RunType=ComputeDescriptorsFromCubes\n#GridFileName\nGrids=grid1.cube, grid2.cube, grid3.cube\nPartitionMethod=on-grid\nEnergies=I, A or E1,E2,E3",
                        "Computes the new energy levels of a system when one or many point charges are added.",
                        "Computes the differences of values of the first two grids provides and assigns them to the third.\n\n Example format for input file : \n\n#Runtype=Help\nRunType=ComputeDifference\n#GridFileName\nGrids=in1.cube, in2.cube, out.cube ",
                        "Computes local integrals of grids on volumes defined by method of choice. A grid is required to define the volumes.\nThe additional grids provided by the user should contain the quantities to be integrated.\n\n **on-grid** : to define volumes using on-grid AIM. Requires electronic density grid.\n **near-grid** : to define volumes using near-grid AIM. Requires electronic density grid.\n **near-grid-refinement : to define volumes using near-grid-refinement AIM. Requires electronic density grid.\n **VDD** : to define volumes by distance to atoms. Can use any type of density.\n **BBS** : Build Basins By SIGN. Requires a grid of density difference. A job is provided in the program to obtain such a grid. An additional input *Cutoff=* is required for BBS that sets a threshold for insignificant values.\n **B2S** : Build 2 basins by SIGN. Same as BBS but only constructs two volumes.\n\n Example format for input file :\n\n#RunType=Help\n#RunType=ComputeIntegrals\n#GridFileName\nGrids=gridDefiningVolumes.cube, grid1ToBeIntegrated.cube, grid2ToBeIntegrated.cube\nPartitionMethod=BBS\nCutoff=1e-10",
                        "Grid-based computations of partial charges of the molecule. We provide 5 ways of computing atomic volumes, the first 3 of which are based on Bader's Atoms in molecule.\n\n **on-grid** : follows Tang's algorithm to find Bader volumes.\n **near-grid** : more precise version of on-grid.\n **near-grid-refinement** : even more precise. Requires more time.\n **VDD** topological method : assigns points to volumes by distance to closest atom.\n **Becke** : uses a regular density grid to interpolate Becke's atomic variable grids.\n\n Example format for input file :\n\n#RunType\n#RunType=Help\nRunType=ComputePartialCharges\n#GridFileName\nGrids=h2o_80_0.gcube \nPartitionMethod=on-grid\n\nW. Tang, E. Sanville, G. Henkelman, A grid-based bader analysis algorithm without lattice bias, Journal of Physics: Condensed Matter 21 (8) (2009) 084204.",
                        "Convert Analytical file.\nSupported file formats are : wfx, fchk, log, molden, gab.\nOutput supported : wfx, molden, gab\n\n Example format for input file : \n\n#RunType=Help\nRunType=ConvertOrbitals\nAnalyticFiles=input.wfx, output.molden",
                        "Prints the result of the Lambda diagnostic test, as described by Peach et al., that judges the reliability of TDDFT excited states calculations. It also allows to validate the grid size configuration by computing overlap integrals between the orbitals involved in the excited states.",
                        "TO BE COMPLETED!",
                        "Create a density grid and save it in .cube format. .wfx , .fchk , .molden , .gab and .log are supported as input files.\nthe user can choose from 3 standard grid sizes:\ncoarse ( 3 pts / Bohr)\nMedium (6 pts / Bohr)\nFine (12 pts / Bohr)\n\nA custom size is also provided in which the user enters the domain data as follows:\nNx, Ny, Nz, Ox, Oy, Oz, T11, T12, T13, T21, T22, T23, T31, T32, T33\nWhere N is the number of points in the ith direction, Oi are the coordinates of the bottom left corner of the cube and Tij are the coeficients of the translation vector.\n\n Example format for input file : \n\n#RunType=Help\nRunType=MakeDensityCube\n#GridFileName\nAnalyticFile=filename.wfx\nSize=Custom\nCustomSizeData=80,80,80,5,5,5,0.15,0,0,0,0.15,0,0,0,0.15\nGrid=save.cube ",
                        "Create a grid and compute the Electron Localisation Function (ELF) using either Savin or Becke method. Grid domain is defined the same as the MakeDensityCube.\nBy default the program will run Savin ELF.\n\n Example format for input file : \n\n#RunType=Help\nRunType=MakeELFCube\n#GridFileName\nAnalyticFile=filename.wfx\nSize=Medium\nELFmethod=Becke\nGrid=save.cube",
                        "Compute a grid of molecular orbitals' values and save it in .cube format. All parameters for the grid domain are the same as MakeDensityCube. Additional input lines are required for the computation of molecular orbitals.\nThe user must specify which orbitals took take into account:\n All : **All**\n Occupied : **Occ**\n Virtual : **Virtual**\n Homo : **Homo**\n Lumo : **Lumo**\n Homo and lumo : **Homo, Lumo**\n Custom : **OrbitalsList=Orbital number specified by user**\nBy default the program will run with all MOs.\n\nThe choice of spin is also given:\n **SpinType=Alpha**\n **SpinType=Beta**\n **SpinType=Alpha, Beta**\n\nIf the user provides a custom list of orbitals the user can provide a list of spins corresponding to each orbital. This is done in **SpinList=alpha, beta, ...**.\nIf SpinList is shorter n length than OrbitalsList the program will fill the rest of the list with the last value read in the list" };
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

void Job::computeChargeNucleiContributions(const std::vector<Atom>& atoms, std::vector<std::vector<double>>& chargeNucleiContributions, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms, double nuclearCutoff)
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

void Job::computeHamiltonianMatrixWithPointCharges(const std::vector<ExcitedState>& states, const std::vector<double>& chargesNucleiContributions, const std::vector<std::vector<std::vector<std::vector<double>>>>& ionicMatrixes, std::vector<std::vector<double>>& psi_i_H_psi_j, std::vector<std::vector<double>>& psi_i_HminusH0_psi_j, std::ostream& outputStream, int verbose)
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
        std::string errorMessage = "Error in Job::computeHamiltonianMatrixWithPointCharges(): the first dimension of ionicMatrixes does not match the dimension of chargesNucleiContributions.";
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

void Job::computeIonicPotentialMatrixesFromBecke(Becke& becke, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& ionicPotentialMatrixes, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms, int kmax, int lebedev_order, int radial_grid_factor)
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

void Job::computeIonicPotentialMatrixesFromGrid(const Orbitals& orbitals, Grid& grid, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& ionicPotentialMatrixes, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms)
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

void Job::computeIonicPotentialMatrixesFromOrbitals(Orbitals& orbitals, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& ionicPotentialMatrixes, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms)
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

void Job::computeIonicPotentialVectorsFromBecke(Becke& becke, std::vector<std::vector<std::vector<std::vector<double>>>>& ionicPotentialVectors, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms, int kmax, int lebedev_order, int radial_grid_factor)
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

void Job::computeIonicPotentialVectorsFromOrbitals(Orbitals& orbitals, std::vector<std::vector<std::vector<std::vector<double>>>>& ionicPotentialVectors, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms)
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

void Job::computeLinearResponseFunctionMatrix(const Orbitals& orbitals, const std::vector<std::vector<std::vector<std::vector<double>>>>& tripleOrbitalIntegralMatrix, std::vector<std::vector<std::vector<double>>>& lrfMatrix)
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

void Job::computeResultsEnergyWithPointCharges(const std::vector<ExcitedState>& states, const std::vector<std::vector<double>>& psi_i_H_psi_j, const std::vector<std::vector<double>>& psi_i_HminusH0_psi_j, const std::string& outputFilePrefix, std::ostream& outputStream, int verbose)
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

void Job::computeResultsLinearResponseWithPointCharges(const std::vector<std::vector<double>>& eigenvalues, const std::vector<std::vector<std::vector<double>>>& ionicPotentialVectors, std::ostream& outputStream, int verbose)
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

void Job::printResultsEnergyWithPointCharges(const std::vector<ExcitedState>& states, const std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& ionicPotentialMatrixes, const std::vector<std::vector<double>>& chargeNucleiContributions, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms, const std::vector<Atom>& atoms, const std::string& outputPrefix, std::ostream& outputStream, int verbose)
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

                computeHamiltonianMatrixWithPointCharges(states, currentChargeNucleiContribution, currentIonicMatrix, psi_i_H_psi_j, psi_i_HminusH0_psi_j, outputStream, verbose);
                computeResultsEnergyWithPointCharges(states, psi_i_H_psi_j, psi_i_HminusH0_psi_j, outputPrefix + "_run" + runNumberStr, outputStream, verbose);

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

        computeHamiltonianMatrixWithPointCharges(states, currentChargeNucleiContribution, currentIonicMatrixes, psi_i_H_psi_j, psi_i_HminusH0_psi_j, outputStream, verbose);
        computeResultsEnergyWithPointCharges(states, psi_i_H_psi_j, psi_i_HminusH0_psi_j, outputPrefix, outputStream, verbose);
    }
}

void Job::printResultsLinearResponseWithPointCharges(const std::vector<std::vector<double>>& eigenvalues, const std::vector<std::vector<std::vector<std::vector<double>>>>& ionicPotentialVectors, const std::vector<double>& charges, const std::vector<std::array<double, 3>>& chargesPositions, bool loopOnAtoms, const std::vector<Atom>& atoms, std::ostream& outputStream, int verbose)
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

                computeResultsLinearResponseWithPointCharges(eigenvalues, currentIonicPotentialVectors, outputStream, verbose);

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
        
        computeResultsLinearResponseWithPointCharges(eigenvalues, currentIonicPotentialVectors, outputStream, verbose);

        logStream << std::defaultfloat << std::endl;
        log(logStream, outputStream);
    }
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
        case RunType::COMPUTE_ENERGY_WITH_POINT_CHARGES:
        {
            run_computeEnergyWithPointCharges();
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
        case RunType::COMPUTE_LINEAR_RESPONSE_WITH_POINT_CHARGES:
        {
            run_computeLinearResponseWithPointCharges();
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
        case RunType::HELP:
        {
            run_help();
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


