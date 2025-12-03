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
#include <Utils/Utils.h>


//----------------------------------------------------------------------------------------------------//
// SPECIFIC PARAMETERS READING FROM INPUT FILE
//----------------------------------------------------------------------------------------------------//

void Job::readAnalyticFilesNames(std::vector<std::string>& analyticFilesNames)
{
    if (!readListType<std::string>(_inputFile, "AnalyticFiles", analyticFilesNames))
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
}

void Job::readCutoff(double& cutoff)
{
    if (!readOneType<double>(_inputFile, "Cutoff", cutoff))
    {
        std::cout << "Warning: the \"Cutoff\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (Cutoff=0)." << std::endl;

        cutoff = 0.0;
    }
}

void Job::readELFMethod(ELFMethod& elfMethod)
{
    std::string strELFMethod;
    if (!readOneString(_inputFile, "ELFMethod", strELFMethod))
    {
        std::cout << "Warning: the \"ELFMethod\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (ELFMethod=Savin)." << std::endl;
        elfMethod = ELFMethod::SAVIN;
    }
    else
    {
        elfMethod = elfMethod_from_string(strELFMethod);
    }

    // Handle unknown ELF method: exit program with error message.
    if (elfMethod == ELFMethod::UNKNOWN)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: ELF method \"" << strELFMethod << "\" unknown." << std::endl;
        errorMessage << "Please check the documentation and the \"ELFMethod\" parameter value in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }
}

void Job::readEnergies(std::vector<double>& energies)
{
    if (!readListType<double>(_inputFile, "Energies", energies))
    {
        std::stringstream errorMessage;
        errorMessage << "Error: could not find energies." << std::endl;
        errorMessage << "Please check if the \"Energies\" parameter is defined and set in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }
}

void Job::readGridFilesNames(std::vector<std::string>& gridFilesNames)
{
    if (!readListType<std::string>(_inputFile, "Grids", gridFilesNames))
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
}

void Job::readOrbitalsNumbers(std::vector<int>& orbitalsNumbers)
{
    if (!readListType<int>(_inputFile, "OrbitalsNumbers", orbitalsNumbers))
    {
        std::stringstream errorMessage;
        errorMessage << "Error: could not find orbitals numbers." << std::endl;
        errorMessage << "Please check if the \"OrbitalsNumbers\" parameter is defined and set in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }
}

void Job::readOrbitalsSpins(std::vector<SpinType>& orbitalsSpins)
{
    std::vector<std::string> strOrbitalsSpins;
    if (!readListType<std::string>(_inputFile, "OrbitalsSpins", strOrbitalsSpins))
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

        orbitalsSpins.push_back(spin);
    }
}

void Job::readOrbitalType(OrbitalType& orbitalType)
{
    std::string strOrbitalType;
    if (!readOneString(_inputFile, "OrbitalType", strOrbitalType))
    {
        std::cout << "Warning: the \"OrbitalType\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (OrbitalType=All)." << std::endl;
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
}

void Job::readPartitionMethod(PartitionMethod& partitionMethod)
{
    std::string strMethod;
    if (!readOneString(_inputFile, "PartitionMethod", strMethod))
    {
        std::cout << "Warning: the \"PartitionMethod\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (PartitionMethod=On-Grid)." << std::endl;
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
}

void Job::readRunType(RunType& runType)
{
    std::string strRunType;
    if (!readOneString(_inputFile, "RunType", strRunType))
    {
        std::cout << "Warning: the \"RunType\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (RunType=Help)." << std::endl;
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
}

void Job::readSize(GridSize& gridSize, CustomSizeData& customSizeData)
{
    std::string strGridSize;
    if (!readOneString(_inputFile, "Size", strGridSize))
    {
        std::cout << "Warning: the \"Size\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (Size=Medium)." << std::endl;

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
}

void Job::readSpinList(std::vector<SpinType>& spinList)
{
    std::vector<std::string> strSpinList;
    if (!readListType<std::string>(_inputFile, "SpinList", strSpinList))
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
}

void Job::readSpinType(SpinType& spinType)
{
    std::string strSpinType;
    if (!readOneString(_inputFile, "SpinType", strSpinType))
    {
        std::cout << "Warning: the \"SpinType\" parameter is not specified in the provided input file (" << _inputFileName << ")." << std::endl;
        std::cout << "The program will use the default value (SpinType=Alpha-Beta)." << std::endl;

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
}

void Job::readTransitionsFileName(std::string& transitionsFileName)
{
    if (!readOneString(_inputFile, "TransitionsFile", transitionsFileName))
    {
        std::stringstream errorMessage;
        errorMessage << "Error: could not find transition file name." << std::endl;
        errorMessage << "Please check if the \"TransitionsFile\" parameter is defined and set in the provided input file (" << _inputFileName << ").";

        print_error(errorMessage.str());

        std::exit(1);
    }
}


//----------------------------------------------------------------------------------------------------//
// SECTION TO BE NAMED
//----------------------------------------------------------------------------------------------------//

void Job::setJobList()
{
    _jobsList = { "Help",
                  "computePartialCharges",
                  "computeDescriptors",
                  "computeIntegrals",
                  "computeGridDifference",
                  "MakeDensityCube",
                  "MakeOrbitalsCube",
                  "MakeELFCube",
                  "ConvertOrbitals" };
    
    _jobDescription = { "Details are given for the available jobs run by this program.\nExample input files for each job are also given. In this format, comment lines are specified by # at the start of the line",
                        "Grid-based computations of partial charges of the molecule. We provide 5 ways of computing atomic volumes, the first 3 of which are based on Bader's Atoms in molecule.\n\n **on-grid** : follows Tang's algorithm to find Bader volumes.\n **near-grid** : more precise version of on-grid.\n **near-grid-refinement** : even more precise. Requires more time.\n **VDD** topological method : assigns points to volumes by distance to closest atom.\n **Becke** : uses a regular density grid to interpolate Becke's atomic variable grids.\n\n Example format for input file :\n\n#RunType\n#RunType=Help\nRunType=ComputePartialCharges\n#GridFileName\nGrids=h2o_80_0.gcube \nPartitionMethod=on-grid\n\nW. Tang, E. Sanville, G. Henkelman, A grid-based bader analysis algorithm without lattice bias, Journal of Physics: Condensed Matter 21 (8) (2009) 084204.",
                        "Computation of chemical descriptors from analytic or cube files using on-grid, near-grid, near-grid-refinement and Becke. Frontier Molecular Orbitals(FMO) and finite difference(FD) are methods also provided for the computation. FMO requires 1 analytic file (.log, .wfx, .molden,...). FD requires 3 analytic files. The other methods require cube files of nucleophilic, electrophilic and radical attacks for the molecule. Energies must also be given by the user:\nif two are given, they are assumed to be the ionisation potential and the electronic affinity. If 3 are given they are assumed to be the total energies of each file. \n\n Example format for input file :\n\n#RunType=Help\n#RunType=ComputeDescriptorsFromCubes\n#GridFileName\nGrids=grid1.cube, grid2.cube, grid3.cube\nPartitionMethod=on-grid\nEnergies=I, A or E1,E2,E3",
                        "Compute local integrals of grids on volumes defined by method of choice. A grid is required to define the volumes.\nThe additional grids provided by the user should contain the quantities to be integrated.\n\n **on-grid** : to define volumes using on-grid AIM. Requires electronic density grid.\n **near-grid** : to define volumes using near-grid AIM. Requires electronic density grid.\n **near-grid-refinement : to define volumes using near-grid-refinement AIM. Requires electronic density grid.\n **VDD** : to define volumes by distance to atoms. Can use any type of density.\n **BBS** : Build Basins By SIGN. Requires a grid of density difference. A job is provided in the program to obtain such a grid. An additional input *Cutoff=* is required for BBS that sets a threshold for insignificant values.\n **B2S** : Build 2 basins by SIGN. Same as BBS but only constructs two volumes.\n\n Example format for input file :\n\n#RunType=Help\n#RunType=ComputeIntegrals\n#GridFileName\nGrids=gridDefiningVolumes.cube, grid1ToBeIntegrated.cube, grid2ToBeIntegrated.cube\nPartitionMethod=BBS\nCutoff=1e-10",
                        "Computes the differences of values of the first two grids provides and assigns them to the third.\n\n Example format for input file : \n\n#Runtype=Help\nRunType=ComputeDifference\n#GridFileName\nGrids=in1.cube, in2.cube, out.cube ",
                        "Create a density grid and save it in .cube format. .wfx , .fchk , .molden , .gab and .log are supported as input files.\nthe user can choose from 3 standard grid sizes:\ncoarse ( 3 pts / Bohr)\nMedium (6 pts / Bohr)\nFine (12 pts / Bohr)\n\nA custom size is also provided in which the user enters the domain data as follows:\nNx, Ny, Nz, Ox, Oy, Oz, T11, T12, T13, T21, T22, T23, T31, T32, T33\nWhere N is the number of points in the ith direction, Oi are the coordinates of the bottom left corner of the cube and Tij are the coeficients of the translation vector.\n\n Example format for input file : \n\n#RunType=Help\nRunType=MakeDensityCube\n#GridFileName\nAnalyticFile=filename.wfx\nSize=Custom\nCustomSizeData=80,80,80,5,5,5,0.15,0,0,0,0.15,0,0,0,0.15\nGrid=save.cube ",
                        "Compute a grid of molecular orbitals' values and save it in .cube format. All parameters for the grid domain are the same as MakeDensityCube. Additional input lines are required for the computation of molecular orbitals.\nThe user must specify which orbitals took take into account:\n All : **All**\n Occupied : **Occ**\n Virtual : **Virtual**\n Homo : **Homo**\n Lumo : **Lumo**\n Homo and lumo : **Homo, Lumo**\n Custom : **OrbitalsList=Orbital number specified by user**\nBy default the program will run with all MOs.\n\nThe choice of spin is also given:\n **SpinType=Alpha**\n **SpinType=Beta**\n **SpinType=Alpha, Beta**\n\nIf the user provides a custom list of orbitals the user can provide a list of spins corresponding to each orbital. This is done in **SpinList=alpha, beta, ...**.\nIf SpinList is shorter n length than OrbitalsList the program will fill the rest of the list with the last value read in the list",
                        "Create a grid and compute the Electron Localisation Function (ELF) using either Savin or Becke method. Grid domain is defined the same as the MakeDensityCube.\nBy default the program will run Savin ELF.\n\n Example format for input file : \n\n#RunType=Help\nRunType=MakeELFCube\n#GridFileName\nAnalyticFile=filename.wfx\nSize=Medium\nELFmethod=Becke\nGrid=save.cube",
                        "Convert Analytical file.\nSupported file formats are : wfx, fchk, log, molden, gab.\nOutput supported : wfx, molden, gab\n\n Example format for input file : \n\n#RunType=Help\nRunType=ConvertOrbitals\nAnalyticFiles=input.wfx, output.molden" };
}

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
/******************************************************************************************/
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


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------//

Job::Job():_inputFileName("input.txt")
{
    setJobList();
    openInputFile();
}
/******************************************************************************************/
Job::Job(std::string inputFileName):_inputFileName(inputFileName)
{
    setJobList();
    openInputFile();
}
/******************************************************************************************/
Job::~Job()
{
    _inputFile.close();
}
/******************************************************************************************/
void Job::buildBasins(GridCP &gcp, const std::string &GridFileName, PartitionMethod partitionMethod)
{
    std::ifstream gridf(GridFileName);
    Grid g(gridf, _table);    
    gcp.buildBasins(g, partitionMethod);
}
void Job::buildBasinsBySign(GridCP& gcp, const std::string& GridFileName,  double cutoff, bool two) 
{
    std::ifstream gridf(GridFileName);
    Grid g(gridf, _table);
    gridf.close();
    if (two)
    {
            gcp.build2BasinSign(g);
    }    
    else
    {
            gcp.buildBasinsBySign(g, cutoff);
    }    
}
void Job::computeLocalIntegrals(GridCP& gcp, const std::vector<std::string>& GridFileNames)
{
    for(size_t i=0;i<GridFileNames.size();i++)
    {
        std::ifstream f(GridFileNames[i]);
        Grid g(f,_table);
        f.close();
        gcp.computeIntegrals(g);
        gcp.printCriticalPoints();
    }
}
void Job::ComputeGridDifference(const std::string& GFN1, const std::string& GFN2 , const std::string& NameNew)
{
    std::ifstream in1(GFN1);
    Grid g(in1, _table);

    std::ifstream in2(GFN2);
    Grid h(in2, _table);

    Grid diff = g-h;
    
    in1.close();
    in2.close();

    std::ofstream out(NameNew, std::ios::out);

    if (out.fail())
    {
        std::cout << "Failed to write to file " << NameNew << "..." << std::endl;
    }

    diff.save(out);
    std::cout <<"Grid has been saved to : " << NameNew << std::endl;
}

void Job::printCriticalPoints()
{
    std::cerr << "Function Job::printCriticalPoints() not implemented yet." << std::endl;
}

std::vector<double> Job::computePartialCharges(const std::string &gridfname, PartitionMethod partitionMethod)
{
    std::vector<double> charges(3);

    if (partitionMethod == PartitionMethod::BECKE)
    {
        Factorial fact(100);
        Binomial bino (100, fact);
        std::ifstream gridf(gridfname);
        Grid grid(gridf, _table);

        Becke B(grid);
        B.partial_charge(grid);
        charges = B.get_Partial_Charge();
        B.printCharges();
    }
    else
    {
        std::ifstream gridf(gridfname);
        Grid grid(gridf, _table);
        GridCP gridcp;

        gridcp.buildBasins(grid, partitionMethod);
        gridcp.computeIntegrals(grid);
        gridcp.printCriticalPoints();

        charges = gridcp.computeAIMCharges(grid);
    }

    return charges;
}
std::vector<double> Job::computePartialChargesAndEnergy(std::vector<double>& E, const std::string& ANAFileName)
{
    Becke B;
    readFileFormat<Becke>(B, ANAFileName);
    E.push_back(B.PartialChargesAndEnergy()[0][0]);
    return B.PartialChargesAndEnergy()[1];
}
Structure Job::returnStruct(const std::string& ANAFileName)
{
    Becke B;
    readFileFormat<Becke>(B, ANAFileName);
    return B.get_struct();
}
Descriptors Job::computeDescriptors(const std::string &GridFileName1, const std::string &GridFileName2, const std::string &GridFileName3, double I, double A, PartitionMethod partitionMethod)
{
    std::ifstream grid1(GridFileName1);
    std::ifstream grid2(GridFileName2);
    std::ifstream grid3(GridFileName3);
    Descriptors D(grid1, grid2, grid3, I, A, partitionMethod);
    return D;
}
Descriptors Job::computeDescriptors(const std::string &GridFileName1, const std::string &GridFileName2, const std::string &GridFileName3, std::vector<double> E, PartitionMethod partitionMethod)
{
    std::ifstream grid1(GridFileName1);
    std::ifstream grid2(GridFileName2);
    std::ifstream grid3(GridFileName3);
    Descriptors D(grid1, grid2, grid3, E, partitionMethod);
    return D;
}
void Job::computeDescriptorsFD(const std::string& ANAFileName1, const std::string& ANAFileName2, const std::string& ANAFileName3)
{
    std::vector<double> E(0);
    std::vector<double> Q1 = computePartialChargesAndEnergy(E, ANAFileName1);
    std::vector<double> Q2 = computePartialChargesAndEnergy(E, ANAFileName2);
    std::vector<double> Q3 = computePartialChargesAndEnergy(E, ANAFileName3);

    Structure s = returnStruct(ANAFileName1);

    Descriptors D(s, Q1, Q2, Q3, E);
    std::cout << D;
}
template<typename T> void Job::readFileFormat(T& AnaClass, const std::string& ANAFileName)
{
    if (ANAFileName.find(".wfx") != std::string::npos)
    {
        std::cout << "Reading data from " << ANAFileName << "... Please wait." << std::endl;
        AnaClass = computeOrbOrBecke<WFX,T>(ANAFileName);
    }
    else if (ANAFileName.find(".fchk") != std::string::npos)
    {
        std::cout << "Reading data from " << ANAFileName << "... Please wait." << std::endl;
        AnaClass = computeOrbOrBecke<FCHK,T>(ANAFileName);
    }
    else if (ANAFileName.find(".molden") != std::string::npos)
    {
        std::cout << "Reading data from " << ANAFileName << "... Please wait." << std::endl;
        AnaClass = computeOrbOrBecke<MOLDENGAB,T>(ANAFileName);
    }
    else if (ANAFileName.find(".gab") != std::string::npos)
    {
        std::cout << "Reading data from " << ANAFileName << "... Please wait." << std::endl;
        AnaClass = computeOrbOrBecke<MOLDENGAB,T>(ANAFileName);
    }
    else if (ANAFileName.find(".log") != std::string::npos)
    {
        std::cout << "Reading data from " << ANAFileName << "... Please wait." << std::endl;
        AnaClass = computeOrbOrBecke<LOG,T>(ANAFileName);
    }
    else
    {
        std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        std::cerr << "Sorry, unknown file format for analytic file. Please check input file. " << std::endl;
        std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        
        std::exit(1);
    }
}
template<typename T,typename U> U Job::computeOrbOrBecke(const std::string& analyticFileName)
{
    Factorial fact(100);
    Binomial bino(100, fact);
    std::ifstream analyticFile(analyticFileName);
    T anaClass(analyticFile);
    analyticFile.close();
    U OrbOrBecke(anaClass, bino, _table);
    return OrbOrBecke;
}
//TypeFlag specifies the type of grid you wnat to make. For now there are 3 types available. electronic density, ELF and orbitals. Others can be added in the else ifs. additional parameters shoud be added before the default values.
void Job::createCube(Orbitals& orb, const Domain& d, const std::string& cubeFileName, int TypeFlag, const ELFMethod elfMethod, std::vector<int> nums, std::vector<int> typesSpin)
{
    Grid g;
    if (TypeFlag == 0)
    {
        g=orb.makeGrid(d);

    }
    else if (TypeFlag == 1)
    {
        g = orb.makeOrbGrid(d, nums, typesSpin);
    }
    else
    {
        if (elfMethod == ELFMethod::BECKE)
        {
            g = orb.makeELFgrid(d, 0);
        }
        else // SAVIN
        {
            g=orb.makeELFgrid(d); 
        }

    }
    std::ofstream out(cubeFileName);
    g.save(out);
    out.close();
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
        std::vector<double> tx = {2 * xmax / (N1 - 1), 0, 0};
        std::vector<double> ty = {0, 2 * ymax / (N2 - 1), 0};
        std::vector<double> tz = {0, 0, 2 * zmax / (N3 - 1)};

        std::vector<std::vector<double>>  t = {tx, ty, tz};

        d.set_all(Nval, N1, N2, N3, xmax, ymax, zmax, t);
    }
    else
    {
        std::vector<std::vector<double>> t(0);

        std::vector<double> dx = {customSizeData[6], customSizeData[7], customSizeData[8]};
        std::vector<double> dy = {customSizeData[9], customSizeData[10], customSizeData[11]};
        std::vector<double> dz = {customSizeData[12], customSizeData[13], customSizeData[14]};

        t.push_back(dx);
        t.push_back(dy);
        t.push_back(dz);

        d.set_all(Nval, int(customSizeData[0]), int(customSizeData[1]), int(customSizeData[2]), customSizeData[3], customSizeData[4], customSizeData[5], t);
    }
    return d;
}    
/******************************************************************************************/

void Job::setOrbitals(Orbitals& o, const int numberOfOrbitals, vector<int>& orbitalsNumbers, vector<int>& orbitalsSpins, const OrbitalType orbitalType, SpinType spinType, const std::vector<SpinType>& spinList)
{
    switch (orbitalType)
    {
        case OrbitalType::ALL:
        {
            setAllOrb(orbitalsNumbers, orbitalsSpins, o, spinType, numberOfOrbitals);
            break;
        }
        case OrbitalType::CUSTOM:
        {
            setCustom(orbitalsNumbers, orbitalsSpins, spinList);
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
            setOccOrb(orbitalsNumbers, orbitalsSpins, o, spinType, numberOfOrbitals);
            break;
        }
        case OrbitalType::VIRTUAL:
        {
            setVirtOrb(orbitalsNumbers, orbitalsSpins, o, spinType, numberOfOrbitals);
            break;
        }
        default:
        {
            // should not happen: verified before function called
            std::cerr << "Error: Unknown orbital type encountered in setOrbitals(). This should not happen as the orbital type is verified before this function is called." << std::endl;
        }
    }
}

void Job::setAllOrb(std::vector<int>& orbnums, std::vector<int>& orbspin, Orbitals& o, SpinType spinType, const int& N)
{
    orbnums.resize(N);
    for(int i = 0; i < N; ++i)
    {
        orbnums[i] = i;
    }

    if (!o.AlphaAndBeta() and spinType == SpinType::ALPHA)
    {
        orbspin.resize(N, 0);
    }
    else if (!o.AlphaAndBeta() and spinType == SpinType::ALPHA_BETA)
    {
        orbspin.resize(N, 0);
        orbnums.resize(2 * N, 1);
        for(int i = N; i < 2 * N; ++i)
        {
            orbnums[i] = i - N;
        }
    }
    else
    {
        orbspin.resize(N, 1);
    }
}
void Job::setOccOrb(std::vector<int>& orbnums, std::vector<int>& orbspin, Orbitals& o, SpinType spinType, const int& N)
{
    std::vector<std::vector<double>> occ=o.OccupationNumber();
    if (!o.AlphaAndBeta() and spinType == SpinType::ALPHA) 
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
        orbspin.resize(k, 0);
    }
    else if (!o.AlphaAndBeta() and spinType == SpinType::ALPHA_BETA) 
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
        orbspin.resize(k, 0);
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
        orbspin.resize(j, 1);    
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
        orbspin.resize(j, 1);
    }

}
void Job::setVirtOrb(std::vector<int>& orbnums, std::vector<int>& orbspin, Orbitals& o, SpinType spinType, const int& N)
{
    std::vector<std::vector<double>> occ=o.OccupationNumber();
    if (!o.AlphaAndBeta() and spinType == SpinType::ALPHA)
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
        orbspin.resize(k, 0);
    }
    else if (!o.AlphaAndBeta() and spinType == SpinType::ALPHA_BETA) 
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
        orbspin.resize(k, 0);
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
        orbspin.resize(j, 1);    
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
        orbspin.resize(j, 1);
    }
}
void Job::setHomo(std::vector<int>& orbnums, std::vector<int>& orbspin, Orbitals& o, SpinType spinType)
{
    int i = 0;
    std::vector<std::vector<double>> occ = o.OccupationNumber();
    if (!o.AlphaAndBeta() and spinType == SpinType::ALPHA)
    {
        while(occ[0][i] > 1e-10)
        {
            i++;
        }
        orbnums = {i - 1};
        orbspin = {0};
    }
    else if (!o.AlphaAndBeta() and spinType == SpinType::ALPHA_BETA)
    {
        int j = 0;
        while(occ[0][i] > 1e-10)
        {
            i++;
        }
        while(occ[1][j] > 1e-10)
        {
            j++;
        }
        orbnums = {i - 1, j - 1};
        orbspin = {0, 1};
    }
    else
    {
        int j = 0;
        while(occ[1][j] > 1e-10)
        {
            j++;
        }
        orbnums = {j - 1};
        orbspin = {1};
    }
}
void Job::setLumo(std::vector<int>& orbnums, std::vector<int>& orbspin, Orbitals& o, SpinType spinType)
{
    int i = 0;
    std::vector<std::vector<double>> occ = o.OccupationNumber();
    if (!o.AlphaAndBeta() and spinType == SpinType::ALPHA)
    {
        while(occ[0][i] > 1e-10)
        {
            i++;
        }
        orbnums = {i};
        orbspin = {0};
    }
    else if (!o.AlphaAndBeta() and spinType == SpinType::ALPHA_BETA)
    {
        int j = 0;
        while(occ[0][i] > 1e-10)
        {
            i++;
        }
        while(occ[1][j] > 1e-10)
        {
            j++;
        }
        orbnums = {i, j};
        orbspin = {0, 1};
    }
    else
    {
        int j = 0;
        while(occ[1][j] > 1e-10)
        {
            j++;
        }
        orbnums = {j};
        orbspin = {1};
    }
}
void Job::setHomoLumo(std::vector<int>& orbnums, std::vector<int>& orbspin, Orbitals& o, SpinType spinType) 
{
    int i = 0;
    std::vector<std::vector<double>> occ = o.OccupationNumber();
    if (!o.AlphaAndBeta() and spinType == SpinType::BETA)
    {
        while(occ[0][i] > 1e-10)
        {
            i++;
        }        
        orbnums = {i - 1, i};
        orbspin = {0, 0};
    }
    else if (!o.AlphaAndBeta() and spinType == SpinType::ALPHA_BETA)
    {
        int j = 0;
        while(occ[0][i] > 1e-10)
        {
            i++;
        }
        while(occ[1][j] > 1e-10)
        {
            j++;
        }
        orbnums = {i - 1, j - 1, i, j};
        orbspin = {0, 1, 0, 1};
    }
    else
    {
        int j = 0;
        while(occ[1][j] > 1e-10)
        {
            j++;
        }
        orbnums = {j - 1, j};
        orbspin = {1, 1};
    }
}
void Job::setCustom(std::vector<int>& orbnums, std::vector<int>& orbspin, const std::vector<SpinType>& spinList) 
{
    for(size_t i=0; i < spinList.size();i++)
    {
        if (spinList[i] == SpinType::ALPHA)
        {
            orbspin.push_back(0);
        }
        else if (spinList[i] == SpinType::BETA)
        {
            orbspin.push_back(1);
        }
    }
    for(size_t i=0; i<orbnums.size(); i++)
    {
        orbnums[i] -= 1;
    }
    if (orbspin.size()<orbnums.size())
    {    
        int last=orbspin.back();
        for(size_t i=orbspin.size(); i<orbnums.size(); i++)
        {
            orbspin.push_back(last);
        }
    }
}


//----------------------------------------------------------------------------------------------------//
// RUN METHOD
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
        case RunType::CONVERT_ORBITALS:
        {
            run_convertOrbitals();
            break;
        }
        case RunType::COMPUTE_PARTIAL_CHARGES:
        {
            run_computePartialCharges();
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


//----------------------------------------------------------------------------------------------------//
// JOBS
//----------------------------------------------------------------------------------------------------//

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
            computeDescriptorsFD(analyticFilesNames[0], analyticFilesNames[1], analyticFilesNames[2]);
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
            readFileFormat<Orbitals>(o, analyticFilesNames[0]);
            o.PrintDescriptors();
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
            std::cout << "Reading Ionisation potential I = " << energies[0] << " and Electronic affinity A = " << energies[1] << std::endl;
            
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
    
    ComputeGridDifference(gridFilesNames[0], gridFilesNames[1], gridFilesNames[2]);
    
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
        readFileFormat<Orbitals>(o, analyticFilesNames[0]);
        o.Save(analyticFilesNames[1]);
    }
    else
    {
        std::cout << "Warning: input and output files have the same format (" << inputFileExtension << "). Nothing to be done." << std::endl;
    }
}

void Job::run_help()
{
    printListOfRunTypes();
}

void Job::run_lambdaDiagnostic()
{
    //Read analytic file
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
    readFileFormat<Orbitals>(orbitals, analyticFilesNames[0]);


    // Setting orbitals
    int numberOfOrbitals = orbitals.NumberOfMo();
    std::vector<int> orbitalsSpins;
    std::vector<int> orbitalsNumbers;

    setOrbitals(orbitals, numberOfOrbitals, orbitalsNumbers, orbitalsSpins, OrbitalType::ALL, SpinType::ALPHA_BETA);


    // Building domain and grid
    std::cout << "Building domain and grid, please wait..." << std::endl;

    Domain domain = buildDomainForCube(orbitals, gridSize, customSizeData, orbitalsNumbers.size());
    Grid orbitalsGrid = orbitals.makeOrbGrid(domain, orbitalsNumbers, orbitalsSpins);


    // Lecture du fichier des transitions
    std::vector<ExcitedState> excitedStates;
    ExcitedState::readTransitionsFile(transitionsFileName, excitedStates);
    std::cout << "Number of excited states read: " << excitedStates.size() << std::endl << std::endl;

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
    readFileFormat<Orbitals>(o, analyticFilesNames[0]);
    

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
    readFileFormat<Orbitals>(o, analyticFilesNames[0]);


    // Setting orbitals
    int numberOfOrbitals = o.NumberOfMo();
    std::vector<int> orbitalsSpins;
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

    createCube(o, d, gridFileName[0], 1, ELFMethod::UNKNOWN, orbitalsNumbers, orbitalsSpins);

    std::cout << "Data saved to file: " << gridFileName[0] << std::endl;
}

void Job::run_makeELFCube()
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
    readFileFormat<Orbitals>(o, analyticFilesNames[0]);


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

    createCube(o, d, gridFilesName[0], 2, elfMethod);

    std::cout << "ELF cube saved to file " << gridFilesName[0] << '.' << std::endl;
}


