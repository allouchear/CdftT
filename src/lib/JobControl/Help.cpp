#include <string>
#include <unordered_map>

#include <JobControl/Help.hpp>


//----------------------------------------------------------------------------------------------------//
// STATIC FIELDS
//----------------------------------------------------------------------------------------------------//

std::unordered_map<RunType, std::string> Help::_s_availableJobs = 
{
    {
        RunType::COMPUTE_CONDENSED_LINEAR_RESPONSE,
        "ComputeCondensedLinearResponse"
    },
    {
        RunType::COMPUTE_DESCRIPTORS,
        "Computation of chemical descriptors from analytic or cube files using on-grid, near-grid, near-grid-refinement and Becke."
        "Frontier Molecular Orbitals (FMO) and finite difference (FD) are methods also provided for the computation."
        "FMO requires 1 analytic file (.log, .wfx, .molden,...). FD requires 3 analytic files."
        "The other methods require cube files of nucleophilic, electrophilic and radical attacks for the molecule.\n\n"

        "Energies must also be given by the user. If two are given, they are assumed to be the ionization energy and the electron affinity."
        "If 3 are given they are assumed to be the total energies of each file.\n\n"
        
        "Example format for input file :\n\n"

        "#RunType=Help\n"
        "#RunType=ComputeDescriptorsFromCubes\n"
        "#GridFileName\n"
        "Grids=grid1.cube, grid2.cube, grid3.cube\n"
        "PartitionMethod=on-grid\n"
        "Energies=I, A or E1,E2,E3"
    },
    {
        RunType::COMPUTE_ENERGY_WITH_POINT_CHARGES,
        "Computes the new energy levels of a system when one or many point charges are added."
        "This computation can be performed analytically, on a cube grid or on a Becke grid."
        "The number of energy levels computed is equal to the number of states given in the input file that describes the electronic transitions in the unperturbed system.\n\n"

        "If the positions of the charges are not given in the input file, they will be sequentially placed on the system's atom nuclei."
    },
    {
        RunType::COMPUTE_GRID_DIFFERENCE,
        "Computes the differences of values of the first two grids provides and assigns them to the third.\n\n"
        
        "Example format for input file :\n\n"
        
        "#Runtype=Help\n"
        "RunType=ComputeDifference\n"
        "#GridFileName\n"
        "Grids=in1.cube, in2.cube, out.cube"
    },
    {
        RunType::COMPUTE_INTEGRALS,
        "Computes local integrals of grids on volumes defined by method of choice. A grid is required to define the volumes.\n"
        "The additional grids provided by the user should contain the quantities to be integrated.\n\n"
        
        "- on-grid: to define volumes using on-grid AIM. Requires electronic density grid.\n"
        "- near-grid: to define volumes using near-grid AIM. Requires electronic density grid.\n"
        "- near-grid-refinement: to define volumes using near-grid-refinement AIM. Requires electronic density grid.\n"
        "- VDD: to define volumes by distance to atoms. Can use any type of density.\n"
        "- BBS: Build Basins By SIGN. Requires a grid of density difference."
        "Note that there is a job in CdftT to obtain such a grid. An additional input *Cutoff=* is required for BBS that sets a threshold for insignificant values.\n"
        "- B2S: Build 2 basins by SIGN. Same as BBS but only constructs two volumes.\n\n"
        
        "Example format for input file :\n\n"

        "#RunType=Help\n"
        "#RunType=ComputeIntegrals\n"
        "#GridFileName\n"
        "Grids=gridDefiningVolumes.cube, grid1ToBeIntegrated.cube, grid2ToBeIntegrated.cube\n"
        "PartitionMethod=BBS\n"
        "Cutoff=1e-10"
    },
    {
        RunType::COMPUTE_PARTIAL_CHARGES,
        "Grid-based computations of partial charges of the molecule."
        "We provide 5 ways of computing atomic volumes. The first 3 methods are based on Bader's Atoms in molecule.\n\n"
        
        "- on-grid: follows Tang's algorithm to find Bader volumes.\n"
        "- near-grid: more precise version of on-grid.\n"
        "- near-grid-refinement: even more precise. Requires more time.\n"
        "- VDD: topological method : assigns points to volumes by distance to closest atom.\n"
        "- Becke: uses a regular density grid to interpolate Becke's atomic variable grids.\n\n"
        
        "Example format for input file :\n\n"
        
        "#RunType=Help\n"
        "RunType=ComputePartialCharges\n"
        "#GridFileName\n"
        "Grids=h2o_80_0.gcube \n"
        "PartitionMethod=on-grid\n\n"

        "Reference:\n"
        "W. Tang, E. Sanville, G. Henkelman, A grid-based bader analysis algorithm without lattice bias, Journal of Physics: Condensed Matter 21 (8) (2009) 084204."
    },
    {
        RunType::CONVERT_ORBITALS,
        "Convert Analytical file. Supported file formats are: wfx, fchk, log, molden, gab.\n"
        "Supported output formats are: wfx, molden, gab.\n\n"
        
        "Example format for input file : \n\n"
        
        "#RunType=Help\n"
        "RunType=ConvertOrbitals\n"
        "AnalyticFiles=input.wfx, output.molden"
    },
    {
        RunType::HELP,
        "Details are given for the available jobs run by this program.\n"
        "Example input files for each job are also given. In this format, comment lines are specified by a hash character (#) at the start of the line"
    },
    {
        RunType::LAMBDA_DIAGNOSTIC,
        "Prints the result of the Lambda diagnostic test, as described by Peach et al., that judges the reliability of TDDFT excited states calculations."
        "It also allows to validate the grid size configuration by computing overlap integrals between the orbitals involved in the excited states.\n\n"
        
        "Example format for input file : \n\n"

        "#RunType=Help\n"
        "RunType=LambdaDiagnostic\n"
        "AnalyticFile=filename.wfx\n"
        "Size=Custom\n"
        "CustomSizeData=128,128,128,5,5,5,7.87401574803148e-02,0,0,0,7.87401574803148e-02,0,0,0,7.87401574803148e-02"
        "TransitionsFile=transitions.txt"
    },
    {
        RunType::MAKE_DENSITY_CUBE,
        "Creates a density grid and saves it in .cube format. Supported input formats are: .wfx , .fchk , .molden , .gab and .log.\n"
        "Three standard grid sizes are available:\n"
        "- coarse ( 3 pts / Bohr)\n"
        "- medium (6 pts / Bohr)\n"
        "- fine (12 pts / Bohr)\n\n"

        "A custom size is also provided in which the user enters the domain data as follows:\n"
        "Nx, Ny, Nz, Ox, Oy, Oz, T11, T12, T13, T21, T22, T23, T31, T32, T33\n"
        "Where N is the number of points in the ith direction, Oi are the coordinates of the bottom left corner of the cube and Tij are the coeficients of the translation vector.\n\n"

        "Example format for input file : \n\n"

        "#RunType=Help\n"
        "RunType=MakeDensityCube\n"
        "#GridFileName\n"
        "AnalyticFile=filename.wfx\n"
        "Size=Custom\n"
        "CustomSizeData=128,128,128,5,5,5,7.87401574803148e-02,0,0,0,7.87401574803148e-02,0,0,0,7.87401574803148e-02"
        "Grid=save.cube"
    },
    {
        RunType::MAKE_ELF_CUBE,
        "Creates a grid and compute the Electron Localisation Function (ELF) using either Savin or Becke method."
        "The grid's domain is defined the same as the MakeDensityCube.\n"
        "By default the program will run Savin ELF.\n\n"
        
        "Example format for input file : \n\n"
        
        "#RunType=Help\n"
        "RunType=MakeELFCube\n"
        "#GridFileName\n"
        "AnalyticFile=filename.wfx\n"
        "Size=Medium\n"
        "ELFmethod=Becke\n"
        "Grid=save.cube"
    },
    {
        RunType::MAKE_ORBITALS_CUBE,
        "Computes a grid of molecular orbitals values and saves it in a .cube format."
        "All parameters for the grid's domain are the same as MakeDensityCube."
        "Additional input lines are required for the computation of molecular orbitals.\n"
        "The user must specify which orbitals took take into account:\n"
        "- all: all MOs\n"
        "- occ: occupied MOs"
        "- virtual: Virtual MOs\n"
        "- homo: Highest Occupied MO\n"
        "- lumo: Lowest Unoccupied MO\n"
        "- homo-lumo: HOMO and LUMO\n"
        "- custom: user-selected MOs. In this case, two parameters must be included:\n"
        "    - OrbitalsList: comma-separated list of MO numbers (starting from 1)\n"
        "    - OrbitalsSpins: comma-separated list of spins types (Alpha or Beta).\n"
        "By default the program will run with all MOs.\n\n"
        
        "Note that if the spin list is shorter in length than the orbital numbers list, the program will fill the rest of the spin list with the last value read from that list."
    }
};


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//----------------------------------------------------------------------------------------------------//

Help::Help():
    Job()
{ }

//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

void Help::run()
{
    std::cout << "Available jobs (runType=) :" << std::endl << std::endl;

    for(const auto& pair : _s_availableJobs)
    {
        std::cout << "--------------------------------------------------" << std::endl;
        std::cout << to_string(pair.first) << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;
        std::cout << pair.second << std::endl << std::endl;
    }
}