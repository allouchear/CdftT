#ifndef CDFTT_JOB_H_INCLUDED
#define CDFTT_JOB_H_INCLUDED

#include <string>
#include <fstream>
#include <set>
#include <vector>

#include <Becke/Becke.h>
#include <Common/Descriptors.h>
#include <Common/PeriodicTable.h>
#include <Cube/Grid.h>
#include <Cube/GridCP.h>
#include <Orbitals/ExcitedState.hpp>
#include <Orbitals/Orbitals.h>
#include <Utils/Enums.hpp>
#include <Utils/Utils.h>


/**
 * @brief Job class.
 * 
 * Manages user-available jobs for the CdtfT program.
 */
class Job
{
    protected:
        /** @brief Name of the input file to parse. */
        std::string _inputFileName;

        /** @brief Input file stream opened on the input file. */
        std::ifstream _inputFile;


        //----------------------------------------------------------------------------------------------------//
        // READING PARAMETERS FROM INPUT FILE
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Reads the name(s) of the analytic file(s) from the "AnalyticFiles" parameter in the input file.
         *
         * @param[out] analyticFilesNames Reference to a vector where the read filename(s) will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readAnalyticFilesNames(std::vector<std::string>& analyticFilesNames);

        /**
         * @brief Reads the Becke parameters from the "Becke" parameter in the input file.
         *
         * @param[out] beckeParameters Reference to a vector where the read Becke parameters will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readBecke(std::vector<int>& beckeParameters);

        /** @brief Reads the charges of the point charges from the "Charges" parameter in the input file.
         * 
         * @param[out] charges Reference to a vector of doubles where the read charge values will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readCharges(std::vector<double>& charges);

        /** @brief Reads the numeric cutoff used by some partitioning methods from the "Cutoff" parameter in the input file.
         * 
         * @param[out] cutoff Reference to a double where the read cutoff value will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readCutoff(double& cutoff);

        /**
         * @brief Reads the method used to compute the Electron Localization Function (ELF) from the "ELFMethod" parameter in the input file.
         *
         * @param[out] elfMethod Reference to an ELFMethod variable where the parsed method will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readELFMethod(ELFMethod& elfMethod);

        /**
         * @brief Reads the energies from the "Energies" parameter in the input file.
         *
         * @param[out] energies Reference to a vector where the energy values will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readEnergies(std::vector<double>& energies);

        /**
         * @brief Reads the energy point charge methods from the "EnergyPointChargeMethods" parameter in the input file.
         * 
         * @param[out] energyPointChargeMethods Reference to a vector where the parsed methods will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readEnergyPointChargeMethods(std::vector<EnergyPointChargeMethod>& energyPointChargeMethods);

        /**
         * @brief Reads a list of state numbers from the "ExcitedStatesNumbers" parameter in the input file.
         *
         * @param[out] excitedStatesNumbers Reference to a vector where the excited state numbers will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readExcitedStatesNumbers(std::vector<int>& excitedStatesNumbers);

        /**
         * @brief Reads a list of orbital numbers to exclude from the density computation from the "ExcludedOrbitals" parameter in the input file.
         *
         * @param[out] excludedOrbitals Reference to a vector where the excluded orbital numbers will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readExcludedOrbitals(std::vector<int>& excludedOrbitals);

        /**
         * @brief Reads the name(s) of the grid file(s) (.cube files) from the "GridFiles" parameter in the input file.
         *
         * @param[out] gridFilesNames Reference to a vector where the read filename(s) will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readGridFilesNames(std::vector<std::string>& gridFilesNames);

        /**
         * @brief Reads the energy of the ground state from the "GroundStateEnergy" parameter in the input file.
         * 
         * @param[out] energy Reference to a double where the read energy value will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readGroundStateEnergy(double& energy);

        /**
         * @brief Reads the maximum number of excited states to consider from the "MaxNumberOfExcitedStates" parameter in the input file.
         *
         * @param[out] maxNumberOfExcitedStates Reference to an integer where the read value will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readMaxNumberOfExcitedStates(int& maxNumberOfExcitedStates);

        /**
         * @brief Reads the cutoff distance for nuclear contribution from the "NuclearCutoff" parameter in the input file.
         *
         * @param[out] nuclearCutoff Reference to a double where the read cutoff distance will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readNuclearCutoff(double& nuclearCutoff);

        /**
         * @brief Reads the option that indicates if each charge is put on a different nucleus from the "OneChargePerNucleus" parameter in the input file.
         * 
         * @param[out] oneChargePerNucleus Reference to a boolean where the read option value will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readOneChargePerNucleus(bool& oneChargePerNucleus);

        /**
         * @brief Reads the option that indicates if each charge has a single position from the "OnePositionPerCharge" parameter in the input file.
         * 
         * @param[out] onePositionPerCharge Reference to a boolean where the read option value will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readOnePositionPerCharge(bool& onePositionPerCharge);

        /**
         * @brief Reads a list of orbital numbers from the "OrbitalsNumbers" parameter in the input file.
         *
         * @param[out] orbitalsNumbers Reference to a vector where the orbital numbers will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readOrbitalsNumbers(std::vector<int>& orbitalsNumbers);

        /**
         * @brief Reads a list of orbitals spins from the "OrbitalsSpins" parameter in the input file.
         *
         * @param[out] orbitalsSpins Reference to a vector where the parsed spin types will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readOrbitalsSpins(std::vector<SpinType>& orbitalsSpins);

        /**
         * @brief Reads the selected orbital type from the "OrbitalType" parameter in the input file.
         *
         * @param[out] orbitalType Reference to an OrbitalType variable where the parsed selected orbital type will be stored.
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readOrbitalType(OrbitalType& orbitalType);

        /**
         * @brief Reads the output filename prefix from the "OutputPrefix" parameter in the input file.
         *
         * @param[out] outputPrefix Reference to a string where the read output prefix will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readOutputPrefix(std::string& outputPrefix);

        /**
         * @brief Reads the partition method from the "PartitionMethod" parameter in the input file.
         *
         * @param[out] partitionMethod Reference to a PartitionMethod variable where the parsed method will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readPartitionMethod(PartitionMethod& partitionMethod);

        /**
         * @brief Reads 3-coordinate positions from the "Positions" parameter in the input file.
         *
         * @param[out] positions Reference to a vector of arrays of three doubles where the read positions will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readPositions(std::vector<std::array<double, 3>>& positions);

        /**
         * @brief Reads the precision for printed numerical values from the "Precision" parameter in the input file.
         * 
         * @param[out] precision Reference to an integer where the read precision value will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readPrecision(int& precision);

        /**
         * @brief Reads the method used to compute the reduced density matrix (RDM) from the "RDMMethod" parameter in the input file.
         * 
         * @param[out] rdmMethod Reference to an RDMMethod variable where the parsed method will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readRDMMethod(RDMMethod& rdmMethod);

        /**
         * @brief Reads the requested run type (job) from the "RunType" parameter in the input file.
         *
         * @param[out] runType Reference to a RunType variable where the parsed run type will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readRunType(RunType& runType);

        /**
         * @brief Reads the save pseudo-orbitals option from the "SavePseudoOrbitals" parameter in the input file.
         * 
         * @param[out] savePseudoOrbitals Reference to a boolean where the read option value will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readSavePseudoOrbitals(bool& savePseudoOrbitals);

        /**
         * @brief Reads the choice of saving the reduced density matrix from the "SaveReducedDensityMatrix" parameter in the input file.
         * 
         * @param[out] saveReducedDensityMatrix Reference to a boolean where the read option value will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readSaveReducedDensityMatrix(bool& saveReducedDensityMatrix);

        /**
         * @brief Reads the show progress option from the "ShowProgress" parameter in the input file.
         * 
         * @param[out] showProgress Reference to a boolean where the read option value will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readShowProgress(bool& showProgress);

        /**
         * @brief Reads the single charge option from the "SingleCharge" parameter in the input file.
         * 
         * @param[out] singleCharge Reference to a boolean where the read option value will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readSingleCharge(bool& singleCharge);

        /**
         * @brief Reads the selected grid size from the "Size" parameter in the input file.
         * 
         * If the selected grid size is "Custom", this function also reads the custom size parameters from the "CustomSizeData" parameter.
         *
         * @param[out] gridSize Reference to a GridSize variable where the parsed grid size will be stored.
         * @param[out] customSizeData Reference to a CustomSizeData variable where the custom size parameters will be stored (if applicable).
         * 
         * @return True if the "Size" parameter was successfully read, false otherwise.
         */
        bool readSize(GridSize& gridSize, CustomSizeData& customSizeData);

        /**
         * @brief Reads a list of spins from the "SpinList" parameter in the input file.
         *
         * @param[out] spinList Vector filled with parsed spin types.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readSpinList(std::vector<SpinType>& spinList);

        /**
         * @brief Reads the selected spin type from the "SpinType" parameter in the input file.
         *
         * @param[out] spinType Reference to a SpinType variable where the parsed spin type value will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readSpinType(SpinType& spinType);

        /**
         * @brief Reads the transition densities to consider from the "TransitionDensities" parameter in the input file.
         *
         * @param[out] transitionDensities Reference to a vector where the read transition densities will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readTransitionDensities(std::vector<std::array<int, 2>>& transitionDensities);

        /**
         * @brief Reads the name of the file that describes excited states transitions from the "TransitionsFile" parameter in the input file.
         *
         * @param[out] transitionsFileName Reference to a string where the read filename will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readTransitionsFileName(std::string& transitionsFileName);

        /**
         * @brief Reads the verbosity level from the "Verbose" parameter in the input file.
         * 
         * @param[out] verbose Reference to an integer where the read verbosity level will be stored.
         * 
         * @return True if the parameter was successfully read, false otherwise.
         */
        bool readVerbose(int& verbose);


        //----------------------------------------------------------------------------------------------------//
        // STATIC PROTECTED METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Builds a `Domain` suitable for cube creation from `orb` and sizing options.
         *
         * @param[in] orb Orbitals instance used to determine molecular extent.
         * @param[in] gridSize Grid size token (Coarse/Medium/Fine/Custom).
         * @param[in] customSizeData Custom size parameters (used when `gridSize==CUSTOM`).
         * @param[in] Nval Number of values per grid point (used by Domain::set_all).
         * @return Configured `Domain` instance.
         */
        static Domain buildDomainForCube(Orbitals& orb, const GridSize gridSize, const CustomSizeData& customSizeData, const int& Nval);

        /**
         * @brief Builds an Orbitals or Becke helper object from an analytic file.
         *
         * @tparam T Analytic file parser type (WFX, MOLDENGAB, FCHK, LOG, ...).
         * @tparam U Resulting class type (Orbitals or Becke).
         * @param[in] analyticFileName Name of the analytic file.
         * @return Constructed instance of type U, initialised with analytic data.
         */
        template<typename T, typename U>
        static U computeOrbitalsOrBecke(const std::string& analyticFileName);

        /**
         * @brief Detects the analytic file format and initialises an instance of the chosen class.
         *
         * @tparam T Resulting class type (Orbitals or Becke).
         * @param[out] analyticObject Constructed instance of type T, initialised with analytic data.
         * @param[in] analyticFileName Path to analytic file; detection by extension.
         */
        template <typename T>
        static void computeOrbitalsOrBecke(T& analyticObject, const std::string& analyticFileName);

        /**
         * @brief Creates and saves a grid file (.cube) from the passed Orbitals instance, over a defined domain.
         *
         * @param[in] orbitals Orbitals instance providing densities/orbitals.
         * @param[in] domain Domain describing grid geometry.
         * @param[in] cubeFileName Output filename for the cube.
         * @param[in] cubeType Type of cube file to create (density, orbitals, ELF).
         * @param[in] elfMethod ELF method selection (SAVIN/BECKE) when creating ELF.
         * @param[in] nums Orbital indices used for orbital grids.
         * @param[in] typesSpin Spin type for orbital grids.
         */
        static void createCube(Orbitals& orbitals, const Domain& domain, const std::string& cubeFileName, CubeType cubeType, bool showProgress = false, const ELFMethod elfMethod = ELFMethod::UNKNOWN, std::vector<int> nums = {0}, std::vector<SpinType> typesSpin = { SpinType::ALPHA });

        /**
         * @brief Selects all molecular orbitals for the requested spin configuration.
         *
         * @param[out] orbnums Output vector of orbital indices.
         * @param[out] orbspin Output vector of spin types corresponding to `orbnums`.
         * @param[in] o Orbitals instance for occupation info.
         * @param[in] spinType Requested spin selection.
         * @param[in] numberOfOrbitals Number of MOs available.
         */
        static void setAllOrbitals(std::vector<int> &orbnums, std::vector<SpinType> &orbspin, Orbitals &o, SpinType spinType, int numberOfOrbitals);

        /**
         * @brief Applies a custom orbital index list and corresponding spin list.
         *
         * @param[in,out] orbnums Input orbital indices (1-based expected; converted internally).
         * @param[out] orbspin Output vector filled with parsed spin types.
         * @param[in] spinList Input vector of SpinType values corresponding to custom orbitals.
         */
        static void setCustomOrbitals(std::vector<int>& orbnums, std::vector<SpinType>& orbspin, const std::vector<SpinType>& spinList);

        /**
         * @brief Selects occupied molecular orbitals according to occupations and spin selection.
         *
         * @param[out] orbnums Output vector of occupied orbital indices.
         * @param[out] orbspin Output vector of spin types for each selected orbital.
         * @param[in] o Orbitals instance containing occupation numbers.
         * @param[in] spinType Requested spin selection.
         * @param[in] N Number of MOs available.
         */
        static void setOccupiedOrbitals(std::vector<int>& orbnums, std::vector<SpinType>& orbspin, Orbitals& o, SpinType spinType, int N);

        /**
         * @brief Configures `orbitalsNumbers` and `orbitalsSpins` based on selection.
         *
         * @param[in] o Orbitals instance to query occupation/spin information.
         * @param[in] numberOfOrbitals Total number of molecular orbitals available.
         * @param[in,out] orbitalsNumbers Input/output list of orbital indices (modified in-place).
         * @param[in,out] orbitalsSpins Input/output list of spin flags (modified in-place).
         * @param[in] orbitalType Selection mode specifying which orbitals to include.
         * @param[in] spinType Spin selection (Alpha/Beta/Alpha-Beta).
         * @param[in] spinList Optional custom spin list for custom orbital selections.
         */
        static void setOrbitals(Orbitals& o, const int numberOfOrbitals, std::vector<int>& orbitalsNumbers, std::vector<SpinType>& orbitalsSpins, const OrbitalType orbitalType, SpinType spinType, const std::vector<SpinType>& spinList = {});

        /**
         * @brief Selects HOMO orbital(s) according to spin selection.
         *
         * @param[out] orbnums Output vector set to HOMO index(es).
         * @param[out] orbspin Output vector set to corresponding spin(s).
         * @param[in] o Orbitals instance used to query occupations.
         * @param[in] spinType Requested spin selection.
         */
        static void setHomo(std::vector<int>& orbnums, std::vector<SpinType>& orbspin, Orbitals& o, SpinType spinType);

        /**
         * @brief Selects HOMO and LUMO together according to spin selection.
         *
         * @param[out] orbnums Output vector containing HOMO and LUMO indices.
         * @param[out] orbspin Output vector of corresponding spin types.
         * @param[in] o Orbitals instance for occupation info.
         * @param[in] spinType Requested spin selection.
         */
        static void setHomoLumo(std::vector<int>& orbnums, std::vector<SpinType>& orbspin, Orbitals& o, SpinType spinType);

        /**
         * @brief Selects LUMO orbital(s) according to spin selection.
         *
         * @param[out] orbnums Output vector set to LUMO index(es).
         * @param[out] orbspin Output vector set to corresponding spin(s).
         * @param[in] o Orbitals instance used to query occupations.
         * @param[in] spinType Requested spin selection.
         */
        static void setLumo(std::vector<int>& orbnums, std::vector<SpinType>& orbspin, Orbitals& o, SpinType spinType);

        /**
         * @brief Selects virtual (unoccupied) molecular orbitals according to spin selection.
         *
         * @param[out] orbnums Output vector of virtual orbital indices.
         * @param[out] orbspin Output vector of spin types for each selected orbital.
         * @param[in] o Orbitals instance containing occupation numbers.
         * @param[in] spinType Requested spin selection.
         * @param[in] N Number of MOs available.
         */
        static void setVirtualOrbitals(std::vector<int>& orbnums, std::vector<SpinType>& orbspin, Orbitals& o, SpinType spinType, int N);


        //----------------------------------------------------------------------------------------------------//
        // OTHER PRIVATE METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * TODO
         */
        Orbitals computePseudoOrbitalsFromLrfMatrix(const Orbitals& orbitals, const std::vector<std::vector<std::vector<double>>>& lrfMatrix, std::vector<std::vector<double>>& eigenvalues, std::vector<std::vector<std::vector<double>>>& eigenvectors, const std::string& outputPrefix, bool savePseudoOrbitals, std::ostream& outputStream, int verbose, bool showProgress = false);

        /**
         * @brief Opens the configured input file.
         */
        void openInputFile();

        /**
         * @brief Prints critical points information from a Critical Points grid (GridCP). (Not implemented yet.)
         */
        void printCriticalPoints();


    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTORS AND DESTRUCTOR
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Default constructor.
         *
         * Sets `_inputFileName` to "input.txt", initialises job lists and
         * opens the input file stream.
         */
        Job();

        /**
         * @brief Construct a Job instance with a custom input filename.
         * @param inputFileName Path to the input file to use.
         */
        Job(const std::string& inputFileName);

        /**
         * @brief Destructor. Closes the input file stream.
         */
        ~Job();


        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Looks for the "runType" parameter in the input file and runs either the according job.
         * 
         * If the "runType" parameter is not found or if its value does not correspond to a valid job, this method will run the "Help" job.
         */
        virtual void run() = 0;
};

#include <JobControl/Job.tpp>

#endif /* CDFTT_JOB_H_INCLUDED */

