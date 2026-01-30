#ifndef CDFTT_JOB_H_INCLUDED
#define CDFTT_JOB_H_INCLUDED

#include <string>
#include <fstream>
#include <set>
#include <vector>

#include <Common/Descriptors.h>
#include <Common/PeriodicTable.h>
#include <Cube/GridCP.h>
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
    private:
        /** @brief Periodic table used for element lookups. */
        PeriodicTable _table;

        /** @brief Name of the input file to parse. */
        std::string _inputFileName;

        /** @brief Input file stream opened on the input file. */
        std::ifstream _inputFile;

        /** @brief List of available job names (displayed in help). */
        std::vector<std::string> _jobsList;

        /** @brief Descriptions of the available jobs. */
        std::vector<std::string> _jobDescription;


        //----------------------------------------------------------------------------------------------------//
        // SPECIFIC PARAMETERS READING FROM INPUT FILE
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Reads the name(s) of the analytic file(s) from the "AnalyticFiles" parameter in the input file.
         *
         * @param[out] analyticFilesNames Reference to a vector where the read filename(s) will be stored.
         */
        void readAnalyticFilesNames(std::vector<std::string>& analyticFilesNames);

        /** @brief Reads the charge of the point charge from the "Charge" parameter in the input file.
         * 
         * @param[out] charge Reference to a double where the read charge value will be stored.
         */
        void readCharge(double &charge);

        /** @brief Reads the numeric cutoff used by some partitioning methods from the "Cutoff" parameter in the input file.
         * 
         * @param[out] cutoff Reference to a double where the read cutoff value will be stored.
         */
        void readCutoff(double& cutoff);

        /**
         * @brief Reads the method used to compute the Electron Localization Function (ELF) from the "ELFMethod" parameter in the input file.
         *
         * @param[out] elfMethod Reference to an ELFMethod variable where the parsed method will be stored.
         */
        void readELFMethod(ELFMethod& elfMethod);

        /**
         * @brief Reads the energies from the "Energies" parameter in the input file.
         *
         * @param[out] energies Reference to a vector where the energy values will be stored.
         */
        void readEnergies(std::vector<double>& energies);

        /**
         * @brief Reads the name(s) of the grid file(s) (.cube files) from the "GridFiles" parameter in the input file.
         *
         * @param[out] gridFilesNames Reference to a vector where the read filename(s) will be stored.
         */
        void readGridFilesNames(std::vector<std::string>& gridFilesNames);

        /**
         * @brief Reads the cutoff distance for nuclear contribution from the "NuclearCutoff" parameter in the input file.
         *
         * @param[out] nuclearCutoff Reference to a double where the read cutoff distance will be stored.
         */
        void readNuclearCutoff(double& nuclearCutoff);

        /**
         * @brief Reads a list of orbital numbers from the "OrbitalsNumbers" parameter in the input file.
         *
         * @param[out] orbitalsNumbers Reference to a vector where the orbital numbers will be stored.
         */
        void readOrbitalsNumbers(std::vector<int>& orbitalsNumbers);

        /**
         * @brief Reads a list of orbitals spins from the "OrbitalsSpins" parameter in the input file.
         *
         * @param[out] orbitalsSpins Reference to a vector where the parsed spin types will be stored.
         */
        void readOrbitalsSpins(std::vector<SpinType>& orbitalsSpins);

        /**
         * @brief Reads the selected orbital type from the "OrbitalType" parameter in the input file.
         *
         * @param[out] orbitalType Reference to an OrbitalType variable where the parsed selected orbital type will be stored.
         */
        void readOrbitalType(OrbitalType& orbitalType);

        /**
         * @brief Reads the partition method from the "PartitionMethod" parameter in the input file.
         *
         * @param[out] partitionMethod Reference to a PartitionMethod variable where the parsed method will be stored.
         */
        void readPartitionMethod(PartitionMethod& partitionMethod);

        /**
         * @brief Reads a 3-coordinate position from the "Position" parameter in the input file.
         *
         * @param[out] position Reference to an array of three doubles where the read position will be stored.
         */
        void readPosition(std::array<double, 3>& position);

        /**
         * @brief Reads the requested run type (job) from the "RunType" parameter in the input file.
         *
         * @param[out] runType Reference to a RunType variable where the parsed run type will be stored.
         */
        void readRunType(RunType& runType);

        /**
         * @brief Reads the selected grid size from the "GridSize" parameter in the input file.
         * 
         * If the selected grid size is "Custom", this function also reads the custom size parameters from the "CustomSizeData" parameter.
         *
         * @param[out] gridSize Reference to a GridSize variable where the parsed grid size will be stored.
         * @param[out] customSizeData Reference to a CustomSizeData variable where the custom size parameters will be stored (if applicable).
         */
        void readSize(GridSize& gridSize, CustomSizeData& customSizeData);

        /**
         * @brief Reads a list of spins (`SpinList`) from input.
         *
         * @param[out] spinList Vector filled with parsed spin types.
         */
        void readSpinList(std::vector<SpinType>& spinList);

        /**
         * @brief Reads the selected spin type from the "SpinType" parameter in the input file.
         *
         * @param[out] spinType Reference to a SpinType variable where the parsed spin type value will be stored.
         */
        void readSpinType(SpinType& spinType);

        /**
         * @brief Reads the name of the file that describes excited states transitions from the "TransitionsFile" parameter in the input file.
         *
         * @param[out] transitionsFileName Reference to a string where the read filename will be stored.
         */
        void readTransitionsFileName(std::string& transitionsFileName);


        //----------------------------------------------------------------------------------------------------//
        // JOBS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Sets the list of available jobs and their descriptions.
         */
        void setJobList();

        /**
         * @brief Runs the job associated with the "Help" runtype.
         */
        void run_help();

        /**
         * @brief Runs the job associated with the "ComputeDescriptors" runtype.
         */
        void run_computeDescriptors();

        /**
         * @brief Runs the job associated with the "ComputeEnergyWithPointCharge" runtype.
         */
        void run_computeEnergyWithPointCharge();

        /**
         * @brief Runs the job associated with the "ComputeGridDifference" runtype.
         */

        void run_computeGridDifference();

        /**
         * @brief Runs the job associated with the "ComputeIntegrals" runtype.
         */
        void run_computeIntegrals();

        /**
         * @brief Runs the job associated with the "ComputePartialCharges" runtype.
         */
        void run_computePartialCharges();

        /**
         * @brief Runs the job associated with the "LambdaDiagnostic" runtype.
         */
        void run_lambdaDiagnostic();

        /**
         * @brief Runs the job associated with the "MakeDensityCube" runtype.
         */
        void run_makeDensityCube();

        /**
         * @brief Runs the job associated with the "MakeOrbitalsCube" runtype.
         */
        void run_makeOrbitalsCube();

        /**
         * @brief Runs the job associated with the "MakeELFCube" runtype.
         */
        void run_makeELFCube();

        /**
         * @brief Runs the job associated with the "ConvertOrbitals" runtype.
         */
        void run_convertOrbitals();


        /**
         * @brief Builds basins from a grid using `partitionMethod`.
         *
         * @param[out] gcp GridCP reference that will be populated with the constructed basins.
         * @param[in] GridFileName Path to the input grid (.cube file) used to build basins.
         * @param[in] partitionMethod Selected partition method.
         */
        void buildBasins(GridCP& gcp, const std::string& GridFileName, PartitionMethod partitionMethod);

        /**
         * @brief Builds basins by sign of grid values (BBS/B2S partition methods).
         *
         * @param[out] gcp GridCP reference that will be populated with sign-based basins.
         * @param[in] GridFileName Path to the grid file used as input.
         * @param[in] cutoff Numerical cutoff below which values are considered zero.
         * @param[in] two If true, builds exactly two basins.
         */
        void buildBasinsBySign(GridCP& gcp, const std::string& GridFileName,  double cutoff, bool two);
        
        /**
         * @brief Integrates the provided grids over the domain of a Critical Points grid (GridCP).
         *
         * @param[in,out] gcp GridCP reference that defines the integration domain and stores results.
         * @param[in] GridFileNames Filenames of grids to integrate.
         */
        void computeLocalIntegrals(GridCP& gcp, const std::vector<std::string>& GridFileNames);

        /**
         * @brief Prints critical points information from a Critical Points grid (GridCP). (Not implemented yet.)
         */
        void printCriticalPoints();

        /**
         * @brief Computes the difference between two grids and saves the result to an output file.
         *
         * @param[in] minuendGridFilename Left (minuend) input grid filename.
         * @param[in] subtrahendGridFilename Right (subtrahend) input grid filename.
         * @param[in] outputGridFilename Output filename for the difference grid.
         */
        void ComputeGridDifference(const std::string& minuendGridFilename, const std::string& subtrahendGridFilename, const std::string& outputGridFilename);

        /**
         * @brief Computes partial atomic charges from a grid using the specified method.
         *
         * @param[in] gridFilename Input grid filename used for partitioning and integration.
         * @param[in] partitionMethod Partition method to use (AIM, VDD, Becke, ...).
         * @return Vector of computed partial charges.
         */
        std::vector<double> computePartialCharges(const std::string& gridFilename, PartitionMethod partitionMethod);

        /**
         * @brief Computes partial charges and energy values from an analytic file.
         *
         * @param[out] energies Reference to vector where energy results will be stored.
         * @param[in] analyticFileName Analytic file path.
         * @return Vector of computed partial charges.
         */
        std::vector<double> computePartialChargesAndEnergy(std::vector<double>& energies, const std::string& analyticFileName);
            
        /**
         * @brief Opens the configured input file.
         */
        void openInputFile();
            
        /**
         * @brief Prints the list of available run types and their descriptions.
         */
        void printListOfRunTypes();

        /**
         * @brief Computes descriptors from three grid files using ionization energy and electron affinity.
         *
         * @param[in] gridFileName1 Name of the electrophilic grid file.
         * @param[in] gridFileName2 Name of the nucleophilic grid file.
         * @param[in] gridFileName3 Name of the radical grid file.
         * @param[in] ionizationEnergy Ionization energy.
         * @param[in] electronAffinity Electron affinity.
         * @param[in] partitionMethod Selected partition method.
         * @return Descriptors object containing computed descriptors.
         */
        Descriptors computeDescriptors(const std::string& gridFileName1, const std::string& gridFileName2, const std::string& gridFileName3, double ionizationEnergy, double electronAffinity, PartitionMethod partitionMethod);

        /**
         * @brief Computes descriptors from three grid files using these file energies.
         *
         * @param[in] gridFileName1 Name of the first grid file.
         * @param[in] gridFileName2 Name of the second grid file.
         * @param[in] gridFileName3 Name of the third grid file.
         * @param[in] energies Vector of energies respectively corresponding to the files.
         * @param[in] partitionMethod Selected partition method.
         * @return Descriptors object containing computed descriptors.
         */
        Descriptors computeDescriptors(const std::string& gridFileName1, const std::string& gridFileName2, const std::string& gridFileName3, const std::vector<double>& energies, PartitionMethod partitionMethod);

        /**
         * @brief Extracts the molecular structure from an analytic file.
         *
         * @param[in] analyticFileName Name of the analytic file.
         * @return Parsed Structure instance.
         */
        Structure returnStruct(const std::string& analyticFileName);

        /**
         * @brief Builds an Orbitals or Becke helper object from an analytic file.
         *
         * @tparam T Analytic file parser type (WFX, MOLDENGAB, FCHK, LOG, ...).
         * @tparam U Resulting class type (Orbitals or Becke).
         * @param[in] analyticFileName Name of the analytic file.
         * @return Constructed instance of type U, initialised with analytic data.
         */
        template<typename T, typename U> U computeOrbitalsOrBecke(const std::string& analyticFileName);

        /**
         * @brief Detects the analytic file format and initialises an instance of the chosen class.
         *
         * @tparam T Resulting class type (Orbitals or Becke).
         * @param[out] analyticObject Constructed instance of type T, initialised with analytic data.
         * @param[in] analyticFileName Path to analytic file; detection by extension.
         */
        template<typename T> void computeOrbitalsOrBecke(T& analyticObject, const std::string& analyticFileName);

        /**
         * @brief Creates and saves a grid file (.cube) from the passed Orbitals instance, over a defined domain.
         *
         * @param[in] orbitals Orbitals instance providing densities/orbitals.
         * @param[in] domain Domain describing grid geometry.
         * @param[in] cubeFileName Output filename for the cube.
         * @param[in] TypeFlag 0=density, 1=orbitals, else ELF.
         * @param[in] elfMethod ELF method selection (SAVIN/BECKE) when creating ELF.
         * @param[in] nums Orbital indices used for orbital grids.
         * @param[in] typesSpin Spin flags for orbital grids.
         */
        void createCube(Orbitals& orbitals, const Domain& domain, const std::string& cubeFileName, int TypeFlag, const ELFMethod elfMethod = ELFMethod::SAVIN, std::vector<int> nums = {0}, std::vector<int> typesSpin = {0});

        /**
         * @brief Builds a `Domain` suitable for cube creation from `orb` and sizing options.
         *
         * @param[in] orb Orbitals instance used to determine molecular extent.
         * @param[in] gridSize Grid size token (Coarse/Medium/Fine/Custom).
         * @param[in] customSizeData Custom size parameters (used when `gridSize==CUSTOM`).
         * @param[in] Nval Number of values per grid point (used by Domain::set_all).
         * @return Configured `Domain` instance.
         */
        Domain buildDomainForCube(Orbitals& orb, const GridSize gridSize, const CustomSizeData& customSizeData, const int& Nval);


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
        void setOrbitals(Orbitals& o, const int numberOfOrbitals, std::vector<int>& orbitalsNumbers, std::vector<int>& orbitalsSpins, const OrbitalType orbitalType, SpinType spinType, const std::vector<SpinType>& spinList = {});

        /**
         * @brief Selects all molecular orbitals for the requested spin configuration.
         *
         * @param[out] orbnums Output vector of orbital indices.
         * @param[out] orbspin Output vector of spin flags corresponding to `orbnums`.
         * @param[in] o Orbitals instance for occupation info.
         * @param[in] spinType Requested spin selection.
         * @param[in] numberOfOrbitals Number of MOs available.
         */
        void setAllOrb(std::vector<int> &orbnums, std::vector<int> &orbspin, Orbitals &o, SpinType spinType, const int& numberOfOrbitals);

              /**
               * @brief Selects occupied molecular orbitals according to occupations and spin selection.
               *
               * @param[out] orbnums Output vector of occupied orbital indices.
               * @param[out] orbspin Output vector of spin flags for each selected orbital.
               * @param[in] o Orbitals instance containing occupation numbers.
               * @param[in] spinType Requested spin selection.
               * @param[in] N Number of MOs available.
               */
           void setOccOrb(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, SpinType spinType, const int& N);
        
              /**
               * @brief Selects virtual (unoccupied) molecular orbitals according to spin selection.
               *
               * @param[out] orbnums Output vector of virtual orbital indices.
               * @param[out] orbspin Output vector of spin flags for each selected orbital.
               * @param[in] o Orbitals instance containing occupation numbers.
               * @param[in] spinType Requested spin selection.
               * @param[in] N Number of MOs available.
               */
           void setVirtOrb(std::vector<int> &orbnums, std::vector<int> &orbspin, Orbitals &o, SpinType spinType, const int &N);

              /**
               * @brief Selects LUMO orbital(s) according to spin selection.
               *
               * @param[out] orbnums Output vector set to LUMO index(es).
               * @param[out] orbspin Output vector set to corresponding spin(s).
               * @param[in] o Orbitals instance used to query occupations.
               * @param[in] spinType Requested spin selection.
               */
           void setLumo(std::vector<int> &orbnums, std::vector<int> &orbspin, Orbitals &o, SpinType spinType);

              /**
               * @brief Selects HOMO orbital(s) according to spin selection.
               *
               * @param[out] orbnums Output vector set to HOMO index(es).
               * @param[out] orbspin Output vector set to corresponding spin(s).
               * @param[in] o Orbitals instance used to query occupations.
               * @param[in] spinType Requested spin selection.
               */
           void setHomo(std::vector<int>& orbnums, std::vector<int>& orbspin, Orbitals& o, SpinType spinType);

        /**
         * @brief Selects HOMO and LUMO together according to spin selection.
         *
         * @param[out] orbnums Output vector containing HOMO and LUMO indices.
         * @param[out] orbspin Output vector of corresponding spin flags.
         * @param[in] o Orbitals instance for occupation info.
         * @param[in] spinType Requested spin selection.
         */
        void setHomoLumo(std::vector<int> &orbnums, std::vector<int> &orbspin, Orbitals &o, SpinType spinType);

              /**
               * @brief Applies a custom orbital index list and corresponding spin list.
               *
               * @param[in,out] orbnums Input orbital indices (1-based expected; converted internally).
               * @param[out] orbspin Output vector filled with parsed spin flags.
               * @param[in] spinList Input vector of SpinType values corresponding to custom orbitals.
               */
           void setCustom(std::vector<int>& orbnums, std::vector<int>& orbspin, const std::vector<SpinType>& spinList);

              /**
               * @brief Computes descriptors with the finite-difference Becke method from analytic files.
               *
               * @param[in] ANAFileName1 First analytic input file.
               * @param[in] ANAFileName2 Second analytic input file.
               * @param[in] ANAFileName3 Third analytic input file.
               */
           void computeDescriptorsFD(const std::string& ANAFileName1, const std::string& ANAFileName2, const std::string& ANAFileName3);


    public:
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
        Job(std::string inputFileName);

            /**
             * @brief Destructor closes the input file stream.
             */
        ~Job();

            /**
             * @brief Parse input and run the selected job.
             */
        void run();
};

#endif /* CDFTT_JOB_H_INCLUDED */

