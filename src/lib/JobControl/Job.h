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

        /**
         * @brief Sets the list of available jobs and their descriptions.
         */
        void setJobList();


        //----------------------------------------------------------------------------------------------------//
        // SPECIFIC PARAMETERS READING FROM INPUT FILE
        //----------------------------------------------------------------------------------------------------//

        void readAnalyticFilesNames(std::vector<std::string>& analyticFilesNames);
        void readCutoff(double& cutoff);
        void readELFMethod(ELFMethod& elfMethod);
        void readEnergies(std::vector<double>& energies);
        void readGridFilesNames(std::vector<std::string>& gridFilesNames);
        void readOrbitalsNumbers(std::vector<int>& orbitalsNumbers);
        void readOrbitalsSpins(std::vector<SpinType>& orbitalsSpins);
        void readOrbitalType(OrbitalType& orbitalType);
        void readPartitionMethod(PartitionMethod& partitionMethod);
		void readRunType(RunType& runType);
        void readSize(GridSize& gridSize, CustomSizeData& customSizeData);
        void readSpinList(std::vector<SpinType>& spinList);
        void readSpinType(SpinType& spinType);
		void readTransitionsFileName(std::string& transitionsFileName);


        //----------------------------------------------------------------------------------------------------//
        // JOBS
        //----------------------------------------------------------------------------------------------------//

        void run_help();
        void run_computePartialCharges();
        void run_computeDescriptors();
        void run_computeIntegrals();
        void run_computeGridDifference();
		void run_lambdaDiagnostic();
        void run_makeDensityCube();
        void run_makeOrbitalsCube();
        void run_makeELFCube();
        void run_convertOrbitals();

        void buildBasins(GridCP& gcp, const std::string& GridFileName, PartitionMethod partitionMethod);

        //! BBS
        /*! takes a grid and a cube file and builds basin according to the sign of the values of the grid. if cutoff=0 two basins will be built. Otherwise, more will be created. cutoff sets the numerical limit of "zero" */
        void buildBasinsBySign(GridCP& gcp, const std::string& GridFileName,  double cutoff, bool two);
        
            //! Compute integrals
            /*! takes a critical point grid and integrates the read Grid Files' values over the domain of GridCP. prints the results*/
        void computeLocalIntegrals(GridCP& gcp, const std::vector<std::string>& GridFileNames);
            
            //! Print critical points
            /*! Should print the critical points of gridCP. Not yet implemented.*/
        void printCriticalPoints();

            //! Compute the difference of grids
            /*! computes GFN1-GFN2 for all values of the grids and saves the result in the chosen output name NameNew*/
        void ComputeGridDifference(const std::string& GFN1, const std::string& GFN2, const std::string& NameNew);

            //! Compute Partial Charges
            /*! computes the partial charges from the atoms of a grid by method of Atoms in molecule(i=0,1,2), VDD(i=3), Becke (i=4)*/ 
        std::vector<double> computePartialCharges(const std::string &gridfname, PartitionMethod partitionMethod);
            //! Compute Partial Charges
            /*! computes the partial charges of atoms from an analytical file*/
        std::vector<double> computePartialChargesAndEnergy(std::vector<double>& E, const std::string& ANAFileName);
            
            //! Open input file
            /*! Opens the input file. shows error message if it fails.*/
        void openInputFile();
            
            //! Print list of run types
            /*! Prints the list of jobs and their description. Used in RunType=Help*/ 
        void printListOfRunTypes();

            //! Read a line of input file
            /*! Reads one line of the input file and identifies the tag (eg: RunType) and reads the value of the tag following it. The value is stored into a string. This removes any preceding or trailing spaces*/
        //bool readOneString(const string& tag, string& value);

            //! Read one value of type T
            /*! Calls readOnestring to convert the string value to a value of type T. Removes any separators and replaces them with spaces*/
        //template<typename T> bool readOneType(const string& tag, T& x);

            //! Read a list of Type T
            /*! Calls readOnestring to convert the string value to a vector of values of type T. Removes any separators and replaces them with spaces*/
        //template<typename T> bool readListType(const string& tag, vector<T>& x);

            //! Compute Descriptors from grids
            /*! Computes chemical descriptors from 3 grids (electrophilic, nucleophilic, radical) and passing ionisation potential and electron affinity.*/
        Descriptors computeDescriptors(const std::string& GridFileName1, const std::string& GridFileName2, const std::string& GridFileName3, double I, double A, PartitionMethod partitionMethod);

            //! Compute descriptors from grids
            /*! Computes chemical descriptors from 3 grids(electrophilic, nucleophilic, radical) and passes total energies of the 3 files*/
        Descriptors computeDescriptors(const std::string &GridFileName1, const std::string &GridFileName2, const std::string &GridFileName3, std::vector<double> E, PartitionMethod partitionMethod);
            //! return Structure
            /*! returns the Structure contained in an analytical class*/
        Structure returnStruct(const std::string& ANAFileName);

            //! Construct Orbitals or becke object
            /*! Constructs Orbitals or Becke (U) from template arguments WFX, MOLDENGAB, FCHK, LOG (T) from analytic file wfx, molden, gab, log, fchk*/
        template<typename T,typename U> U computeOrbOrBecke(const std::string& analyticFileName);
            //! Construct Orbitals or Becke
            /*! Initializes Orbitals or Becke conditionally based on what the detected file type */
        template<typename T> void readFileFormat(T& AnaClass, const std::string& FileName);

            //! Create and save a grid onto a cube file
            /*! Takes a Domain and creates a Grid and saves it in cubeFileName. TypeFlag determines the type of value the grid will contain. Additional parameters may be required for various density types*/
        void createCube(Orbitals& orb, const Domain& d, const std::string& cubeFileName, int TypeFlag, const ELFMethod elfMethod = ELFMethod::SAVIN, std::vector<int> nums={0}, std::vector<int> typesSpin={0});
            //! Create a Domain from Orbitals
            /*! Initializes a Domain. Sizes can be coarse fine medium or custom. Custom requires csizes to be specified.*/ 
        Domain buildDomainForCube(Orbitals& orb, const GridSize gridSize,const CustomSizeData& customSizeData, const int& Nval);


        void setOrbitals(Orbitals& o, const int numberOfOrbitals, std::vector<int>& orbitalsNumbers, std::vector<int>& orbitalsSpins, const OrbitalType orbitalType, SpinType spinType, const std::vector<SpinType>& spinList = {});

            //! Modifies a vector
            /*! Modifies the vectors orbnums and orbspin appropriately to include all Molecular Orbitals for Orbitals::phis()*/
        void setAllOrb(std::vector<int> &orbnums, std::vector<int> &orbspin, Orbitals &o, SpinType spinType, const int& numberOfOrbitals);
            //! Modifies a vector
            /*! Modifies the vectors orbnums and orbspin appropriately to include occupied Molecular Orbitals for Orbitals::phis()*/
        void setOccOrb(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, SpinType spinType, const int& N);
        
            //! Modifies a vector
            /*! Modifies the vectors orbnums and orbspin appropriately to include virtual Molecular Orbitals for Orbitals::phis()*/
        void setVirtOrb(std::vector<int> &orbnums, std::vector<int> &orbspin, Orbitals &o, SpinType spinType, const int &N);

            //! Modifies a vector
            /*! Modifies the vectors orbnums and orbspin appropriately to include LUMO Molecular Orbitals for Orbitals::phis()*/
        void setLumo(std::vector<int> &orbnums, std::vector<int> &orbspin, Orbitals &o, SpinType spinType);

            //! Modifies a vector
            /*! Modifies the vectors orbnums and orbspin appropriately to include HOMO Molecular Orbitals for Orbitals::phis()*/
        void setHomo(std::vector<int>& orbnums, std::vector<int>& orbspin, Orbitals& o, SpinType spinType);
        //! Modifies a vector
        /*! Modifies the vectors orbnums and orbspin appropriately to include HOMO-LUMO Molecular Orbitals for Orbitals::phis()*/
        void setHomoLumo(std::vector<int> &orbnums, std::vector<int> &orbspin, Orbitals &o, SpinType spinType);

            //! Modifies a vector
            /*! Modifies the vectors orbnums and orbspin appropriately to include a custom list of Molecular Orbitals for Orbitals::phis()*/
        void setCustom(std::vector<int>& orbnums, std::vector<int>& orbspin, const std::vector<SpinType>& spinList);
            //! Compute descriptors 
            /*! compute descriptors from 3 analytical files following the finite difference method of Becke*/
        void computeDescriptorsFD(const std::string& ANAFileName1, const std::string& ANAFileName2, const std::string& ANAFileName3);


    public:
            //!Default constructor
            /*! set _inputFileName to input.cdft and open it*/
        Job();

            //! constructor
            /*! parameter as input file name and open it */
        Job(std::string inputFileName);

            //! destructor
            /*! close input file */
        ~Job();

            //! method
            /*! run job */
        void run();
};

#endif /* CDFTT_JOB_H_INCLUDED */

