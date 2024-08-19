#ifndef __CDFTT_JOB_H_INCLUDED__
#define __CDFTT_JOB_H_INCLUDED__
#include<Common/Descriptors.h>
#include<string>
#include<fstream>
#include<vector>
#include <Orbitals/Orbitals.h>
#include <Common/PeriodicTable.h>

	//! Job class
	/*! Class that manages the user-available jobs to be run by CDFTT*/
class Job
{
	private:
		PeriodicTable _table;
		std::string _inputFileName;
		std::ifstream _inputFile;
		std::vector<std::string> _jobsList;
		std::vector<std::string> _jobDescription;

			//! Set Job list
			/*! sets the list of possible jobs and the description of these jobs*/
		void setJobList();

			//! Build integrable basins
			/*! takes a grid and a cube file and builds basins (integrable sub domains) according to the method specified( AIM, Becke, ...)*/
		void buildBasins(GridCP& gcp, const string& GridFileName, int method);

			//! BBS
			/*! takes a grid and a cube file and builds basin according to the sign of the values of the grid. if cutoff=0 two basins will be built. Otherwise, more will be created. cutoff sets the numerical limit of "zero" */
		void buildBasinsBySign(GridCP& gcp, const string& GridFileName,  double cutoff, bool two);
		
			//! Compute integrals
			/*! takes a critical point grid and integrates the read Grid Files' values over the domain of GridCP. prints the results*/
		void computeLocalIntegrals(GridCP& gcp, const vector<string>& GridFileNames);
			
			//! Print critical points
			/*! Should print the critical points of gridCP. Not yet implemented.*/
		void printCriticalPoints();

			//! Compute the difference of grids
			/*! computes GFN1-GFN2 for all values of the grids and saves the result in the chosen output name NameNew*/
		void ComputeGridDifference(const string& GFN1, const string& GFN2, const string& NameNew);

			//! Compute Partial Charges
			/*! computes the partial charges from the atoms of a grid by method of Atoms in molecule(i=0,1,2), VDD(i=3), Becke (i=4)*/ 
		std::vector<double> computePartialCharges(const string& gridfname, int method);

			//! Compute Partial Charges
			/*! computes the partial charges of atoms from an analytical file*/
		vector<double> computePartialChargesAndEnergy(vector<double>& E, const string& ANAFileName);
			
			//! Open input file
			/*! Opens the input file. shows error message if it fails.*/
		void openInputFile();
			
			//! Print list of run types
			/*! Prints the list of jobs and their description. Used in RunType=Help*/ 
		void printListOfRunTypes();

			//! Read a line of input file
			/*! Reads one line of the input file and identifies the tag (eg: RunType) and reads the value of the tag following it. The value is stored into a string. This removes any preceding or trailing spaces*/
		bool readOneString(const string& tag, string& value);

			//! Read one value of type T
			/*! Calls readOnestring to convert the string value to a value of type T. Removes any separators and replaces them with spaces*/
		template<typename T> bool readOneType(const string& tag, T& x);

			//! Read a list of Type T
			/*! Calls readOnestring to convert the string value to a vector of values of type T. Removes any separators and replaces them with spaces*/
		template<typename T> bool readListType(const string& tag, vector<T>& x);

			//! Compute Descriptors from grids
			/*! Computes chemical descriptors from 3 grids (electrophilic, nucleophilic, radical) and passing ionisation potential and electron affinity.*/
		Descriptors computeDescriptors(const string& GridFileName1, const string& GridFileName2, const string& GridFileName3, double I, double A, int AIMmethod);

			//! Compute descriptors from grids
			/*! Computes chemical descriptors from 3 grids(electrophilic, nucleophilic, radical) and passes total energies of the 3 files*/
		Descriptors computeDescriptors(const string& GridFileName1, const string& GridFileName2, const string& GridFileName3, vector<double>E,  int AIMmethod);

			//! return Structure
			/*! returns the Structure contained in an analytical class*/
		Structure returnStruct(const string& ANAFileName);

			//! Construct Orbitals or becke object
			/*! Constructs Orbitals or Becke (U) from template arguments WFX, MOLDENGAB, FCHK, LOG (T) from analytic file wfx, molden, gab, log, fchk*/
		template<typename T,typename U> U computeOrbOrBecke(const string& analyticFileName);

			//! Construct Orbitals or Becke
			/*! Initializes Orbitals or Becke conditionally based on what the detected file type */
		template<typename T> void readFileFormat(T& AnaClass, const string& FileName);

			//! Create and save a grid onto a cube file
			/*! Takes a Domain and creates a Grid and saves it in cubeFileName. TypeFlag determines the type of value the grid will contain. Additional parameters may be required for various density types*/
		void createCube(Orbitals& orb, const Domain& d, const string& cubeFileName, int TypeFlag, const string& ELFtype ="", vector<int> nums={0}, vector<int> typesSpin={0});

			//! Create a Domain from Orbitals
			/*! Initializes a Domain. Sizes can be coarse fine medium or custom. Custom requires csizes to be specified.*/ 
		Domain DomainForCube(Orbitals& orb, const string& size,const vector<double>& csizes, const int& Nval);

			//! Modifies a vector
			/*! Modifies the vectors orbnums and orbspin appropriately to include all Molecular Orbitals for Orbitals::phis()*/
		void setAllOrb(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, bool A, bool B, bool AnB, const int& N);

			//! Modifies a vector
			/*! Modifies the vectors orbnums and orbspin appropriately to include occupied Molecular Orbitals for Orbitals::phis()*/
		void setOccOrb(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, bool A, bool B, bool AnB, const int& N);
		
			//! Modifies a vector
			/*! Modifies the vectors orbnums and orbspin appropriately to include virtual Molecular Orbitals for Orbitals::phis()*/
		void setVirtOrb(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, bool A, bool B, bool AnB, const int& N);

			//! Modifies a vector
			/*! Modifies the vectors orbnums and orbspin appropriately to include LUMO Molecular Orbitals for Orbitals::phis()*/
		void setLumo(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, bool A, bool B, bool AnB);

			//! Modifies a vector
			/*! Modifies the vectors orbnums and orbspin appropriately to include HOMO Molecular Orbitals for Orbitals::phis()*/
		void setHomo(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, bool A, bool B, bool AnB);

			//! Modifies a vector
			/*! Modifies the vectors orbnums and orbspin appropriately to include HOMO-LUMO Molecular Orbitals for Orbitals::phis()*/
		void setHomoLumo(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, bool A, bool B, bool AnB);

			//! Modifies a vector
			/*! Modifies the vectors orbnums and orbspin appropriately to include a custom list of Molecular Orbitals for Orbitals::phis()*/
		void setCustom(vector<int>& orbnums, vector<int>& orbspin, vector<string>& tmplist);

			//! Compute descriptors 
			/*! compute descriptors from 3 analytical files following the finite difference method of Becke*/
		void computeDescriptorsFD(const string& ANAFileName1, const string& ANAFileName2, const string& ANAFileName3);
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

#endif /* __CDFTT_JOB_H_INCLUDED__ */

