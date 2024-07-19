#ifndef __CDFTT_JOB_H_INCLUDED__
#define __CDFTT_JOB_H_INCLUDED__
#include<Common/Descriptors.h>
#include<string>
#include<fstream>
#include<vector>
#include <Orbitals/Orbitals.h>
#include <Common/PeriodicTable.h>

class Job
{
	private:
		PeriodicTable _table;
		std::string _inputFileName;
		std::ifstream _inputFile;
		std::vector<std::string> _jobsList;
		std::vector<std::string> _jobDescription;
		void setJobList();
		void buildBasins(GridCP& gcp, const string& GridFileName, int method);
		void buildBasinsBySign(GridCP& gcp, const string& GridFileName,  double cutoff, bool two);
		void computeLocalIntegrals(GridCP& gcp, const vector<string>& GridFileNames);
		void printCriticalPoints();
		void ComputeGridDifference(const string& GFN1, const string& GFN2, const string& NameNew);
		std::vector<double> computePartialCharges(const string& gridfname, int method);
		void openInputFile();
		void printListOfRunTypes();
		bool readOneString(const string& tag, string& value);
		template<typename T> bool readOneType(const string& tag, T& x);
		template<typename T> bool readListType(const string& tag, vector<T>& x);
		Descriptors computeDescriptors(const string& GridFileName1, const string& GridFileName2, const string& GridFileName3, double I, double A, int AIMmethod);
		Descriptors computeDescriptors(const string& GridFileName1, const string& GridFileName2, const string& GridFileName3, vector<double>E,  int AIMmethod);
		template<typename T> Orbitals computeOrbitals(const string& analyticFileName);
		void createCube(Orbitals& orb, const Domain& d, const string& CubeFileName, bool orbFlag, vector<int>, vector<int>);
		Domain DomainForCube(Orbitals& orb, const string& size,const vector<double>& csizes, const int& Nval);
		void setAllOrb(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, bool A, bool B, bool AnB, const int& N);
		void setOccOrb(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, bool A, bool B, bool AnB, const int& N);
		void setVirtOrb(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, bool A, bool B, bool AnB, const int& N);
		void setLumo(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, bool A, bool B, bool AnB);
		void setHomo(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, bool A, bool B, bool AnB);
		void setHomoLumo(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, bool A, bool B, bool AnB);
		void setCustom(vector<int>& orbnums, vector<int>& orbspin, vector<string>& tmplist);
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

