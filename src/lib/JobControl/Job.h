#ifndef __CDFTT_JOB_H_INCLUDED__
#define __CDFTT_JOB_H_INCLUDED__

#include<string>
#include<fstream>
#include<vector>
#include <Common/PeriodicTable.h>

class Job
{
	private:
		PeriodicTable _table;
		std::string _inputFileName;
		std::ifstream _inputFile;
		std::vector<std::string> _jobsList;
		void setJobList();
		void buildBasins();
		void computeLocalIntegrals();
		void printCriticalPoints();
		std::vector<double> computeAIMCharges(const string& gridfname, int method);
		void openInputFile();
		void printListOfRunTypes();
		bool readOneString(const string& tag, string& value);
		template<typename T> bool readOneType(const string& tag, T& x);
		template<typename T> bool readListType(const string& tag, vector<T>& x);
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

