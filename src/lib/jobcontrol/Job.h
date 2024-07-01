#ifndef __CDFTT_JOB_H_INCLUDED__
#define __CDFTT_JOB_H_INCLUDED__

#include<string>
#include<fstream>
#include<vector>
#include <common/PeriodicTable.h>

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
		std::vector<double> computeAIMCharges(const string& gridfname);
		void openInputFile();
		void printListOfRunTypes();
		bool readOneString(const string& tag, string& value);
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

