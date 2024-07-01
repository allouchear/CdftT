using namespace std;
#include <jobcontrol/Job.h>
#include <numeric/GridCP.h>
#include <numeric/functions.h>
#include <numeric/Grid.h>
#include <common/Element.h>
#include <common/Descriptors.h>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <regex>
#ifdef ENABLE_OMP
#include <omp.h>
#endif

/******************************************************************************************/
std::string to_lower(const std::string& str)
{
	std::string res=str;
	for(size_t i=0;i<res.length();i++)
		res[i] = tolower(res[i]);
    return res;
}
std::string to_upper(const std::string& str)
{
	std::string res=str;
	for(size_t i=0;i<res.length();i++)
		res[i] = toupper(res[i]);
    return res;
}
/*************************************************************************************/
bool Job::readOneString(const string& tag, string& value)
{
	if(tag.length()<1) return false;

	string TAG = to_upper(tag);
	_inputFile.clear();
	_inputFile.seekg(0);

	string t;
	string t2;
	value="";
	while(!_inputFile.eof())
  	{
    		getline(_inputFile,t);
    		if(_inputFile.fail()) break;
		t = std::regex_replace(t, std::regex("^ +"), ""); // Easy removing leading,
		//t = std::regex_replace(t, std::regex(" +$"), ""); // removing only trailing spaces
		//t = regex_replace(t, std::regex(" +"), " "); // removing only extra spaces
		//t = std::regex_replace(t, std::regex("^ +| +$|( ) +"), "$1"); // Easy removing leading, trailing and extra spaces from a std::string in one line
		if(t[0]=='#') continue;
		t2=to_upper(t);
		std::size_t pos =  t2.find(TAG);
		if (pos == std::string::npos) continue;
		if(t2.find("=")!= std::string::npos) 
		{
			
			pos =  t2.find("=");
			value = t.substr(pos+1);
		}
		else 
		{
			pos =  t2.find(" ");
			value = t.substr(pos+1);
		}
		if(value.length()>0)
		{
			value = std::regex_replace(value, std::regex(" +"), " ");
			value = std::regex_replace(value, std::regex("^ +"), "");
			return true;
		}
	}
	return false;
}
/******************************************************************************************/
void Job::setJobList()
{
	_jobsList = {"computeAIMCharges" ,"Help"}; 
}
void Job::printListOfRunTypes()
{
	cout<<"available runType :"<<endl;
	for(size_t i=0;i<_jobsList.size();i++) 
			cout<<_jobsList[i]<<endl;
}
/******************************************************************************************/
void Job::openInputFile()
{
	_inputFile.open(_inputFileName);
	if(_inputFile.fail())
	{
		cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		cerr<<"Sorry, I cannot open the input file : "<<_inputFileName<<endl;
		cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		exit(1);
	}
}
/******************************************************************************************/
Job::Job():_inputFileName("input.cdft")
{
	setJobList();
	openInputFile();
}
/******************************************************************************************/
Job::Job(string inputFileName):_inputFileName(inputFileName)
{
	setJobList();
	_inputFile.open(_inputFileName);
	openInputFile();
}
/******************************************************************************************/
Job::~Job()
{
	_inputFile.close();
}
/******************************************************************************************/
void Job::buildBasins()
{
}
void Job::computeLocalIntegrals()
{
}
void Job::printCriticalPoints()
{
}
vector<double> Job::computeAIMCharges(const string& gridfname)
{
	ifstream gridf(gridfname);
	Grid grid(gridf, _table);

	GridCP gridcp;
	cout<<"Begin Sign"<<endl;
	gridcp.buildBasins(grid,0);
	gridcp.computeIntegrals(grid);
	gridcp.printCriticalPoints();
	vector<double> charges=gridcp.computeAIMCharges(grid);
	return charges;
}
/******************************************************************************************/
void Job::run() 
{
	string runType;
	string gridFileName;
	if(!readOneString("RunType",runType)) runType = "computeAIMCharges";
	if(!readOneString("Grid",gridFileName)) runType = "grid.cube";
	cout<<"----------------------------------------------------------"<<endl;
	cout<<"runType="<<runType<<endl;
	cout<<"----------------------------------------------------------"<<endl;
	printf("\n");

	if(to_upper(runType)=="HELP")
		printListOfRunTypes();
	else
	{
		if ( to_upper(runType) == to_upper("computeAIMCharges"))
		{
			vector<double> charges = computeAIMCharges(gridFileName);
		}
	}
}
