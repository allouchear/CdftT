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
template<typename T> bool Job::readOneType(const string& tag, T& x)
{ 
	string& value;
	if(readOneString(tag, value))
	{
		size_t i=0;
		while(i<value.size())
		{
			if(value[i]==',' or value[i]==';')
			{
				value[i]=' ';
			}
			i++;
		}
		stringstream ss(value);
		ss>>x;
		if(ss.fail()) 
		{	
			return false;
		}	
		return true;
	}
	return false;
}
/******************************************************************************************/
template<typename T> bool Job::readListType(const string& tag, vector<T>& x)
{
	string& value;
	if(readOneString(tag, value))
	{
		size_t i=0;
		while(i<value.size())
		{
			if(value[i]==',' or value[i]==';')
			{
				value[i]=' ';
			}
			i++;
		}

		T a;
		stringstream ss(value);

		while(!ss.fail())
		{
			ss>>a;
			if(ss.fail()) 
			{	
				return false;
			}
			x.push_back(a);
		}
		return true;
	}
	return false;
}
/******************************************************************************************/
void Job::setJobList()
{
	_jobsList = {"computeAIMCharges" ,"Help","computeDescriptorsFromCubes"}; 
}
void Job::printListOfRunTypes()
{
	cout<<"Available runType :"<<endl;
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
vector<double> Job::computeAIMCharges(const string& gridfname, int method)
{
	ifstream gridf(gridfname);
	Grid grid(gridf, _table);

	GridCP gridcp;
	gridcp.buildBasins(grid,method);
	gridcp.computeIntegrals(grid);
	gridcp.printCriticalPoints();
	vector<double> charges=gridcp.computeAIMCharges(grid);
	return charges;
}
Descriptors computeDescriptors(const string& GridFileName1, const string& GridFileName2, const string& GridFileName3, double I, double A, int AIMmethod)
{
	ifstream grid1(GridFileName1);
	ifstream grid2(GridFileName2);
	ifstream grid3(GridFileName3);
	Descriptors D(grid1, grid2, grid3, I, A, AIMmethod);
	return D;
}	
/******************************************************************************************/
void Job::run() 
{
	string runType;
	string gridFileName;
	if(!readOneString("RunType",runType)) runType = "HELP";
	cout<<"----------------------------------------------------------"<<endl;
	cout<<"runType="<<runType<<endl;
	cout<<"----------------------------------------------------------"<<endl;
	printf("\n");

	if(to_upper(runType)=="HELP")
		printListOfRunTypes();
	else if ( to_upper(runType) == to_upper("computeAIMCharges"))
	{
		if(!readOneString("Grid",gridFileName)) gridFileName = "grid.cube";
		int Method;
		string Mdef;
		if(!readOneString("AIMmethod",Mdef))
		{
			Method=0;
		}
		else if(to_upper(Mdef) == to_upper("on-grid"))
		{
			Method=0;
		}
		else if(to_upper(Mdef) == to_upper("near-grid"))
		{
			Method=1;
		}
		else if(to_upper(Mdef) == to_upper("near-grid-refinement"))
		{
			Method=2;
		}
		else if(to_upper(Mdef) == to_upper("VDD"))
		{
			Method=3;
		}
		else
		{
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cerr<<"Sorry, AIM method type unknown. Please check AIMmethod input. "<<endl;
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			exit(1);
		}
		vector<double> charges = computeAIMCharges(gridFileName,Method);
	}
	else if ( to_upper(runType) == to_upper("ComputeDescriptorsFromCube"))
	{
		vector<string> gridFileNames;
		if(!readListType<string>("Grid",gridFileNames)) 
		{
			gridFileNames[0]="grid1.cube";
			gridFileNames[1]="grid2.cube";
			gridFileNames[2]="grid3.cube";
		}
		if(gridFileNames.size()>3)
		{
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cerr<<"Warning: too many grids. Proceeding with the first three. "<<endl;
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		}

	}
}
