using namespace std;
#include <JobControl/Job.h>
#include <Cube/GridCP.h>
#include <Cube/functions.h>
#include <Cube/Grid.h>
#include <Common/Element.h>
#include <Common/Descriptors.h>
#include <Becke/Becke.h>
#include <Orbitals/Orbitals.h>
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
		t = std::regex_replace(t, std::regex(" +$"), ""); // removing only trailing spaces
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
	string value;
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
	string value;
	if(readOneString(tag, value))
	{
		for(char& ch:value)
		{
			if(ch==',' or ch==';')
			{
				ch=' ';
			}
		}
		T a;
		stringstream ss(value);
		while(!ss.eof())
		{
			ss>>a;
			if(ss.fail() and !ss.eof()) 
			{	
				return false;
			}
			x.push_back(a);
			if (ss.peek() == EOF) 
			{
				break;
			}
		}
		return true;
	}
	return false;
}
/******************************************************************************************/
void Job::setJobList()
{
	_jobsList = {"Help Me Pls!!!", "computePartialCharges", "computeDescriptorsFromCubes", "computeIntegrals", "computeGridDifference", "MakeDensityCube", "MakeOrbitalsCube"}; 
	_jobDescription = {"Details are given for the available jobs run by this program.\nExample input fles for each job are also given. In this format, comment lines are specified by # at the start of the line",
		"Grid based computations of partial charges of the molecule. We provide 5 ways of computing atomic volumes, the first 3 of which are based on Bader's Atoms in molecule.\n\n **on-grid** : follows Tang's algorithm to find Bader volumes.\n **near-grid** : more precise version of on-grid.\n **near-grid-refinement** : even more precise. Requires more time.\n **VDD** topological method : assigns points to volumes by distance to closest atom.\n **Becke** : uses a regular density grid to interpolate Becke's atomic variable grids.\n\n Example format for input file :\n\n#RunType\n#RunType=Help\nRunType=ComputePartialCharges\n#GridFileName\nGrids=h2o_80_0.gcube \nPartitionMethod=on-grid\n\nW. Tang, E. Sanville, G. Henkelman, A grid-based bader analysis algorithm without lattice bias, Journal of Physics: Condensed Matter 21 (8) (2009) 084204.",
	       	"Computation of chemical descriptors from cube files using the same techniques as computePartialCharges job. Requires cube files of nucleophilic, electrophilic and radical attacks for the molecule. Energies must also be given by the user:\nif two are given, they are assumed to be the ionisation potential and the electronic affinity. If 3 are given they are assumed to be the total energies of each cube file.\n\n Example format for input file :\n\n#RunType=Help\n#RunType=ComputeDescriptorsFromCubes\n#GridFileName\nGrids=grid1.cube, grid2.cube, grid3.cube\nPartitionMethod=on-grid\nEnergies=I, A or E1,E2,E3",
		"Compute local integrals of grids on volumes defined by method of choice. A grid is required to define the volumes.\nThe additional grids provided by the user should contain the quantities to be integrated.\n\n **on-grid** : to define volumes using on-grid AIM. Requires electronic density grid.\n **near-grid** : to define volumes using near-grid AIM. Requires electronic density grid.\n **near-grid-refinement : to define volumes using near-grid-refinement AIM. Requires electronic density grid.\n **VDD** : to define volumes by distance to atoms. Can use any type of density.\n **BBS** : Build Basins By SIGN. Requires a grid of density difference. A job is provided in the program to obtain such a grid. An additional input *Cutoff=* is required for BBS that sets a threshold for insignificant values.\n **B2S** : Build 2 basins by SIGN. Same as BBS but only constructs two volumes.\n\n Example format for input file :\n\n#RunType=Help\n#RunType=ComputeIntegrals\n#GridFileName\nGrids=gridDefiningVolumes.cube, grid1ToBeIntegrated.cube, grid2ToBeIntegrated.cube\nPartitionMethod=BBS\nCutoff=1e-10", "Computes the differences of values of the first two grids provides and assigns them to the third.\n\n Example format for input file : \n\n#Runtype=Help\nRunType=ComputeDifference\n#GridFileName\nGrids=in1.cube, in2.cube, out.cube ",
	"Create a density grid and save it in .cube format. .wfx , .fchk , .molden , .gab and .log are supported as input files.\nthe user can choose from 3 standard grid sizes:\ncoarse ( 3 pts / Bohr)\nMedium (6 pts / Bohr)\nFine (12 pts / Bohr)\n\nA custom size is also provided in which the user enters the domain data as follows:\nNx, Ny, Nz, Ox, Oy, Oz, T11, T12, T13, T21, T22, T23, T31, T32, T33\nWhere N is the number of points in the ith direction, Oi are the coordinates of the bottom left corner of the cube and Tij are the coeficients of the translation vector.\n\n Example format for input file : \n\n#RunType=Help\nRunType=MakeDensityCube\n#GridFileName\nAnalyticFile=filename.wfx\nSize=Custom\nCustomSizeData=80,80,80,5,5,5,0.15,0,0,0,0.15,0,0,0,0.15\nGrid=save.cube ",
	"Compute a grid of molecular orbitals' values and save it in .cube format. All parameters for the grid domain are the same as MakeDensityCube. Additional input lines are required for the computation of molecular orbitals.\nThe user must specify which orbitals took take into account:\n All : **All**\n Occupied : **Occ**\n Virtual : **Virtual**\n Homo : **Homo**\n Lumo : **Lumo**\n Homo and lumo : **Homo, Lumo**\n Custom : **OrbitalsList=Orbital number specified by user**\nBy default the program will run with all MOs.\n\nThe choice of spin is also given:\n **SpinType=Alpha**\n **SpinType=Beta**\n **SpinType=Alpha, Beta**\n\nIf the user provides a custom list of orbitals the user can provide a list of spins corresponding to each orbital. This is done in **SpinList=alpha, beta, ...**.\nIf SpinList is shorter n length than OrbitalsList the program will fill the rest of the list with the last value read in the list"};
}
void Job::printListOfRunTypes()
{
	cout<<"Available runType :"<<endl;
	for(size_t i=0;i<_jobsList.size();i++)
	{
		cout<<"----------------------------------------------------------------------------------"<<endl;
		cout<<_jobsList[i]<<endl;
		cout<<"----------------------------------------------------------------------------------"<<endl;
		cout<<_jobDescription[i]<<endl<<endl;
	}
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
void Job::buildBasins(GridCP& gcp, const string& GridFileName, int method)
{
	ifstream gridf(GridFileName);
	Grid g(gridf, _table);	
	gcp.buildBasins(g, method);

}
void Job::buildBasinsBySign(GridCP& gcp, const string& GridFileName,  double cutoff, bool two) 
{
	ifstream gridf(GridFileName);
	Grid g(gridf, _table);
	gridf.close();
	if(two)
	{
			gcp.build2BasinSign(g);
	}	
	else
	{
			gcp.buildBasinsBySign(g, cutoff);
	}	
}
void Job::computeLocalIntegrals(GridCP& gcp, const vector<string>& GridFileNames)
{
	for(size_t i=0;i<GridFileNames.size();i++)
	{
		ifstream f(GridFileNames[i]);
		Grid g(f,_table);
		f.close();
		gcp.computeIntegrals(g);
		gcp.printCriticalPoints();
	}
}
void Job::ComputeGridDifference(const string& GFN1, const string& GFN2 , const string& NameNew)
{
	ifstream in1(GFN1);
	Grid g(in1, _table);
	ifstream in2(GFN2);
	Grid h(in2, _table);
	Grid diff=g-h;
	in1.close();
	in2.close();
	ofstream out(NameNew, ios::out);
	if(out.fail()) cout<<"failed to write to file"<<NameNew<<"..."<<endl;
	diff.save(out);
	cout<<"Grid has been saved to : " <<NameNew<<endl;
}
void Job::printCriticalPoints()
{
}
vector<double> Job::computePartialCharges(const string& gridfname, int method)
{
	vector<double> charges(3);
	if(method==4)
	{
		Factorial fact(100);
		Binomial bino (100, fact);
		ifstream gridf(gridfname);
		Grid grid(gridf, _table);
		Becke B(grid);
		B.partial_charge(grid);
		charges=B.get_Partial_Charge();
		B.printCharges();
	}
	else
	{
		ifstream gridf(gridfname);
		Grid grid(gridf, _table);
		GridCP gridcp;
		gridcp.buildBasins(grid,method);
		gridcp.computeIntegrals(grid);
		gridcp.printCriticalPoints();
		charges=gridcp.computeAIMCharges(grid);
	}
	return charges;
}
Descriptors Job::computeDescriptors(const string& GridFileName1, const string& GridFileName2, const string& GridFileName3, double I, double A, int AIMmethod)
{
	ifstream grid1(GridFileName1);
	ifstream grid2(GridFileName2);
	ifstream grid3(GridFileName3);
	Descriptors D(grid1, grid2, grid3, I, A, AIMmethod);
	return D;
}	
Descriptors Job::computeDescriptors(const string& GridFileName1, const string& GridFileName2, const string& GridFileName3, vector<double>E,  int AIMmethod)
{
	ifstream grid1(GridFileName1);
	ifstream grid2(GridFileName2);
	ifstream grid3(GridFileName3);
	Descriptors D(grid1, grid2, grid3, E, AIMmethod);
	return D;
}
template<typename T> Orbitals Job::computeOrbitals(const string& analyticFileName)
{
	Factorial fact(100);
	Binomial bino(100, fact);
	ifstream analyticFile(analyticFileName);
	T anaClass(analyticFile);
	analyticFile.close();
	Orbitals orb(anaClass, bino, _table);
	return orb;
}
void Job::createCube(Orbitals& orb, const Domain& d, const string& cubeFileName, bool orbFlag, vector<int> nums={0}, vector<int> typesSpin={0})
{
	Grid g;
	if(orbFlag)
	{
		g=orb.makeOrbGrid(d,nums, typesSpin);

	}
	else
	{
		g=orb.makeGrid(d);
	}
	ofstream out(cubeFileName);
	g.save(out);
	out.close();
}
Domain Job::DomainForCube(Orbitals& orb, const string& size,const vector<double>& csizes,const int& Nval)
{
	Domain d;
	if(to_upper(size)==to_upper("Coarse"))
	{
		double size=d.sizeUpMol(orb.get_struct(), 1.5);
		double xmax = size;
		double ymax = size;
		double zmax = size;
		int N1 = floor(size*6);
		int N2 = floor(size*6);
		int N3 = floor(size*6);
		vector<double> tx = {2*xmax/(N1-1),0,0};
		vector<double> ty = {0,2*ymax/(N2-1),0};
		vector<double> tz = {0,0,2*zmax/(N3-1)};
		vector<vector<double>>  t = {tx, ty, tz};
		d.set_all(Nval, N1, N2, N3, xmax, ymax, zmax, t);
	}
	else if(to_upper(size)==to_upper("Medium"))
	{
		
		double size=d.sizeUpMol(orb.get_struct(), 1.5);
		double xmax = size;
		double ymax = size;
		double zmax = size;
		int N1 = floor(size*12);
		int N2 = floor(size*12);
		int N3 = floor(size*12);
		vector<double> tx = {2*xmax/(N1-1),0,0};
		vector<double> ty = {0,2*ymax/(N2-1),0};
		vector<double> tz = {0,0,2*zmax/(N3-1)};
		vector<vector<double>>  t = {tx, ty, tz};
		d.set_all(Nval, N1, N2, N3, xmax, ymax, zmax, t);
	}
	else if(to_upper(size)==to_upper("Fine"))
	{
		double size=d.sizeUpMol(orb.get_struct(), 1.5);
		double xmax = size;
		double ymax = size;
		double zmax = size;
		int N1 = floor(size*24);
		int N2 = floor(size*24);
		int N3 = floor(size*24);
		vector<double> tx = {2*xmax/(N1-1),0,0};
		vector<double> ty = {0,2*ymax/(N2-1),0};
		vector<double> tz = {0,0,2*zmax/(N3-1)};
		vector<vector<double>>  t = {tx, ty, tz};
		d.set_all(Nval, N1, N2, N3, xmax, ymax, zmax, t);

	}
	else
	{
		vector<vector<double>> t(0);	
		vector<double> dx={csizes[6],csizes[7],csizes[8]};
		vector<double> dy={csizes[9],csizes[10],csizes[11]};
		vector<double> dz={csizes[12],csizes[13],csizes[14]};
		t.push_back(dx);
		t.push_back(dy);
		t.push_back(dz);

		d.set_all(Nval, int(csizes[0]), int(csizes[1]), int(csizes[2]), csizes[3], csizes[4], csizes[5], t);
	}
	return d;
}	
/******************************************************************************************/
void Job::setAllOrb(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, bool A, bool B, bool AnB, const int& N)
{
	orbnums.resize(N);
	for(int i=0; i<N;i++)
	{
		orbnums[i]=i;
	}
	if(!o.AlphaAndBeta() and A and !B)
	{
		orbspin.resize(N, 0);
	}
	else if(!o.AlphaAndBeta() and AnB)
	{
		orbspin.resize(N, 0);
		orbnums.resize(2*N,1);
		for(int i=N; i<2*N;i++)
		{
			orbnums[i]=i-N;
		}
	}
	else
	{
		orbspin.resize(N,1);
	}
}
void Job::setOccOrb(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, bool A, bool B, bool AnB,const int& N)
{
	vector<vector<double>> occ=o.OccupationNumber();
	if(!o.AlphaAndBeta() and A and !B) 
	{
		int k=0;
		for(int i=0; i<N;i++)
		{
			if(occ[0][i] > 1e-10)
			{
				orbnums[k]=i;
				k++;
			}
		}
		orbnums.resize(k);
		orbspin.resize(k,0);
	}
	else if(!o.AlphaAndBeta() and AnB) 
	{
		int k=0;
		orbnums.resize(2*N,0);
		for(int i=0; i<N;i++)
		{
			if(occ[0][i] > 1e-10)
			{
				orbnums[k]=i;
				k++;
			}
		}
		orbspin.resize(k,0);
		int j=k;
		for(int i =0;i<N;i++)
		{
			if(occ[1][i] > 1e-10)
			{
				orbnums[j]=i;
				j++;
			}
		}
		orbnums.resize(j);
		orbspin.resize(j,1);	
	}
	else
	{
		int j=0;
		for(int i=0;i<N;i++)
		{
			if(occ[1][i] > 1e-10)
			{
				orbnums[j]=i;
				j++;
			}
		}
		orbnums.resize(j);
		orbspin.resize(j,1);
	}

}
void Job::setVirtOrb(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, bool A, bool B, bool AnB,const int& N)
{
	vector<vector<double>> occ=o.OccupationNumber();
	if(!o.AlphaAndBeta() and A and !B)
	{
		int k=0;
		for(int i=0; i<N;i++)
		{
			if(occ[0][i] < 1e-10)
			{
				orbnums[k]=i;
				k++;
			}
		}
		orbnums.resize(k);
		orbspin.resize(k,0);
	}
	else if(!o.AlphaAndBeta() and AnB) 
	{
		int k=0;
		orbnums.resize(2*N,0);
		for(int i=0; i<N;i++)
		{
			if(occ[0][i] < 1e-10)
			{
				orbnums[k]=i;
				k++;
			}
		}
		orbspin.resize(k,0);
		int j=k;
		for(int i =0;i<N;i++)
		{
			if(occ[1][i] < 1e-10)
			{
				orbnums[j]=i;
				j++;
			}
		}
		orbnums.resize(j);
		orbspin.resize(j,1);	
	}
	else
	{
		int j=0;
		for(int i=0;i<N;i++)
		{
			if(occ[1][i] < 1e-10)
			{
				orbnums[j]=i;
				j++;
			}
		}
		orbnums.resize(j);
		orbspin.resize(j,1);
	}
}
void Job::setHomo(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, bool A, bool B, bool AnB)
{
	int i=0;
	vector<vector<double>> occ=o.OccupationNumber();
	if(!o.AlphaAndBeta() and A and !B)
	{
		while(occ[0][i] > 1e-10)
		{
			i++;
		}
		orbnums={i-1};
		orbspin={0};
	}
	else if(!o.AlphaAndBeta() and AnB)
	{
		int j=0;
		while(occ[0][i] > 1e-10)
		{
			i++;
		}
		while(occ[1][j] >1e-10)
		{
			j++;
		}
		orbnums={i-1,j-1};
		orbspin={0, 1};
	}
	else
	{
		int j=0;
		while(occ[1][j] >1e-10)
		{
			j++;
		}
		orbnums={j-1};
		orbspin={1};
	}
}
void Job::setLumo(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, bool A, bool B, bool AnB)
{
	int i=0;
	vector<vector<double>> occ=o.OccupationNumber();
	if(!o.AlphaAndBeta() and A and !B)
	{
		while(occ[0][i] > 1e-10)
		{
			i++;
		}
		orbnums={i};
		orbspin={0};
	}
	else if(!o.AlphaAndBeta() and AnB)
	{
		int j=0;
		while(occ[0][i] > 1e-10)
		{
			i++;
		}
		while(occ[1][j] >1e-10)
		{
			j++;
		}
		orbnums={i,j};
		orbspin={0, 1};
	}
	else
	{
		int j=0;
		while(occ[1][j] >1e-10)
		{
			j++;
		}
		orbnums={j};
		orbspin={1};
	}
}
void Job::setHomoLumo(vector<int>& orbnums, vector<int>& orbspin, Orbitals& o, bool A, bool B, bool AnB) 
{
	int i=0;
	vector<vector<double>> occ=o.OccupationNumber();
	if(!o.AlphaAndBeta() and B and !A)
	{
		while(occ[0][i] > 1e-10)
		{
			i++;
		}		
		orbnums={i-1, i};
		orbspin={0,0};
	}
	else if(!o.AlphaAndBeta() and AnB)
	{
		int j=0;
		while(occ[0][i] > 1e-10)
		{
			i++;
		}
		while(occ[1][j] >1e-10)
		{
			j++;
		}
		orbnums={i-1,j-1, i,j};
		orbspin={0, 1, 0, 1};
	}
	else
	{
		int j=0;
		while(occ[1][j] >1e-10)
		{
			j++;
		}
		orbnums={j-1, j};
		orbspin={1, 1};
	}
}
void Job::setCustom(vector<int>& orbnums, vector<int>& orbspin, vector<string>& tmplist) 
{
	for(size_t i=0; i<tmplist.size();i++)
	{
		if(to_upper(tmplist[i])==to_upper("alpha"))
		{
			orbspin.push_back(0);
		}
		else if(to_upper(tmplist[i])==to_upper("beta"))
		{
			orbspin.push_back(1);
		}
		else
		{
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cerr<<"No such spin. Please check input file. "<<endl;
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			exit(1);
		}
	}
	for(size_t i=0; i<orbnums.size(); i++)
	{
		orbnums[i] -= 1;
	}
	if(orbspin.size()<orbnums.size())
	{	
		int last=orbspin.back();
		for(size_t i=orbspin.size(); i<orbnums.size(); i++)
		{
			orbspin.push_back(last);
		}
	}
}
/******************************************************************************************/
void Job::run() 
{
	string runType;
	if(!readOneString("RunType",runType)) runType = "HELP";
	cout<<"----------------------------------------------------------"<<endl;
	cout<<"runType="<<runType<<endl;
	cout<<"----------------------------------------------------------"<<endl;
	//printf("\n");

	if(to_upper(runType)=="HELP")
	{
		printListOfRunTypes();
	}
	else if ( to_upper(runType) == to_upper("computePartialCharges"))
	{
		vector<string> gridFileName;
		if(!readListType<string>("Grids",gridFileName)) gridFileName.push_back("grid.cube");
		int Method;
		string Mdef;
		if(!readOneString("PartitionMethod",Mdef))
		{
			Method=0;
			cout<<"----------------------------------------------------------"<<endl;
			cout<<"volume partition method : "<<Mdef<<endl;
			cout<<"----------------------------------------------------------"<<endl;
		}
		else if(to_upper(Mdef) == to_upper("on-grid"))
		{
			Method=0;
			cout<<"----------------------------------------------------------"<<endl;
			cout<<"volume partition method : "<<Mdef<<endl;
			cout<<"----------------------------------------------------------"<<endl;
		}
		else if(to_upper(Mdef) == to_upper("near-grid"))
		{
			Method=1;
			cout<<"----------------------------------------------------------"<<endl;
			cout<<"volume partition method : "<<Mdef<<endl;
			cout<<"----------------------------------------------------------"<<endl;
		}
		else if(to_upper(Mdef) == to_upper("near-grid-refinement"))
		{
			Method=2;
			cout<<"----------------------------------------------------------"<<endl;
			cout<<"volume partition method : "<<Mdef<<endl;
			cout<<"----------------------------------------------------------"<<endl;
		}
		else if(to_upper(Mdef) == to_upper("VDD"))
		{
			Method=3;
			cout<<"----------------------------------------------------------"<<endl;
			cout<<"volume partition method : "<<Mdef<<endl;
			cout<<"----------------------------------------------------------"<<endl;
		}
		else if(to_upper(Mdef) == to_upper("Becke"))
		{
			Method=4;
			cout<<"----------------------------------------------------------"<<endl;
			cout<<"volume partition method : "<<Mdef<<endl;
			cout<<"----------------------------------------------------------"<<endl;
		}
		else
		{
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cerr<<"Sorry, Partition method type unknown. Please check PartitionMethod input. "<<endl;
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			exit(1);
		}
		vector<double> charges = computePartialCharges(gridFileName[0],Method);
	}
	else if ( to_upper(runType) == to_upper("ComputeDescriptorsFromCube"))
	{
		vector<string> gridFileNames;
		vector<double> E;
		int Method;
		string Mdef;
		if(!readListType<string>("Grids",gridFileNames)) 
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

		if(!readOneString("PartitionMethod",Mdef))
		{
			Method=0;
			cout<<"----------------------------------------------------------"<<endl;
			cout<<"volume partition method : "<<Mdef<<endl;
			cout<<"----------------------------------------------------------"<<endl;
		}
		else if(to_upper(Mdef) == to_upper("on-grid"))
		{
			Method=0;
			cout<<"----------------------------------------------------------"<<endl;
			cout<<"volume partition method : "<<Mdef<<endl;
			cout<<"----------------------------------------------------------"<<endl;
		}
		else if(to_upper(Mdef) == to_upper("near-grid"))
		{
			Method=1;
			cout<<"----------------------------------------------------------"<<endl;
			cout<<"volume partition method : "<<Mdef<<endl;
			cout<<"----------------------------------------------------------"<<endl;
		}
		else if(to_upper(Mdef) == to_upper("near-grid-refinement"))
		{
			Method=2;
			cout<<"----------------------------------------------------------"<<endl;
			cout<<"volume partition method : "<<Mdef<<endl;
			cout<<"----------------------------------------------------------"<<endl;
		}
		else if(to_upper(Mdef) == to_upper("VDD"))
		{
			Method=3;
			cout<<"----------------------------------------------------------"<<endl;
			cout<<"volume partition method : "<<Mdef<<endl;
			cout<<"----------------------------------------------------------"<<endl;
		}
		else if(to_upper(Mdef) == to_upper("Becke"))
		{
			Method=4;
			cout<<"----------------------------------------------------------"<<endl;
			cout<<"volume partition method : "<<Mdef<<endl;
			cout<<"----------------------------------------------------------"<<endl;
		}
		else
		{
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cerr<<"Sorry, volume partitioning method type unknown. Please check AIMmethod input. "<<endl;
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			exit(1);
		}

		if(!readListType<double>("Energies", E))
		{
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cerr<<"Sorry:faulty Energies input."<<endl;
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			exit(1);
		}
		else if(E.size()==2)
		{
			cout<< "Reading Ionisation potential I= "<< E[0]<<" and Electronic affinity A = "<< E[1]<<endl;
			Descriptors D=computeDescriptors(gridFileNames[0], gridFileNames[1], gridFileNames[2], E[0], E[1],Method);
			cout<<D;
		}
		else if(E.size()==3)
		{
			cout<<" Reading Energies of files"<<endl;
			Descriptors D=computeDescriptors(gridFileNames[0], gridFileNames[1], gridFileNames[2], E, Method);
			cout<<D;
		}
		else
		{
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cerr<<"Sorry:faulty Energies input :not enough/too many values." << endl;
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			exit(1);
		}
	}
	else if(to_upper(runType) == to_upper("ComputeIntegrals"))
	{
		int Method;
		vector<string> GridFileNames;
		string Mbuild;
		double C;
		GridCP gcp;
		if(!readOneString("PartitionMethod",Mbuild))
		{
			Method=0;
			cout<<"----------------------------------------------------------"<<endl;
			cout<<"volume build method : "<<Mbuild<<endl;
			cout<<"----------------------------------------------------------"<<endl;
		}
		else if(to_upper(Mbuild) == to_upper("on-grid"))
		{
			Method=0;
			cout<<"----------------------------------------------------------"<<endl;
			cout<<"volume build method : "<<Mbuild<<endl;
			cout<<"----------------------------------------------------------"<<endl;
		}
		else if(to_upper(Mbuild) == to_upper("near-grid"))
		{
			Method=1;
			cout<<"----------------------------------------------------------"<<endl;
			cout<<"volume build method : "<<Mbuild<<endl;
			cout<<"----------------------------------------------------------"<<endl;
		}
		else if(to_upper(Mbuild) == to_upper("near-grid-refinement"))
		{
			Method=2;
			cout<<"----------------------------------------------------------"<<endl;
			cout<<"volume build method : "<<Mbuild<<endl;
			cout<<"----------------------------------------------------------"<<endl;
		}
		else if(to_upper(Mbuild) == to_upper("VDD"))
		{
			Method=3;
			cout<<"----------------------------------------------------------"<<endl;
			cout<<"volume build method : "<<Mbuild<<endl;
			cout<<"----------------------------------------------------------"<<endl;
		}
		else if(to_upper(Mbuild) == to_upper("BBS"))
		{
			Method=4;
			cout<<"----------------------------------------------------------"<<endl;
			cout<<"volume build method : "<<Mbuild<<endl;
			cout<<"----------------------------------------------------------"<<endl;
		}
		else if(to_upper(Mbuild) == to_upper("B2S"))
		{
			Method=5;
			cout<<"volume build method : "<<Mbuild<<endl;
		}
		else
		{
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cerr<<"Sorry, volume build method unknown. Please check input file. "<<endl;
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			exit(1);
		}

		if(!readListType<string>("Grids", GridFileNames))
		{
			for(size_t i=0;i<GridFileNames.size(); i++)
			{
				GridFileNames[i]="grid" +to_string(i)+".cube";
			}
		}
		
		if(!readOneType<double>("Cutoff", C))
		{
			C=0;
		}
		if(Method==4)
		{
			buildBasinsBySign(gcp, GridFileNames[0], C, false);
		}
		else if(Method==5)
		{
			buildBasinsBySign(gcp, GridFileNames[0], C, true);
		}
		else
		{
			buildBasins(gcp, GridFileNames[0], Method);
		}
		cout<<GridFileNames[0]<<" taken as input for Basin build... Proceeding... "<<endl;
		computeLocalIntegrals(gcp, GridFileNames); 
	}
	else if(to_upper(runType) == to_upper("ComputeGridDifference"))
	{
		vector<string> GridFileNames;
		cout<<"yes"<<endl;
		if(!readListType<string>("Grids", GridFileNames))
		{
			cout<<"yes"<<endl;
			GridFileNames[0]="grid1.cube";
			GridFileNames[1]="grid2.cube";
			GridFileNames[2]="grid3.cube";
		}
		ComputeGridDifference(GridFileNames[0], GridFileNames[1], GridFileNames[2]);		
	}
	else if(to_upper(runType) == to_upper("MakeDensityCube"))
	{	
		vector<double> csizes(0);
		string size;
		string ANAFileName;
		vector<string> gridFileName;
		Orbitals o;
		if(!readOneString("AnalyticFile", ANAFileName))
		{
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cerr<<"Sorry, no analytic file found. Please check input file. "<<endl;
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			exit(1);
		}
		else if(ANAFileName.find(".wfx")!=string::npos)
		{
			cout<<"Reading data from"<<ANAFileName<<" ...Please wait"<<endl;
			o=computeOrbitals<WFX>(ANAFileName);
		}
		else if(ANAFileName.find(".fchk")!=string::npos)
		{
			cout<<"Reading data from"<<ANAFileName<<" ...Please wait"<<endl;
			o=computeOrbitals<FCHK>(ANAFileName);
		}
		else if(ANAFileName.find(".molden")!=string::npos)
		{
			cout<<"Reading data from"<<ANAFileName<<" ...Please wait"<<endl;
			o=computeOrbitals<MOLDENGAB>(ANAFileName);
		}
		else if(ANAFileName.find(".gab")!=string::npos)
		{
			cout<<"Reading data from"<<ANAFileName<<" ...Please wait"<<endl;
			o=computeOrbitals<MOLDENGAB>(ANAFileName);
		}
		else if(ANAFileName.find(".log")!=string::npos)
		{
			cout<<"Reading data from"<<ANAFileName<<" ...Please wait"<<endl;
			o=computeOrbitals<LOG>(ANAFileName);
		}
		else
		{
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cerr<<"Sorry, unknown file format for analytic file. Please check input file. "<<endl;
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			exit(1);
		}
		if(!readOneString("Size", size))
		{
			cout<<" No grid size chosen... Making medium grid (6pts / Bohr) "<<endl;
			size="Medium";
		}
		if(to_upper(size) == to_upper("Custom") and !readListType<double>("CustomSizeData", csizes))
		{
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cerr<<"Sorry, Custom sized grid requested but no CustomSizeData specified. Please check input file. "<<endl;
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			exit(1);
		}
		if(to_upper(size) == to_upper("Custom") and csizes.size()!=15)
		{
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cerr<<"Sorry, missing values for grid customsizes. Please check input file. "<<endl;
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			exit(1);
		}
		if(!readListType<string>("Grid", gridFileName))
		{
			cout<<" No grid name chosen... Saving as grid.cube"<<endl;
			gridFileName[0]="grid.cube";
		}
		cout<<"Writing... Please wait"<<endl;
		Domain d=DomainForCube(o, size, csizes, 1);
		createCube(o,d,gridFileName[0],false);
		cout<<"Data saved to file : "<<gridFileName[0]<<endl;
	}
	else if(to_upper(runType) == to_upper("MakeOrbitalsCube"))
	{	
		vector<double> csizes(0);
		string size;
		string ANAFileName;
		vector<string> OrbType(0);
		vector<string> gridFileName;
		Orbitals o;
		if(!readOneString("AnalyticFile", ANAFileName))
		{
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cerr<<"Sorry, no analytic file found. Please check input file. "<<endl;
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			exit(1);
		}
		else if(ANAFileName.find(".wfx")!=string::npos)
		{
			cout<<"Reading data from"<<ANAFileName<<" ...Please wait"<<endl;
			o=computeOrbitals<WFX>(ANAFileName);
		}
		else if(ANAFileName.find(".fchk")!=string::npos)
		{
			cout<<"Reading data from"<<ANAFileName<<" ...Please wait"<<endl;
			o=computeOrbitals<FCHK>(ANAFileName);
		}
		else if(ANAFileName.find(".molden")!=string::npos)
		{
			cout<<"Reading data from"<<ANAFileName<<" ...Please wait"<<endl;
			o=computeOrbitals<MOLDENGAB>(ANAFileName);
		}
		else if(ANAFileName.find(".gab")!=string::npos)
		{
			cout<<"Reading data from"<<ANAFileName<<" ...Please wait"<<endl;
			o=computeOrbitals<MOLDENGAB>(ANAFileName);
		}
		else if(ANAFileName.find(".log")!=string::npos)
		{
			cout<<"Reading data from"<<ANAFileName<<" ...Please wait"<<endl;
			o=computeOrbitals<LOG>(ANAFileName);
		}
		else
		{
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cerr<<"Sorry, unknown file format for analytic file. Please check input file. "<<endl;
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			exit(1);
		}
		if(!readOneString("Size", size))
		{
			cout<<" No grid size chosen... Making medium grid (6pts / Bohr) "<<endl;
			size="Medium";
		}
		if(to_upper(size) == to_upper("Custom") and !readListType<double>("CustomSizeData", csizes))
		{
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cerr<<"Sorry, Custom sized grid requested but no CustomSizeData specified. Please check input file. "<<endl;
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			exit(1);
		}
		if(to_upper(size) == to_upper("Custom") and csizes.size()!=15)
		{
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cerr<<"Sorry, missing values for grid customsizes. Please check input file. "<<endl;
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			exit(1);
		}
		if(!readListType<string>("Grid", gridFileName))
		{
			cout<<" No grid name chosen... Saving as grid.cube"<<endl;
			gridFileName[0]="grid.cube";
		}

		int N=o.NumberOfMo();
		vector<int> orbspin(0);
		vector<int> orbnums(N);
		vector<string> SpinType(0);

		if(!readListType<string>("SpinType", SpinType))
		{
			cout<<"No spin type selected... Computing both alpha and beta." <<endl;
		}
		if(SpinType.size()<2) SpinType.resize(2,"");
		bool A=(to_upper(SpinType[0])==to_upper("alpha") or to_upper(SpinType[1])==to_upper("alpha"));
		bool B=(to_upper(SpinType[0])==to_upper("alpha") or to_upper(SpinType[1])==to_upper("alpha"));
		bool AnB=(A and B);

		if(!readListType<string>("OrbType", OrbType))
		{
			cout<<"No orbitals selected for computation... Taking ALL"<<endl;
			setAllOrb(orbnums, orbspin, o, A, B, AnB, N);
		}
		else if(to_upper(OrbType[0])==to_upper("All"))
		{
			orbnums.resize(N);
			setAllOrb(orbnums, orbspin, o, A, B, AnB, N);
		}
		else if(to_upper(OrbType[0])==to_upper("Occupied"))
		{
			setOccOrb(orbnums, orbspin, o, A, B, AnB, N);
		}
		else if(to_upper(OrbType[0])==to_upper("Virtual"))
		{
			setVirtOrb(orbnums, orbspin, o, A, B, AnB, N);
		}
		else if(OrbType.size()==1 and to_upper(OrbType[0])==to_upper("Homo"))
		{
			setHomo(orbnums, orbspin, o, A, B, AnB);
		}
		else if(OrbType.size()==1 and to_upper(OrbType[0])==to_upper("lumo"))
		{
			setLumo(orbnums, orbspin, o, A, B, AnB);
		}
		else if(OrbType.size()>1 and ((to_upper(OrbType[0])==to_upper("lumo") and to_upper(OrbType[1])==to_upper("Homo")) or (to_upper(OrbType[1])==to_upper("lumo") and to_upper(OrbType[0])==to_upper("Homo"))))
		{
			setHomoLumo(orbnums, orbspin, o, A, B, AnB);
		}
		else if(to_upper(OrbType[0]) == to_upper("Custom"))
		{
			vector<string> tmplist(0);
			orbnums.resize(0);
			if(!readListType<int>("OrbitalsList", orbnums))
			{
				cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
				cerr<<"No orbitals specified for computation. Please check input file. "<<endl;
				cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
				exit(1);
			}
			if(!readListType<string>("SpinList", tmplist))
			{
				cout<<"WARNING!!! No spin list provided for the orbitals occupation... Alpha spin taken."<<endl;
				orbspin.resize(orbnums.size(), 0);
			}
			else
			{
				setCustom(orbnums, orbspin, tmplist);
			}
		}
		else
		{
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cerr<<"Unknown orbitals type. Please check input file. "<<endl;
			cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			exit(1);
		}
		cout<<"Writing... Please wait "<<endl;
		Domain d=DomainForCube(o, size, csizes, orbnums.size());
		createCube(o,d,gridFileName[0],true, orbnums, orbspin);
		cout<<"Data saved to file : "<<gridFileName[0]<<endl;
	}
}
