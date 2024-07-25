#include<iostream>
#include <cmath>
#include <sys/time.h>
#include "Timer.h"
#include <JobControl/Job.h>
#include <string>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{
	string fname="inputhelp.txt";
	if(argc<2)
	{
		cout<<"You have to give the name of the input file"<<endl;
		cout<<"Example : "<<argv[0]<<" "<<"input.cdft"<<endl;
		cout<<"Here a list of all available methods"<<endl;
		ofstream f(fname);
		f<<"RunType=HELP"<<endl;
		f.close();
	}
	else 
		fname=argv[1];
	Timer timer;
	Job job(fname);
	job.run();
	cout<<"Time in ms "<<timer.get()<<endl;
	return 0;
}
