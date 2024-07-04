#include<iostream>
#include <cmath>
#include <sys/time.h>
#include "Timer.h"
#include <JobControl/Job.h>

using namespace std;

int main()
{
	Timer timer;
	Job job;
	job.run();
	cout<<"Time in ms "<<timer.get()<<endl;
	return 0;
}
