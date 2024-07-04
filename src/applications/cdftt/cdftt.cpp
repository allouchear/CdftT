#include<iostream>
#include <cmath>
#include <numeric/functions.h>
#include <numeric/Grid.h>
#include <numeric/GridCP.h>
#include <common/Element.h>
#include <sys/time.h>
#include "Timer.h"
#include <common/Descriptors.h>
#include <jobcontrol/Job.h>

using namespace std;

int main()
{
	Timer timer;
	Job job;
	job.run();
	cout<<"Time in ms "<<timer.get()<<endl;
	return 0;
}
