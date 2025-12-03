#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/time.h>

#include "Timer.h"

#include <JobControl/Job.h>
#include <Utils/Utils.h>


int main(int argc, char* argv[])
{
    std::string fname = "inputhelp.txt";
    
    if(argc < 2)
    {
        std::stringstream errorMessage;
        errorMessage << "Error: input file name not provided. Please provide the input file name." << std::endl << std::endl;
        errorMessage << "Example:" << std::endl;
        errorMessage << "$ " << argv[0] << " " << "inputhelp.txt" << std::endl << std::endl;

        print_error(errorMessage.str());

        std::ofstream f(fname);
        f << "RunType=HELP" << std::endl;
        f.close();
    }
    else
    {
        fname = argv[1];
    }

    Timer timer;
    Job job(fname);
    job.run();

    std::cout << "Time in ms " << timer.get() << std::endl;

    return 0;
}
