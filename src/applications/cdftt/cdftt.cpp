#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <JobControl/JobManager.hpp>


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

        std::ofstream exampleFile(fname);
        if (exampleFile)
        {
            exampleFile << "RunType=HELP" << std::endl;
            exampleFile.close();
        }
    }
    else
    {
        fname = argv[1];
    }

    JobManager jobManager(fname);
    jobManager.runJobs();

    return 0;
}
