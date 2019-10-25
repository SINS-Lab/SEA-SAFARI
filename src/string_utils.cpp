#include "string_utils.h"

std::vector<std::string> split(std::string input)
{
    std::vector<std::string> args;

    std::istringstream iss(input);

    //Split the string into an array.
    do
    {
        std::string subs;
        iss >> subs;
        args.push_back(subs);
    }
    while (iss);

    return args;
}

double* to_double_array(std::vector<std::string> input, int start, int end)
{
    double *ret = new double[end - start + 1];
    for(int i = start; i<=end; i++)
    {
        //Index on the returned array should start at 0, rather than start.
        ret[i-start] = atof(input[i].c_str());
    }
    return ret;
}