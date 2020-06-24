#include "string_utils.h"
#include <string.h>

bool ArgValue::as_bool()
{
    return raw_value == "t";
}

double ArgValue::as_double()
{
    return atof(raw_value.c_str());
}

int ArgValue::as_int()
{
    return atoi(raw_value.c_str());
}

std::string ArgValue::as_string()
{
    return raw_value;
}

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

bool starts_with(std::string string, const char* start)
{
    return strncmp(string.c_str(), start, 1) == 0;
}

std::map<std::string, ArgValue> get_arguments(int argc, char* argv[])
{
    std::map<std::string, ArgValue> map;

    std::string key;
    std::string val;

    //Index 0 is executable name, we ignore it.
    for(int i = 1; i<argc; i++)
    {
        //We found a key!
        if(strncmp(argv[i], "-", 1)==0)
        {
            key = argv[i];
            //Assume value follows, or a "t" if end of args
            val = i+1 < argc ? argv[i+1] : "t";
            //Check that it isn't another key
            if(starts_with(val, "-")) val = "t";
            //Increment again if we actually had a value.
            if(val != "t")
            {
                i++;
            }
            map[key] = val;
        }
    }
    return map;
}

void findAndReplaceAll(std::string &data, std::string toSearch, std::string replaceStr)
{
    // Get the first occurrence
    size_t pos = data.find(toSearch);

    // Repeat till end is reached
    while (pos != std::string::npos)
    {
        // Replace this occurrence of Sub String
        data.replace(pos, toSearch.size(), replaceStr);
        // Get the next occurrence from the current position
        pos = data.find(toSearch, pos + replaceStr.size());
    }
}