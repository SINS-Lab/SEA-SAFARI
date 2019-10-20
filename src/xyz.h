#ifndef XYZ_H_INCLUDED
#define XYZ_H_INCLUDED

#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>

struct XYZ_Single
{
    int number;
    //Number of values per row, note everything in the row
    //must be a double as currently implemented.
    int num_per_row;
    std::string comment;
    double **values;
    std::string *atoms;

    void load(std::ifstream& input, int number);
    void save(std::ofstream& output);
};

struct XYZ
{
    std::vector<XYZ_Single> xyzs;

    void load(std::ifstream& input);
    void save(std::ofstream& output);
};

#endif // XYZ_H_INCLUDED