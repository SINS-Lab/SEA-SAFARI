#pragma once
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>

struct XYZ_Single
{
    //Number of particles in this XYZ
    int number;
    //Number of values per row, note everything in the row
    //must be a double as currently implemented.
    int num_per_row;
    //Comment for the XYZ
    std::string comment;
    //Arrays of positions for the particles.
    //Note that the first 3 are x,y,z, the remainder
    //can be anything.
    double **values;
    //List of atomic symbols for the paricles.
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