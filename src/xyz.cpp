#include "xyz.h"
#include "string_uitls.h"

void XYZ_Single::load(std::ifstream& input, int number)
{
    //Store the number
    this->number = number;
    //First line is comment
    getline(input, comment);

    values = new double*[number];
    atoms = new std::string[number];

    std::string line;

    //For number lines, we read the line.
    for(int i = 0; i<number; i++)
    {
        getline(input, line);
        std::vector<std::string> args = split(line);
        atoms[i] = args[0];
        num_per_row = args.size()-1;
        values[i] = to_double_array(args, 1, num_per_row);
    }
}

void XYZ_Single::save(std::ofstream& output)
{
    output << number << std::endl;
    output << comment << std::endl;
    for(int i = 0; i<number; i++)
    {
        output << atoms[i] << " ";
        int num = num_per_row;
        for(int j = 0; j<num; j++)
        {
            output << values[i][j];
            if(j < num-1) output << " ";
        }
        output << std::endl;
    }
}

void XYZ::load(std::ifstream& input)
{
    std::string line;
    while(getline(input, line))
    {
        int num = atoi(line.c_str());
        XYZ_Single xyz;
        xyz.load(input, num);
        xyzs.push_back(xyz);
    }
}

void XYZ::save(std::ofstream& output)
{
    for(XYZ_Single xyz: xyzs)
    {
        xyz.save(output);
    }
}