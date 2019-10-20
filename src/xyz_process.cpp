#include "xyz.h"
#include <cstdio>

int main(int argc,char* argv[])
{
    if(argc==1) 
    {
        std::cout << "Please include argument of file to parse" << std::endl;
        return 1;
    }
    std::string file_name = argv[1];
    XYZ xyz;
    std::ifstream input;
    input.open(file_name);
    xyz.load(input);

    return 0;
}