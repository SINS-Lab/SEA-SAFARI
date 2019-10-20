#ifndef STRING_UTILS_H_INCLUDED
#define STRING_UTILS_H_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>

std::vector<std::string> split(std::string input);
double* to_double_array(std::vector<std::string> input, int start, int end);

#endif // STRING_UTILS_H_INCLUDED