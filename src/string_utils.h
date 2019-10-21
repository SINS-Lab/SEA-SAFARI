#ifndef STRING_UTILS_H_INCLUDED
#define STRING_UTILS_H_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>

/**
 * Splits the given string input into a vector of strings.
 * The string is split by whitespace characters.
 * 
 * @param input - string to split
 * @return vector of sub-strings
 */ 
std::vector<std::string> split(std::string input);

/**
 * Converts the given vector of strings into an array of doubles.
 * Note that the loop is from i=start to i=end, so if you want to
 * loop over the entire vector, send end as input.size()-1
 * 
 * @param input - the vector of strings
 * @param start - start index in the vector
 * @param end - the end index in the vector.
 * @return an array containing the strings converted to doubles.
 */ 
double* to_double_array(std::vector<std::string> input, int start, int end);

#endif // STRING_UTILS_H_INCLUDED