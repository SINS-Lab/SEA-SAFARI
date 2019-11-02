#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <map>

class ArgValue
{
private:
    std::string raw_value;
public:
    ArgValue(){raw_value="";}
    ArgValue(const std::string& value){raw_value = value;}
    bool as_bool();
    double as_double();
    std::string as_string();
    int as_int();
    /**
     * This returns whether there is a raw value,
     * use as_bool() for checking if the value is
     * valid for true.
     */ 
    operator bool() { return raw_value != "";}
};

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

/**
 * Returns whether string starts with start
 */ 
bool starts_with(std::string string, const char* start);

/**
 * This returns a map of arguments given, It assumes the following:
 * 
 * Argument keys start with -, and are a single entry.
 * Argument values directly follow the key, and do not start with -
 * if no value is given, the value is taken to be "t", ie a boolean flag
 */ 
std::map<std::string, ArgValue> get_arguments(int argc, char* argv[]);