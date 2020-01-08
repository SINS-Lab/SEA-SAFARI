#pragma once
/**
 * This file contains some general function stubs for specifc
 * signaling flags, such as printing rest of debug file on exit.
 */ 

/**
 * This will call exit(EXIT_FAILURE),
 * after finishing the debug file, and printing
 * the reason to it.
 */ 
void exit_fail(std::string reason);