# Sea-Safari

This software is for simulating scattering of low and hyperthermal energy ions off solid surfaces. It is a continuation of SAFARI, which filled the same role, however was limited to specific crystal faces.

This version of Safari allows using an arbitrary input crystal, and is optimized for working on larger sized crystals than the old one, as a result, it runs much faster for larger crystals, however uses significantly more RAM.

Unlike SAFARI, this version only outputs to an ascii file, and also includes the entire lattice in the xyz output, instead of just the interacting atoms.

## Setup and Testing:

A Simple setup/test would be to download and run ./setup_and_test.sh. This will do the following:

- Clone this repository to a directory called "SAFARI" *Warning, it will also remove that directory if it exits!*
- Clone the SAFARI-ANALYSIS repository to "SAFARI-ANALYSIS" as well *Warning, it will also remove that directory if it exits!*
- build SAFARI, and copy executables to a "runs" directory and a "tests" directory, as well as to "SAFARI-ANALYSIS"
- Run the standard SAFARI Test runs in "./tests"

This requires the following:
- make
- g++

For using on Windows, use the Windows Subsystem for Linux, and make sure you have make and g++ installed on there first!

## Analysing results:

Analysing results is generally done using the safari_detect gui found in SAFARI-ANALYSIS, which the above setup will automatically download. This requires the following:

- Python3
- numpy
- matplotlib
- scipy
- mendeleev

Commandline analysis scripts can also be made if needed, as all of the output is in a relatively straight-forwards plain-text format.
