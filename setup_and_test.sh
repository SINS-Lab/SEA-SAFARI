#!/bin/bash
#
# This will clone the github repository, then make SAFARI,
# Then it will copy the executable to two locations:
#
# runs and tests
#
# in tests, it will then run the standardized test set for
# SAFARI, it confirm that things are functional
#
#

CopyExecutables() {

    mkdir -p "$1"

    cp SAFARI/bin/Release/Sea-Safari "$1/Sea-Safari"
    cp SAFARI/utility_scripts/XYZ "$1/XYZ"
}

CopyRunScripts() {

    mkdir -p "$1"

    cp SAFARI/utility_scripts/safari_input.py "$1/safari_input.py"
    cp SAFARI/utility_scripts/run_all.py "$1/run_all.py"
    cp SAFARI/utility_scripts/run_generator.py "$1/run_generator.py"
    
    cp SAFARI/bin/Release/sample.input "$1/template.input"
}

if [[ -d "./utility_scripts" && -d "./bin" && -d "./src" ]]
then
    echo "Running from source!"
    cd ..
else
    if [ -d "./SAFARI" ] 
    then
        echo "Cleaning up Previous installation"
        rm -rf "./SAFARI"
    else
        echo "No Previous installation."
    fi
    git clone https://github.com/SINS-Lab/SEA-SAFARI.git SAFARI
fi

# Clone the SAFARI-ANALYSIS as well
if [ -d "./SAFARI-ANALYSIS" ] 
then
    echo "Cleaning up Previous installation"
    rm -rf "./SAFARI-ANALYSIS"
else
    echo "No Previous installation."
fi
git clone https://github.com/SINS-Lab/SAFARI-ANALYSIS.git SAFARI-ANALYSIS

cd SAFARI
make

# Back out of SAFARI to make the rest of them
cd ..

echo ""
echo "Copying required files to runs and tests"
echo ""

CopyExecutables runs
CopyRunScripts runs

CopyExecutables SAFARI-ANALYSIS
CopyExecutables tests

# Copy some extra files for tests as well
cp SAFARI/bin/Release/test_sample.py tests/test_sample.py
cp SAFARI/bin/Release/test_montecarlo.input tests/test_montecarlo.input
cp SAFARI/bin/Release/test_adaptive_grid.input tests/test_adaptive_grid.input
cp SAFARI/bin/Release/test_chain.input tests/test_chain.input

cd tests

echo ""
echo "Running standardized test set."
echo ""

# Run the tests
./test_sample.py

# Return back up to root directory
cd ..


echo ""
echo ""
echo "Finished build and tests, there should be output files in tests/tests, if not, something went wrong!"
echo ""
echo ""
