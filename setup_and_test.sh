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

CopyFiles() {

    mkdir -p "$1"

    cp SAFARI/bin/Release/Sea-Safari "$1/Sea-Safari"
    cp SAFARI/bin/Release/sample.input "$1/sample.input"
    cp SAFARI/analysis/XYZ "$1/XYZ"

    cp SAFARI/analysis/safari_input.py "$1/safari_input.py"
    cp SAFARI/analysis/run_all.py "$1/run_all.py"
    cp SAFARI/analysis/run_generator.py "$1/run_generator.py"

    cp SAFARI/analysis/detect_processor.py "$1/detect_processor.py"
    cp SAFARI/analysis/detect_impact.py "$1/detect_impact.py"
    cp SAFARI/analysis/detect_gui.py "$1/detect_gui.py"
}

if [[ -d "./analysis" && -d "./bin" && -d "./src" ]]
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

cd SAFARI
make

# Back out of SAFARI to make the rest of them
cd ..

echo ""
echo "Copying required files to runs and tests"
echo ""

CopyFiles runs
CopyFiles tests

# Copy some extra files for tests as well
cp SAFARI/bin/Release/run_tests.py tests/run_tests.py
cp SAFARI/bin/Release/sample_ag.input tests/sample_ag.input
cp SAFARI/bin/Release/sample_chain.input tests/sample_chain.input

cd tests

# tests folder inside tests for data output
mkdir -p tests

echo ""
echo "Running standardized test set."
echo ""

# Run the tests
python3 run_tests.py

# Return back up to root directory
cd ..


echo ""
echo ""
echo "Finished install and tests, there should be output files in tests/tests, if not, something went wrong!"
echo ""
echo ""
