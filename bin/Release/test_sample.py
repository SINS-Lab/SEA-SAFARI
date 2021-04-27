#!/usr/bin/env python3

import subprocess   # used to run safari and vmd
import argparse     # parses command line arguments
import platform     # Linux vs Windows check
import os           # os.remove is used for .vmd file
import time         # sleeps before removing file

def test_montecarlo(template):
    command = template.format('test_montecarlo', 'tests/test_montecarlo')
    print("Running Test of Montecarlo\n\n")
    subprocess.run(command, shell=True)
    
def test_adaptive_grid(template):
    command = template.format('test_adaptive_grid', 'tests/test_adaptive_grid')
    print("\n\nRunning Test of Adaptive Grid\n\n")
    subprocess.run(command, shell=True)

def test_chainscat(template):
    command = template.format('test_chain', 'tests/test_chain')
    print("\n\nRunning Test of Chainscat\n\n")
    subprocess.run(command, shell=True)

if __name__ == '__main__':
    # Run a single shot safari for this run,
    # assuming the input file 
    # was already configured properly.
    template = 'Sea-Safari.exe -i {} -o {}'

    #Change command accordingly for linux
    if platform.system() == 'Linux':
        template = template.replace('Sea-Safari.exe', './Sea-Safari')

    if not os.path.exists('./tests'):
        os.makedirs('./tests')

    test_montecarlo(template)
    test_adaptive_grid(template)
    test_chainscat(template)