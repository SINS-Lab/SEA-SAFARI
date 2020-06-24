#!/usr/bin/env python3

import subprocess   # used to run safari and vmd
import argparse     # parses command line arguments
import platform     # Linux vs Windows check
import os           # os.remove is used for .vmd file
import time         # sleeps before removing file

def test_montecarlo(template):
    command = template.format('sample', 'tests/sample')
    print("Running Test of Montecarlo\n\n")
    subprocess.run(command, shell=True)
    
def test_adaptive_grid(template):
    command = template.format('sample_ag', 'tests/sample_ag')
    print("\n\nRunning Test of Adaptive Grid\n\n")
    subprocess.run(command, shell=True)

def test_chainscat(template):
    command = template.format('sample_chain', 'tests/sample_chain')
    print("\n\nRunning Test of Chainscat\n\n")
    subprocess.run(command, shell=True)

def test_single_shot(template):
    command = template.format('cs_cu -s -r', 'tests/cs_cu_ss')
    print("\n\nRunning Test of Restricted Single Shot\n\n")
    subprocess.run(command, shell=True)

    xyz_in = 'tests/cs_cu_ss.xyz'
    fileOut = 'tests/cs_cu_ss_proc.xyz'

    command = ''

    #format input argument for XYZ processor
    command = 'XYZ.exe -i {} -o {}'
    command = command.format(xyz_in, fileOut)

    #Change command accordingly for linux
    if platform.system() == 'Linux':
        command = command.replace('XYZ.exe', './XYZ')

    #Run XYZ processor and wait for it to finish.
    subprocess.run(command, shell=True)

    # MAKE THE FILENAME INCLUDE DIRECTORY
    #Replace \ with / in filenames
    fileOut = fileOut.replace('\\','/')

    command = "mol new {}\n".format(fileOut)

    commands = []

    commands.append("color Display Background white\n")
    commands.append("display depthcue off\n") # This fixes washed out colours on white background
    commands.append("mol default style {Licorice 1.1 20 12}\n")
    commands.append("mol default color MASS\n")
    commands.append(command)
    commands.append("display rendermode GLSL\n")
    commands.append("display update\n")
    commands.append("display update ui\n")
    try:
        with open("commands.vmd", "w") as file:
            file.writelines(commands)
        subprocess.run(["vmd", "-e", "commands.vmd"])
    finally:
        time.sleep(5)
        os.remove("commands.vmd")


if __name__ == '__main__':
    # Run a single shot safari for this run,
    # assuming the input file 
    # was already configured properly.
    template = 'Sea-Safari.exe -i {} -o {}'

    #Change command accordingly for linux
    if platform.system() == 'Linux':
        template = template.replace('Sea-Safari.exe', './Sea-Safari')

    test_montecarlo(template)
    test_adaptive_grid(template)
    test_chainscat(template)