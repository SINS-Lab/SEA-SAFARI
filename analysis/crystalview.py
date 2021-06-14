#!/usr/bin/env python3

import argparse     # Parsing input arguments
import subprocess   # Runs VMD
import platform     # Linux vs Windows check

import os           # os.remove is used for .vmd file
import time         # sleeps before removing file

from mendeleev import element # Lookup atomic symbol from number

elements = {}

def find(number_str):
    if number_str in elements:
        return elements[number_str]
    num = int(float(number_str))
    A = element(num)
    sym = A.symbol
    elements[number_str] = sym
    return sym

def make_xyz(file_in, file_out):
    lines = []
    file = open(file_in, 'r')
    for line in file:
        args = line.split()
        if len(args)==5:
            lines.append(line.split())
    file.close()

    elements = {}


    xyz_lines = []
    # Write the number of particles
    xyz_lines.append('{}\n'.format(len(lines)))
    # Write the blank comment line
    xyz_lines.append('\n'.format(len(lines)))
    # Now we write each particle line
    for line in lines:
        sym = find(line[3])
        xyz_lines.append('{} {} {} {}\n'.format(sym, line[0], line[1], line[2]))
    file = open(file_out, 'w')
    file.writelines(xyz_lines)
    file.close()


if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input .crys or .crys_in file")
    args = parser.parse_args()

    file_out = args.input.replace('.crys_in', '').replace('.crys', '') + '.xyz'

    print(args.input)
    print(file_out)

    make_xyz(args.input, file_out)

    commands = []
    
    commands.append("color Display Background white\n")
    commands.append("display depthcue off\n") # This fixes washed out colours on white background
    commands.append("mol default style {VDW 1.0 10.0}\n")
    commands.append("mol new {}\n".format(file_out))
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