from xyz import XYZ
from xyz import XYZ_Single
import subprocess
import os
import platform
import time as T
import numpy as np
import argparse
import math

def process_velocities(xyz):
    for xyz_single in xyz.xyzs:
        num = xyz_single.number
        avg = 0
        p = 0
        px = 0
        py = 0
        pz = 0
        p_arr = []
        for i in range(num):
            if i == 0:
                p_arr.append(0)
                continue
            px = xyz_single.values[i][3]
            py = xyz_single.values[i][4]
            pz = xyz_single.values[i][5]
            p = math.sqrt(px*px + py*py + pz*pz)
            p_arr.append(p)
            avg = avg + p
        avg = avg / num
        for i in range(num):
            p = p_arr[i]
            # TODO better definition of "fast"
            if p > 2*avg + 1 and i != 0:
                xyz_single.atoms[i] = 'X'
    return

def process_nearby(xyz):
    for xyz_single in xyz.xyzs:
        num = xyz_single.number
        for i in range(num):
            if i != 0 and xyz_single.values[i][7] and xyz_single.values[i][8]:
                xyz_single.atoms[i] = 'X'
    return

def process(xyz, color=""):
    if color == "velocity":
        process_velocities(xyz)
    if color == "nearest":
        process_nearby(xyz)
    return xyz


def process_file(fileIn, fileOut=None, color="", load_vmd=False):
    """
    :param color: string describing the coloring pattern. "" for uniform coloring, "nearest" for coloring of
        nearest 10 atoms, and "momentum" to color the atoms that receive additional momentum differently
    :param load_vmd: boolean describing whether to open vmd upon completetion (true) or not (false)
    """
    if fileOut is None:
        fileOut = fileIn
    
    if not os.path.exists(fileOut):
        #Get the XYZ executable to handle smoothing.
        command = 'XYZ.exe {} {}'
        if platform.system() == 'Linux':
            command = './XYZ {} {}'
        command = command.format(fileIn, fileOut)
        subprocess.run(command, shell=True)
    
    xyz = XYZ()
    xyz.load(fileOut)

    if color != "":
        if color == "nearest":
            fileOut = fileOut[:-4] + "COLOREDNear.xyz"
        elif color == "velocity":
            fileOut = fileOut[:-4] + "COLOREDVel.xyz"
        #We need to edit the xyz to have different atom names
        xyz = process(xyz, color=color)
        xyz.save(fileOut)
    
    if load_vmd:
        # MAKE THE FILENAME INCLUDE DIRECTORY
        #Replace \ with / in filenames
        fileOut = fileOut.replace('\\','/')
        if color != "":
            commands = ["topo readvarxyz {}\n".format(fileOut)]
        else:
            commands = ["mol new {}\n".format(fileOut)]
        commands.append("mol modstyle 0 0 \"VDW\"")
        try:
            with open("commands.vmd", "w") as file:
                file.writelines(commands)
            subprocess.Popen(["vmd", "-e", "commands.vmd"])
        finally:
            T.sleep(5)
            os.remove("commands.vmd")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input file to use")
    parser.add_argument("-o", "--output", help="Output file to use")
    parser.add_argument("-c", "--colour", help="Colouring method to use")
    args = parser.parse_args()

    colour = ''
    if args.colour:
        colour = args.colour
    inputfile = args.input
    outputfile = args.output
    process_file(inputfile, outputfile, color=colour)
