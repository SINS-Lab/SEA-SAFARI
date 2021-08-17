import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input XYZ file")
parser.add_argument("-z", "--z", help="Z rotation for the lattice")
parser.add_argument("-c", "--colour", help="Whether to not colour", action='store_true')
args = parser.parse_args()

commands = []

z = 0
if args.z:
    z = float(args.z)

if not args.colour:
    command = "topo readvarxyz {}\n".format(args.input)
else:
    command = "mol new {}\n".format(args.input)

commands.append("color Display Background white\n")
commands.append("display depthcue off\n") # This fixes washed out colours on white background
commands.append("mol default style {VDW 1.0 10.0}\n")
commands.append(command)
commands.append("mol modcolor 0 0 PosZ\n")
commands.append("axes location off\n")
commands.append("scale to 0.0125\n")
commands.append("rotate y to 0\n")
commands.append("rotate z to {}\n".format(z))
commands.append("rotate x to 0\n")
commands.append("display resize 2000 2000\n")
commands.append("display rendermode GLSL\n")
commands.append("display update\n")
commands.append("display update ui\n")
commands.append("display reposition 2000 2000\n")
commands.append("wait 5\n")
commands.append("render TachyonInternal {}".format(args.input.replace('.xyz', '.ppm')))

with open("commands.vmd", "w") as file:
    file.writelines(commands)
subprocess.run(["vmd", "-e", "commands.vmd"])