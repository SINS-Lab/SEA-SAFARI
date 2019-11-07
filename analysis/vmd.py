import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input XYZ file")
parser.add_argument("-c", "--colour", help="Whether to not colour", action='store_true')
args = parser.parse_args()

commands = []

if not args.colour:
    command = "topo readvarxyz {}\n".format(args.input)
else:
    command = "mol new {}\n".format(args.input)

commands.append("color Display Background white\n")
commands.append("display depthcue off\n") # This fixes washed out colours on white background
commands.append("mol default style {VDW 1.0 30.0}\n")
commands.append(command)
commands.append("axes location off\n")
commands.append("rotate y to 0\n")
commands.append("rotate z to 0\n")
commands.append("rotate x to -30")

with open("commands.vmd", "w") as file:
    file.writelines(commands)
subprocess.run(["vmd", "-e", "commands.vmd"])