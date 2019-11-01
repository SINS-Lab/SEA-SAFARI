import safari_input
import subprocess
import argparse
import platform
import os
import time

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="SAFIO input file")
parser.add_argument("-o", "--output", help="Output file names")
parser.add_argument("-c", "--colour", help="Colour parameter for xyz")
args = parser.parse_args()

# Run a single shot safari for this run,
# assuming the input file 
# was already configured properly.
command = 'Sea-Safari.exe'
if platform.system() == 'Linux':
    command = './Sea-Safari'
subprocess.run(command, shell=True)

xyz_in = args.input.replace('.input', '.xyz')

command = ''

#format input argument for XYZ processor
if args.colour:
    command = 'XYZ.exe -i {} -o {} -c {}'
    command = command.format(xyz_in, args.output, args.colour)
else:
    command = 'XYZ.exe -i {} -o {}'
    command = command.format(xyz_in, args.output)

#Change command accordingly for linux
if platform.system() == 'Linux':
    command = command.replace('XYZ.exe', './XYZ')

#Run XYZ processor and wait for it to finish.
subprocess.run(command, shell=True)
fileOut = args.output

# MAKE THE FILENAME INCLUDE DIRECTORY
#Replace \ with / in filenames
fileOut = fileOut.replace('\\','/')

if args.colour:
    command = "topo readvarxyz {}\n".format(fileOut)
else:
    command = "mol new {}\n".format(fileOut)

commands = [command]
commands.append("mol modstyle 0 0 \"VDW\"")
try:
    with open("commands.vmd", "w") as file:
        file.writelines(commands)
    subprocess.run(["vmd", "-e", "commands.vmd"])
finally:
    time.sleep(5)
    os.remove("commands.vmd")

