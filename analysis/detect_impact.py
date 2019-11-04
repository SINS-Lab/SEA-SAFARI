import safari_input
import subprocess
import argparse
import platform
import os
import time

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="SAFIO input file")
parser.add_argument("-o", "--output", help="Output file names")
parser.add_argument("-x", "--x_coord", help="x - impact coordinate")
parser.add_argument("-y", "--y_coord", help="y - impact coordinate")
parser.add_argument("-r", "--restricted", help="whether the xyz is only nearish particles", action='store_true')
parser.add_argument("-c", "--colour", help="Colour parameter for xyz")
args = parser.parse_args()

# Run a single shot safari for this run,
# assuming the input file 
# was already configured properly.
command = 'Sea-Safari.exe -i {} -o {} -s -x {} -y {}'

if args.restricted:
    command = 'Sea-Safari.exe -i {} -o {} -s -x {} -y {} -r'

#Change command accordingly for linux
if platform.system() == 'Linux':
    command = command.replace('Sea-Safari.exe', './Sea-Safari')

#Format the command
command = command.format(args.input, args.output, args.x_coord, args.y_coord)

subprocess.run(command, shell=True)

xyz_in = args.output + '.xyz'
fileOut = xyz_in

command = ''

#format input argument for XYZ processor
if args.colour:
    fileOut = fileOut.replace('.xyz', '_{}.xyz'.format(args.colour))
    command = 'XYZ.exe -i {} -o {} -c {}'
    command = command.format(xyz_in, fileOut, args.colour)
else:
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

