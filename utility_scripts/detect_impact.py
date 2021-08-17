import subprocess   # used to run safari and vmd
import argparse     # parses command line arguments
import platform     # Linux vs Windows check
import os           # os.remove is used for .vmd file
import time         # sleeps before removing file

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="SAFIO input file")
parser.add_argument("-o", "--output", help="Output file names")
parser.add_argument("-x", "--x_coord", help="x - impact coordinate")
parser.add_argument("-y", "--y_coord", help="y - impact coordinate")
parser.add_argument("-r", "--restricted", help="whether the xyz is only nearish particles", action='store_true')
parser.add_argument("-c", "--colour", help="Colour parameter for xyz")
parser.add_argument("-s", "--seed", help="Ion index/thermal seed for replicating the thermalization of the surface")
args = parser.parse_args()

# Run a single shot safari for this run,
# assuming the input file 
# was already configured properly.
command = 'Sea-Safari.exe -i {} -o {} -s -x {} -y {} --seed {}'

if args.restricted:
    command = 'Sea-Safari.exe -i {} -o {} -s -x {} -y {} --seed {} -r'

#Change command accordingly for linux
if platform.system() == 'Linux':
    command = command.replace('Sea-Safari.exe', './Sea-Safari')

#Format the command
command = command.format(args.input, args.output, args.x_coord, args.y_coord, args.seed)

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

commands = []

commands.append("color Display Background white\n")
commands.append("display depthcue off\n") # This fixes washed out colours on white background
commands.append("mol default style {VDW 1.0 10.0}\n")
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

