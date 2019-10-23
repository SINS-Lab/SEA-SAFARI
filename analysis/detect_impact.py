import safari_input
import subprocess
import xyz_postprocess as xyz_p
import argparse
import platform
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="SAFIO input file")
parser.add_argument("-o", "--output", help="Output file names")
parser.add_argument("-c", "--colour", help="Colour parameter for xyz")
args = parser.parse_args()

print("test a")
if not os.path.exists(args.output):
    command = 'Sea-Safari.exe'
    if platform.system() == 'Linux':
        command = './Sea-Safari'
    #subprocess.run(command, shell=True)
print("test b")

xyz_in = args.input.replace('.input', '.xyz')

if args.colour:
    xyz_p.process_file(xyz_in, args.output, color=args.colour, load_vmd=True)
else:
    command = 'XYZ.exe {} {}';
    if platform.system() == 'Linux':
        command = './XYZ {} {}'
    command = command.format(xyz_in, args.output);
    subprocess.run(command, shell=True)
    fileOut = args.output
    # MAKE THE FILENAME INCLUDE DIRECTORY
    #Replace \ with / in filenames
    fileOut = fileOut.replace('\\','/')
    commands = ["mol new {}\n".format(fileOut)]
    commands.append("mol modstyle 0 0 \"VDW\"")
    try:
        with open("commands.vmd", "w") as file:
            file.writelines(commands)
        subprocess.Popen(["vmd", "-e", "commands.vmd"])
    finally:
        T.sleep(5)
        os.remove("commands.vmd")

