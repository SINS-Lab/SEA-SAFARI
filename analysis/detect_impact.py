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

safio = safari_input.SafariInput(args.input)

command = 'Sea-Safari.exe'
if platform.system() == 'Linux':
    command = './Sea-Safari'

sub = subprocess.run(command, shell=True)

xyz_in = args.input.replace('.input', '.xyz')

if args.colour:
    xyz_p.process_file(xyz_in, args.output, color=args.colour, load_vmd=True)
else:
    xyz_p.process_file(xyz_in, args.output, load_vmd=True)

