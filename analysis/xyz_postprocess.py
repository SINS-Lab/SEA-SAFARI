from xyz import XYZ
from xyz import XYZ_Single
import subprocess
import os
import platform
import time as T
import numpy as np
import math


class Particle:

    nextid = 1

    def __init__(self):
        self.atom = None
        self.pos = [0, 0, 0]
        self.momentum = [0, 0, 0]
        self.velocity = [0, 0, 0]
        self.mass = 0
        #ion won't have id in array, so this sets us as ion.
        self.id = -1

    def fromXYZ(self, atom, value):
        self.atom = atom
        self.pos = [value[0], value[1], value[2]]
        self.momentum = [value[3], value[4], value[5]]
        self.mass = value[6]
        self.velocity = [value[3]/value[6],
                         value[4]/value[6],
                         value[5]/value[6]]
        #Everyone except ion will load this in, this is >=0
        if len(value)>7:
            self.id = value[7];

def xyzFromParticles(time, pset):
    xyz_single = XYZ_Single()
    xyz_single.number = len(pset)
    xyz_single.comment = str(time)+'\n'
    for particle in pset:
        xyz_single.atoms.append(particle.atom)
        value = []
        value.append(particle.pos[0])
        value.append(particle.pos[1])
        value.append(particle.pos[2])
        value.append(particle.momentum[0])
        value.append(particle.momentum[1])
        value.append(particle.momentum[2])
        value.append(particle.mass)
        xyz_single.values.append(value)
    return xyz_single

def get_velocity_parameters(particles):
    helper_sum = 0
    for state in particles:
        for particle in state:
            if particle.id != -1:
                helper_sum += np.linalg.norm(particle.velocity)
    average = helper_sum/(len(particles)*len(particles[0]) - 1)
    helper_sum = 0
    for state in particles:
        for particle in state:
            if particle.id != -1:
                helper_sum += (average - np.linalg.norm(particle.velocity))**2
    variation_average = helper_sum/(len(particles)*len(particles[0]) - 1)
    print(math.sqrt(variation_average))
    return average, math.sqrt(variation_average)

def process(xyz, color=""):
    # Array of elapsed times
    times = []
    # This is an array of arrays of Particles, which correspond to times
    particles = []
    n = 0
    for xyz_single in xyz.xyzs:
        # Particle array for this timestep
        pset = []
        comment = xyz_single.comment.split()
        times.append(float(comment[0]))
        for i in range(xyz_single.number):
            particle = Particle()
            # Load in the particle from xyz, the ion will stay as id == -1
            # Everyone else will find their lattice site ID and use that.
            particle.fromXYZ(xyz_single.atoms[i], xyz_single.values[i])
            pset.append(particle)
            n = n + 1
        particles.append(pset)

    if color == "velocity":
        mean, std = get_velocity_parameters(particles)
        for state in particles:
            for atom in state:
                if np.linalg.norm(atom.velocity) - mean > std and atom.id != -1:
                    atom.atom = "X"
    xyz = XYZ()
    for i in range(len(times)):
        time = times[i]
        pset = particles[i]
        xyz_single = xyzFromParticles(time, pset)
        xyz.xyzs.append(xyz_single)
    return xyz


def process_file(fileIn, fileOut=None, color="", load_vmd=False):
    """
    :param color: string describing the coloring pattern. "" for uniform coloring, "nearest" for coloring of
        nearest 10 atoms, and "momentum" to color the atoms that receive additional momentum differently
    :param load_vmd: boolean describing whether to open vmd upon completetion (true) or not (false)
    """
    if fileOut is None:
        fileOut = fileIn
        
    #Get the XYZ executable to handle smoothing.
    command = 'XYZ.exe {} {}';
    if platform.system() == 'Linux':
        command = './XYZ {} {}'
    command = command.format(fileIn, fileOut);
    subprocess.run(command, shell=True)
    
    xyz = XYZ()
    xyz.load(fileOut)
    if color != "":
        #We need to edit the xyz to have different atom names
        xyz = process(xyz, color=color)
        xyz.save(fileOut)
    if color == "nearest":
        fileOut = fileOut[:-4] + "COLOREDNear.xyz"
    elif color == "velocity":
        fileOut = fileOut[:-4] + "COLOREDVel.xyz"
    
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
    process_file('sample.xyz', 'sample_mod.xyz')
