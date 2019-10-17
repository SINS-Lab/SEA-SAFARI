from mpl_toolkits.mplot3d import Axes3D  # This is needed for the plot.
import matplotlib.pyplot as plt
import numpy as np
import helpers

class Atom:
    def __init__(self, mass, charge):
        self.mass = mass
        self.charge = charge
        
    def __str__(self):
        return str(self.charge)+'\t'+str(self.mass)

class Basis:
    def __init__(self):
        self.atoms = []
        self.sites = []
    
    def append(self, atom, site):
        self.atoms.append(atom)
        self.sites.append(site)

def fccBasis(atom):
    basis = Basis()
    # Corners
    basis.append(atom, np.array([ 0.0, 0.0, 0.0]))
    basis.append(atom, np.array([ 0.0, 1.0, 0.0]))
    basis.append(atom, np.array([ 0.0, 0.0, 1.0]))
    basis.append(atom, np.array([ 0.0, 1.0, 1.0]))
    
    basis.append(atom, np.array([ 1.0, 0.0, 0.0]))
    basis.append(atom, np.array([ 1.0, 1.0, 0.0]))
    basis.append(atom, np.array([ 1.0, 0.0, 1.0]))
    basis.append(atom, np.array([ 1.0, 1.0, 1.0]))

    #atom = Atom(atom.mass/2, atom.charge/2)
    # Face Cemtres
    basis.append(atom, np.array([ 0.0, 0.5, 0.5]))
    basis.append(atom, np.array([ 1.0, 0.5, 0.5]))
    
    basis.append(atom, np.array([ 0.5, 0.0, 0.5]))
    basis.append(atom, np.array([ 0.5, 1.0, 0.5]))
    
    basis.append(atom, np.array([ 0.5, 0.5, 0.0]))
    basis.append(atom, np.array([ 0.5, 0.5, 1.0]))
    return basis

def bccBasis(atom):
    basis = Basis()
    # Corners
    basis.append(atom, np.array([ 0.0, 0.0, 0.0]))
    basis.append(atom, np.array([ 0.0, 1.0, 0.0]))
    basis.append(atom, np.array([ 0.0, 0.0, 1.0]))
    basis.append(atom, np.array([ 0.0, 1.0, 1.0]))
    
    basis.append(atom, np.array([ 1.0, 0.0, 0.0]))
    basis.append(atom, np.array([ 1.0, 1.0, 0.0]))
    basis.append(atom, np.array([ 1.0, 0.0, 1.0]))
    basis.append(atom, np.array([ 1.0, 1.0, 1.0]))
    
    # Body Cemtres
    basis.append(atom, np.array([ 0.5, 0.5, 0.5]))
    return basis

def gen(show, size, dir, axis, basis):
    
    for point in basis.sites:
        np.multiply(point, size, out=point)
    
    original = np.array(basis.sites)
    cube = []
    
    lattice = [size,size,size]
    
    im = np.array(axis)
    ix = np.array(dir)
    R = helpers.rotate(ix, im)
    
    file = open('basis.input', 'w')
    
    file.write(str(lattice[0])+' '+str(lattice[1])+' '+str(lattice[2])+'\n')
    file.write(str(dir[0])+' '+str(dir[1])+' '+str(dir[2])+'\n')
    file.write(str(axis[0])+' '+str(axis[1])+' '+str(axis[2])+'\n')
    
    for i in range(len(original)):
        v = original[i]
        p = np.asarray(np.asarray((np.matmul(R, v))[0])[0]);
        cube.append(p)
    
    if show:
        fig = plt.figure()
    if show:
        ax = fig.add_subplot(111, projection='3d')
    
    for p in original:
        if show:
            ax.scatter(p[0], p[1], p[2], c='g')
    n = 0
    for p in cube:
        if show:
            ax.scatter(p[0], p[1], p[2], c='r')
        file.write(str(p[0])+'\t'+str(p[1])+'\t'+str(p[2])+'\t'+str(basis.atoms[n])+'\n')
        n = n+1
    
        
    file.close()
    if show:
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        
        ax.set_xlim(-size, size)
        ax.set_ylim(-size, size)
        ax.set_zlim(-size, size)
        plt.show()


if __name__ == '__main__':
    size = 4.08
    dir = [0,0,1]
    axis = [0,0,1]
    atom = Atom(196.966570,79)
    gen(True, size, dir, axis, fccBasis(atom))