from scipy.stats import maxwell
import scipy.linalg as linal
import numpy as np
import math
import crystalgen
import basisgen
import sys

class Particles:
    def __init__(self, *args, **kwargs):
        self.positions = np.zeros((0,3))
        self.momenta = np.zeros((0,3))
        self.masses = np.ones((1))
        self.charges = np.ones((1))
        self.springs = np.ones((0,3))
        self.r0 = np.zeros((0,3))
        self.coupling = True
        self.steps = True
        self.latticeMult = 1
        self.couplingMult = 1
    
    def T(self):
        if not self.steps:
            return 0
        p = (linal.norm(self.momenta, axis=1,keepdims=True))
        T = p * p / (2 * self.masses)
        return np.sum(T) / len(self.masses)

    def p_dot_coupled(self, dt):
        k = self.couplingMult
        n = len(self.positions)
        r0 = np.empty((n,3))
        zero = np.array([0,0,0])
        one = np.array([1,1,1])
        # Instead, this should use a matrix to do the same thing.
        # Need to see how to go about making the NxNx3 matrix that
        # this loop emulates.
        for i in range(n):
            r0[:] = self.positions[i]
            r = self.positions - r0
            r[i] = one
            s = (linal.norm(r, axis=1,keepdims=True))
            r = r/(s*s*s)
            dp = -k * r
            dp[i] = zero
            self.momenta = self.momenta + dp * dt


    def r_dot(self):
        p = self.momenta
        m = self.masses
        dr = p/m
        return dr

    def p_dot(self):
        r = self.positions - self.r0
        k = self.springs * self.latticeMult
        dp = -np.multiply(k, r)
        return dp

    def step(self, dt):
        if not self.steps:
            return

        self.momenta = self.momenta + self.p_dot() * dt
        self.positions = self.positions + self.r_dot() * dt
        if self.coupling:
            self.p_dot_coupled(dt)

    def randomizePositions(self, scale=0.1):
        if not self.steps:
            self.positions = self.r0
            return
        n = len(self.r0)
        rx = maxwell.rvs(size=n)
        ry = maxwell.rvs(size=n)
        rz = maxwell.rvs(size=n)
        r = np.array([rx,ry,rz])
        sign = lambda x: math.copysign(1, x)
        r = np.transpose(r)
        self.positions = self.positions + r * scale

    def genMomenta(self, masses, scale=0.01):
        if not self.steps:
            return
        n = len(masses)
        rx = maxwell.rvs(size=n)
        ry = maxwell.rvs(size=n)
        rz = maxwell.rvs(size=n)
        r = np.array([rx,ry,rz])
        sign = lambda x: math.copysign(1, x)
        r = np.transpose(r)
        for i in range(len(r)):
            x = np.random.random_sample() - 0.5
            y = np.random.random_sample() - 0.5
            z = np.random.random_sample() - 0.5
            r[i][0] = r[i][0]*sign(x)
            r[i][1] = r[i][1]*sign(y)
            r[i][2] = r[i][2]*sign(z)
        r = np.multiply(masses, r)
        r = np.multiply(r, scale)
        return r

    def save(self, inputfile):
        lines = np.hstack((self.positions, self.charges, self.masses, self.momenta, self.springs, self.r0 ));
        file = open(inputfile, 'w')

        for line in lines:
            file.write( str(line[0])+'\t'+str(line[1])+'\t'+str(line[2])+'\t' \
                       +str(line[3])+'\t' \
                       +str(line[4])+'\t' \
                       +str(line[5])+'\t'+str(line[6])+'\t'+str(line[7])+'\t' \
                       +str(line[8])+'\t'+str(line[9])+'\t'+str(line[10])+'\t'\
                       +str(line[11])+'\t'+str(line[12])+'\t'+str(line[13])+'\n')
        file.close()
        return

    def load(self, inputfile):
        file = open(inputfile, 'r')
        positions = []
        momenta = []
        masses = []
        charges = []
        springs = []
        r0 = []
        hasMomenta = False
        hasSprings = False
        for line in file:
            arr = line.split()
            xs = float(arr[0])
            ys = float(arr[1])
            zs = float(arr[2])

            positions.append([xs,ys,zs])

            charges.append([float(arr[3])])
            
            if len(arr) >= 5:
                masses.append([float(arr[4])])
            else:
                masses.append([1])
            if len(arr) >= 8:
                hasMomenta = True
                px = float(arr[5])
                py = float(arr[6])
                pz = float(arr[7])
                momenta.append([px,py,pz])
            if len(arr) >= 11:
                x = float(arr[8])
                y = float(arr[9])
                z = float(arr[10])
                springs.append([x,y,z])
            else:
                springs.append([100,100,100])
            if len(arr) >= 14:
                x = float(arr[11])
                y = float(arr[12])
                z = float(arr[13])
                r0.append([x,y,z])
            else:
                r0.append([xs,ys,zs])
        file.close()
        self.positions = np.array(positions)
        self.r0 = np.array(r0)
        self.masses = np.array(masses)
        self.charges = np.array(charges)
        self.springs = np.array(springs)
        self.momenta = np.array(momenta)
        if not hasMomenta:
            self.momenta = self.genMomenta(masses)
            self.randomizePositions()
        return

if __name__ == "__main__" :
    particles = Particles()
    size = 4.08
    dir = [0,0,1]
    axis = [7,8,8]
    atom = basisgen.Atom(196.966570,79)
   # crystalgen.gen(size, dir, axis, basisgen.fccBasis(atom), 10, 0.1, -3.5*size)
    particles.load('crystal.input')
    
    particles.step(.01)
    for i in range(100):
        particles.step(.01)
        if i%20==0:
            print(particles.T())
            particles.save('crystal.input')
        print(".", end='')
        if(i%10==9):
            print('')
    print(particles.positions)
    print(particles.momenta)
    
    particles.save('crystal.input')