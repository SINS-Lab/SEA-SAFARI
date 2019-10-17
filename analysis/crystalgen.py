import helpers
import numpy as np 
import basisgen

def equals(site, old, tol=.5):
    dist = (site[0]-old[0])*(site[0]-old[0]) + \
           (site[1]-old[1])*(site[1]-old[1]) + \
           (site[2]-old[2])*(site[2]-old[2])
    return dist < tol

def clearDuplicates(crystal):
    d = 0
    for i in range(len(crystal)):
        old = crystal[i]
        n = i + 1
        for j in range(i+1, len(crystal)):
            site = crystal[n]
            n = n + 1
            if equals(site, old):
                d = d + 1
                n = n - 1
                crystal.remove(site)
            if n >= len(crystal):
                break
        if i + 1 >= len(crystal):
            break
    if d != 0:
        print(str(d)+' Dupes found.')
        #clearDuplicates(crystal)

def clearOutOfBounds(crystal, xmin, xmax, ymin, ymax):
    culled = []
    for i in range(len(crystal)):
        site = crystal[i]
        if site[0] < xmin or site[0] > xmax or site[1] < ymin or site[1] > ymax:
            continue
        culled.append(site)
    output = open('culled.input', 'w')
    output2 = open('crystal.input', 'w')
    output3 = open('crystal.xyz', 'w')
    output3.write(str(len(culled)-1)+'\n')
    for site in culled:
        output.write(str(round(site[0] - xmin,4))+'\t'+str(round(site[1] - ymin,4))+ '\t'+ str(round(site[2],4)) \
        +'\t1\n')
        output2.write(str(site[0])+'\t'+str(site[1])+ '\t'+ str(site[2]) \
        +'\t'+str(site[3])+'\n')
        output3.write('A\t'+str(site[0])+'\t'+str(site[1])+ '\t'+ str(site[2])+'\n')
    output.close()
    output2.close()
    output3.close()
    return culled


def contains(crystal, site, tol=.5):
    for old in crystal:
        if equals(site, old, tol):
            return True
    return False

def gen(_size, _dir, _axis, _basis, _n, zTop, zBottom):
    
    basisgen.gen(False, _size, _dir, _axis, _basis)
    
    file = open('basis.input', 'r')
    
    onlyTop = False
    
    if onlyTop:
        zBottom = -2.5*_size
    
    n = 0;
    
    # Repeat distance for crystal
    lattice = [0,0,0]
    
    dir = [0,0,1]
    axis = [0,0,1]
    
    # List of basis coordinates
    basis = []
    
    # Lists of crystal coordinates to export
    crystal = []
    
    for line in file:
        arr = line.split()
        if n==0 :
            lattice[0] = float(arr[0].replace('D0','').replace('d0',''))
            lattice[1] = float(arr[1].replace('D0','').replace('d0',''))
            lattice[2] = float(arr[2].replace('D0','').replace('d0',''))
            n = n + 1
            continue
        if n==1 :
            dir[0] = int(arr[0].replace('D0','').replace('d0',''))
            dir[1] = int(arr[1].replace('D0','').replace('d0',''))
            dir[2] = int(arr[2].replace('D0','').replace('d0',''))
            n = n + 1
            continue
        if n==2 :
            axis[0] = int(arr[0].replace('D0','').replace('d0',''))
            axis[1] = int(arr[1].replace('D0','').replace('d0',''))
            axis[2] = int(arr[2].replace('D0','').replace('d0',''))
            n = n + 1
            continue
        basis.append([float(arr[0]),float(arr[1]),float(arr[2])])
        n = n + 1
    file.close()
    
    
    im = np.array(axis)
    ix = np.array(dir)
    R = helpers.rotate(ix, im)
    R_inv = np.linalg.inv(R)
    
    ex = np.asarray(np.matmul(R_inv, np.array([1,0,0])))[0]*lattice[0]
    ey = np.asarray(np.matmul(R_inv, np.array([0,1,0])))[0]*lattice[1]
    ez = np.asarray(np.matmul(R_inv, np.array([0,0,1])))[0]*lattice[2]
    
    maxZ = -1e20
    maxZI = 0
    
    minZ = 1e20
    minZI = 0
    
    for i in range(len(basis)):
        if basis[i][2]>maxZ:
            maxZI = i
            maxZ = basis[i][2]
        if basis[i][2]<minZ:
            minZI = i
            minZ = basis[i][2]
    
    n = 0
    m = 0
    diff = [0,0,0]
    
    ns = -_n
    ne = _n

    for x in range(ns, ne):
        for y in range(ns, ne):
            for z in range(ns, ne):
                diff[0] = (ex[0] * x + ex[1] * y + ex[2] * z)
                diff[1] = (ey[0] * x + ey[1] * y + ey[2] * z)
                diff[2] = (ez[0] * x + ez[1] * y + ez[2] * z)
                # Check if entire cell fits
                if (diff[2] + basis[maxZI][2]) <= zTop:
                    for i in range(len(basis)):
                        if onlyTop and i!=maxZI:
                            continue
                        px = diff[0] + basis[i][0]
                        py = diff[1] + basis[i][1]
                        pz = diff[2] + basis[i][2]
                        if pz < zBottom:
                            continue
                        site = [px,py,pz,_basis.atoms[i]]
                        if not site in crystal:
                            crystal.append(site)
    print('Generated Points, Now Clearing Duplicates')
    n = len(crystal)
    print(str(n)+ ' Points to check')
    clearDuplicates(crystal)
    print('Duplicates Cleared')
    output = open('crystal.input', 'w')
    for site in crystal:
        output.write(str(site[0])+'\t'+str(site[1])+ '\t'+ str(site[2]) \
        +'\t'+str(site[3])+'\n')
    output.close()
    return crystal
    

if __name__ == '__main__':
    size = 4.0786
    dir = [0,0,1]
    axis = [1,1,1]
    atom = basisgen.Atom(196.967,79)
    crystal = gen(size, dir, axis, basisgen.fccBasis(atom), 10, 0.1, -1.8*size)
    n = 8.9
    clearOutOfBounds(crystal, -size * n, size * n, -size * n, size * n)