import os
import shutil
import time

# A bunch of functions for fortran IO
def parseVar(input):
    string = str(input)
    if string.endswith('d0') or string.endswith('D0'):
        var = string.replace('D0','').replace('d0','')
        return float(var)
    if string == 't' or string == 'T':
        return True
    if string == 'f' or string == 'F':
        return False
    if isInt(string):
        return int(string)
    if isFloat(string):
        return float(string)
    return string

def parseLine(input):
    vars = input.split()
    ret = []
    for var in vars:
        ret.append(parseVar(var))
    return ret

def toStr(input):
    if isinstance(input, bool):
        return 't' if input else 'f'
    return str(input)

def serializeArr(input):
    output = ''
    n = 0
    for var in input:
        if n == 0:
            output = toStr(var)
        else:
            output = output + ' ' + toStr(var)
        n = n + 1
    return output

def serialize(*input):
    return serializeArr(input)

def isInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def isFloat(s):
    try: 
        float(s)
        return True
    except ValueError:
        return False

class SafariInput:

    def __init__(self, fileIn):
        self.fileIn = fileIn
        
        # Inialize All variables
        self.E0 = 1625.0
        self.THETA0 = 55.0
        self.PHI0 = 45.0
        self.MASS = 22.989
        self.SYMION = 'Na'

        self.EMIN = 0.5
        self.EMAX = 1625.0
        self.ESIZE = 16.0
        self.ASIZE = 5.0

        self.NDTECT = 1
        self.DTECTPAR = [55.0, 45.0, 3.0, 0.0]
        
        self.DELLOW = 1e-8
        self.DELT0 = 10.0
        self.DEMAX = 0.3
        self.DEMIN = 0.0
        self.ABSERR = 5e-6
        self.NPART = 45
        self.RECOIL = True
        self.Z1 = 5.0
        self.MAX_STEPS = 30000
        self.RRMIN = 0.0
        self.RRSTEP = 2.5e-3
        self.ZMIN = 0.0
        self.ZSTEP = 3e-4
        self.MAXDIV = 10
        self.MINDIV = 2

        self.NWRITX = 666
        self.NWRITY = 666

        self.RAX = 4
        self.RAY = 4

        self.NPAR = 4
        self.IPOT = 1
        self.POTPAR = [9125.84, 4.44, 66703.59, 11.24]
        self.NIMPAR = 2
        self.IIMPOT = 1
        self.PIMPAR = [1.26, 2.0]
        self.TEMP = 0.0
        self.SEED = 0.9436337324
        self.Ion_Index = 1
        self.IMAGE = True
        self.SENRGY = -1.5
        self.BDIST = 8.18

        self.AX = 4.09
        self.AY = 4.09
        self.AZ = 4.09

        self.NBASIS = 1
        self.BASIS = [[0.0, 0.0, 0.0, 1]]
        self.NTYPES = 1

        self.ATOMS = [[107.87, 47, 'Ag']]
        self.SPRINGS = [[5.0, 5.0, 5.0]]
        self.CORR = True
        self.ATOMK = 0
        self.RNEIGH = 8.36405
        
        #No idea on good defaults for these
        self.NUMCHA = 10000
        self.XSTART = 0
        self.XSTEP = 0.1
        self.XSTOP = 10
        self.YSTART = 0
        self.YSTEP = 0.1
        self.YSTOP = 10

        self.NBZ = 1
        self.TOL = 0
        self.NBG = 1
        self.GTOL = 0
        self.GMAX = [0]
        self.NG = [0]
        self.NZ = [0]
        self.ZMAX = [0]

        self.face = [0,0,1]
        self.load_crystal = False
        self.loaded_face = [0,0,1]

        self.F_a = 0
        self.F_b = 0

        self.load()
        # Instead here we should check for a default file, and use that
        # if we do not have the requested input file.
        self.fileIn = fileIn.replace('.input', '_mod.input')
        # This copies stuff from the original input to the modified one.
        self.save()
        return
    
    def genInputFile(self, fileIn=None):
        if fileIn is None:
            fileIn = self.fileIn
            fileIn = time.strftime("%Y%m%d_%H%M%S") +'.input'
        self.save(fileIn)
        saf_file = open('safari.input', 'w')
        saf_file.write(fileIn.replace('.input', ''))
        saf_file.close()
        
    def load(self):
        if not os.path.exists(self.fileIn):
            return

        file = open(self.fileIn, 'r')
        # This starts at 1, so it matches the file lines.
        n = 1
        # This is used for offsetting the lines depending
        # on whether anything depends on entries before.
        o = 0

        for line in file:
            # The variable names here are identical to 
            # the ones used in the fortran source.

            # This is a comment in the file
            if line.startswith("#"):
                continue
            
            args = parseLine(line)
            # Number of arguments, used for padding arrays with 0
            num = len(args)
            # Ignore blank lines.
            if num == 0:
                continue

            # E0 THETA0 PHI0 MASS SYMION
            if n == 1:
                self.E0 = args[0]
                self.THETA0 = args[1]
                self.PHI0 = args[2]
                self.MASS = args[3]
                self.SYMION = args[4]
                # Ensure phi is in correct range
                while self.PHI0 > 180:
                    self.PHI0 -= 360
                while self.PHI0 < -180:
                    self.PHI0 += 360
                
            #EMIN,EMAX,ESIZE,ASIZE
            if n == 2:
                self.EMIN = args[0]
                self.EMAX = args[1]
                self.ESIZE = args[2]
                self.ASIZE = args[3]
            # Detector type
            if n == 3:
                self.NDTECT = args[0]
            # Detector Params, 4 of them
            if n == 4:
                self.DTECTPAR = [args[0], args[1], args[2], args[3]]
            # min and max time steps
            if n == 5:
                self.DELLOW = args[0]
                self.DELT0 = args[1]
            # DEMAX,DEMIN,ABSERR
            if n == 6:
                self.DEMAX = args[0]
                self.DEMIN = args[1]
                self.ABSERR = args[2]
            # NPART - ATOMS ON SURFACE
            if n == 7:
                self.NPART = args[0]
            # Recoil
            if n == 8:
                self.RECOIL = args[0]
            # Interaction Height
            if n == 9:
                self.Z1 = args[0]
            # Table Entries
            if n == 10:
                self.MAX_STEPS = args[0]
            # R min and R step
            if n == 11:
                self.RRMIN = args[0]
                self.RRSTEP = args[1]
            # Z min and Z step
            if n == 12:
                self.ZMIN = args[0]
                self.ZSTEP = args[1]
            # Min and Max div, Here we check if o needs incrementing.
            if n == 13:
                if o <= 0:
                    self.MAXDIV = args[0]
                    self.MINDIV = args[1]
                    # This will be decremented after leaving this if
                    o = 4
                else:
                    if o == 3:
                        self.NUMCHA = args[0]
                    if o == 2:
                        self.XSTART = args[0]
                        self.XSTEP = args[1]
                        self.XSTOP = args[2]
                    if o == 1:
                        self.YSTART = args[0]
                        self.YSTEP = args[1]
                        self.YSTOP = args[2]
            # NWRITX, NWRITY
            if n == 14:
                self.NWRITX = args[0]
                self.NWRITY = args[1]
            # RAX and RAY
            if n == 15:
                self.RAX = args[0]
                self.RAY = args[1]
            # NPAR and IPOT
            if n == 16:
                self.NPAR = args[0]
                self.IPOT = args[1]
            # Parameters for IPOT
            if n == 17:
                self.POTPAR = []
                for i in range(self.NPAR):
                    if i<num:
                        self.POTPAR.append(args[i])
                    else:
                        self.POTPAR.append(0)
            # NIMPAR and IIMPOT
            if n == 18:
                self.NIMPAR = args[0]
                self.IIMPOT = args[1]
            # Parameters for IIMPOT, as well as other stuff?
            if n == 19:
                if o <= 0:
                    self.PIMPAR = []
                    for i in range(self.NIMPAR):
                        if i<num:
                            self.PIMPAR.append(args[i])
                        else:
                            self.PIMPAR.append(0)
                    if self.IIMPOT == 2:
                        # This will be decremented after leaving this if
                        o = 7
                elif o == 6:
                    self.NBZ = args[0]
                    self.TOL = args[1]
                elif o == 5:
                    self.ZMAX = []
                    for i in range(self.NBZ):
                        if i<num:
                            self.ZMAX.append(args[i])
                        else:
                            self.ZMAX.append(0)
                elif o == 4:
                    self.NZ = []
                    for i in range(self.NBZ):
                        if i<num:
                            self.NZ.append(args[i])
                        else:
                            self.NZ.append(0)
                elif o == 3:
                    self.NBG = args[0]
                    self.GTOL = args[1]
                elif o == 2:
                    self.GMAX = []
                    for i in range(self.NBG):
                        if i<num:
                            self.GMAX.append(args[i])
                        else:
                            self.GMAX.append(0)
                elif o == 1:
                    self.NG = []
                    for i in range(self.NBG):
                        if i<num:
                            self.NG.append(args[i])
                        else:
                            self.NG.append(0)
            # Temp, Seed, Ion_Index
            if n == 20:
                self.TEMP = args[0]
                self.SEED = args[1]
                self.Ion_Index = args[2]
            # Use Image Charge
            if n == 21:
                self.IMAGE = args[0]
            # SENRGY and BDIST
            if n == 22:
                self.SENRGY = args[0]
                self.BDIST = args[1]
            # AX and AY
            if n == 23:
                self.AX = args[0]
                self.AY = args[1]
                self.AZ = args[2]
            # Load the basis cell
            if n == 24:
                # How many atoms in the basis
                if o <= 0:
                    self.NBASIS = args[0]
                    o = self.NBASIS + 1
                    self.BASIS = []
                # Load each atom
                else:
                    # X, Y, Z, Atom Type
                    site = [args[0], args[1], args[2], args[3]]
                    self.BASIS.append(site)
            # Number of atom types
            if n == 25:
                # Load number of atoms
                if o <= 0:
                    self.NTYPES = args[0]
                    self.ATOMS = []
                    self.SPRINGS = []
                    # 2 lines per atom.
                    o = self.NTYPES * 2 + 1
                # Load atom parameters
                else:
                    # Even entries are the atom
                    if o % 2 == 0:
                        # Mass, Charge, and Symbol
                        atom = [args[0], args[1], args[2]]
                        self.ATOMS.append(atom)
                    #Odd entries are the spring
                    else:
                        # Spring parameters, X, Y and Z
                        spring = [args[0], args[1], args[2]]
                        self.SPRINGS.append(spring)
            # CORR,ATOMK,RNEIGH
            if n == 26:
                self.CORR = args[0]
                self.ATOMK = args[1]
                self.RNEIGH = args[2]
            # load in the face
            if n == 27:
                self.face = [args[0], args[1], args[2]]
                if len(args) > 3:
                    self.load_crystal = args[3]
                if self.load_crystal:
                    self.loaded_face = [args[4], args[5], args[6]]
            if n==28:
                self.F_a = args[0]
                self.F_b = args[1]
                
            # Decrement our sub-line first.
            o = o - 1
            # Only increment this if not in a sub-line section.
            if o <= 0:
                n = n+1
        file.close()
        return

    def isMonteCarlo(self):
        return self.NWRITX == 666 and self.NWRITY == 666

    def setMonteCarlo(self, montecarlo):
        if montecarlo:
            self.NWRITX = 666
            self.NWRITY = 666
        else:
            self.NWRITX = 666
            self.NWRITY = 888

    def isGridScat(self):
        return self.NWRITX == 666 and self.NWRITY == 777

    def setGridScat(self, grid):
        if grid:
            self.NWRITX = 666
            self.NWRITY = 777
        else:
            self.NWRITX = 666
            self.NWRITY = 888

    def save(self, file=None):
        if file is None:
            output = open(self.fileIn, 'w')
        else:
            output = open(file, 'w')
            if self.load_crystal:
                # We need to copy the crystal file over as well.
                crys_file_in = self.fileIn.replace('.input', '.crys_in')
                crys_file_out = file.replace('.input', '.crys_in')
                try:
                    shutil.copy(crys_file_in, crys_file_out)
                except:
                    print("Error copying the crystal file over!")
                    print(crys_file_in)
                    print(crys_file_out)
                                
            
        # Ensure phi is between -180 and 180
        while self.PHI0 > 180:
            self.PHI0 -= 360
        while self.PHI0 < -180:
            self.PHI0 += 360

        file = serialize(self.E0, self.THETA0, self.PHI0, self.MASS, self.SYMION) + '\n' \
            +  serialize(self.EMIN, self.EMAX, self.ESIZE, self.ASIZE)  + '\n' \
            +  serialize(self.NDTECT)  + '\n' \
            +  serializeArr(self.DTECTPAR)  + '\n' \
            +  serialize(self.DELLOW, self.DELT0)  + '\n' \
            +  serialize(self.DEMAX, self.DEMIN, self.ABSERR)  + '\n' \
            +  serialize(self.NPART)  + '\n' \
            +  serialize(self.RECOIL)  + '\n' \
            +  serialize(self.Z1)  + '\n' \
            +  serialize(self.MAX_STEPS)  + '\n' \
            +  serialize(self.RRMIN, self.RRSTEP)  + '\n' \
            +  serialize(self.ZMIN, self.ZSTEP)  + '\n' \
            +  serialize(self.MAXDIV, self.MINDIV)  + '\n'
        
        file = file \
            +  serialize(self.NUMCHA)  + '\n' \
            +  serialize(self.XSTART, self.XSTEP, self.XSTOP)  + '\n' \
            +  serialize(self.YSTART, self.YSTEP, self.YSTOP)  + '\n'
            
        file = file \
            +  serialize(self.NWRITX, self.NWRITY)  + '\n' \
            +  serialize(self.RAX, self.RAY)  + '\n' \
            +  serialize(self.NPAR, self.IPOT)  + '\n' \
            +  serializeArr(self.POTPAR)  + '\n' \
            +  serialize(self.NIMPAR, self.IIMPOT)  + '\n' \
            +  serializeArr(self.PIMPAR)  + '\n'
        
        if self.IIMPOT == 2:
            file = file \
            +  serialize(self.NBZ, self.TOL)  + '\n' \
            +  serializeArr(self.ZMAX)  + '\n' \
            +  serializeArr(self.NZ)  + '\n' \
            +  serialize(self.NBG, self.GTOL)  + '\n' \
            +  serializeArr(self.GMAX)  + '\n' \
            +  serializeArr(self.NG)  + '\n'
            
        file = file \
            +  serialize(self.TEMP, self.SEED, self.Ion_Index)  + '\n' \
            +  serialize(self.IMAGE)  + '\n' \
            +  serialize(self.SENRGY, self.BDIST)  + '\n' \
            +  serialize(self.AX, self.AY, self.AZ)  + '\n' \
            +  serialize(self.NBASIS)  + '\n'
        
        for i in range(self.NBASIS):
            file = file + serializeArr(self.BASIS[i])  + '\n'
        
        file = file + serialize(self.NTYPES) + '\n'
        
        for i in range(self.NTYPES):
            atom = self.ATOMS[i]
            spring = self.SPRINGS[i]
            file = file \
            +  serializeArr(atom)  + '\n' \
            +  serializeArr(spring)  + '\n' \
        
        file = file + serialize(self.CORR, self.ATOMK, self.RNEIGH) + '\n'
        
        file = file + serializeArr(self.face) + ' ' \
                    + serialize(self.load_crystal) + ' ' \
                    + serializeArr(self.loaded_face) + '\n'

        if self.F_b != 0 or self.F_a != 0:
            file = file + serialize(self.F_a, self.F_b) + '\n'
            
        output.write(file)
        output.close()
        return
    
