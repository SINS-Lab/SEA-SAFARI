class XYZ:
    def __init__(self):
        self.xyzs = []
        
    def load(self, filename):
        file = open(filename, 'r')
        line = file.readline()
        var = line
        while var is not None:
            xyz = XYZ_Single()
            num = int(var)
            xyz.load(file, num)
            self.xyzs.append(xyz)
            var = file.readline()
            if (len(var.split()) == 0): 
                break
        file.close()
        
    def save(self, filename):
        file = open(filename, 'w')
        for xyz in self.xyzs:
            xyz.write(file)
        file.close()

class XYZ_Single:
    def __init__(self):
        self.number = 0
        # For safari, this is time since start.
        self.comment = ''
        # This is an array of atomic symbols
        self.atoms = []
        # This is an array of the arrays of values for the atoms above.
        self.values = []
        
    def load(self, file, number):
        self.number = number
        # Reading in like this includes a \n at the end.
        self.comment = file.readline().split()[0]
        for i in range(number):
            line = file.readline().split()
            # this is the array of values for the atom on this line
            value = []
            for j in range(len(line) - 1):
                value.append(float(line[j + 1]))
            self.values.append(value)
            self.atoms.append(line[0])
    
    def write(self, file):
        file.write(str(self.number)+'\n')
        #Comment is read in with the formatting, so includes own /n
        file.write(self.comment)
        for i in range(self.number):
            line = self.atoms[i]
            value = self.values[i]
            for j in range(len(value)):
                line = line + ' ' + str(value[j])
            line = line + '\n'
            file.write(line)
        return