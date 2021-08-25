class PotProvider:

    def __init__(self, r_max, r_step):
        self.r_max = r_max
        self.r_min = r_step
        self.r_step = r_step
        self.abr = "RSQ"

    def V_r(self, A, B, r):
        return 1/r

    def dV_dr(self, A, B, r):
        return 1/r**2

    def name(self):
        return self.abr

    def abbrv(self):
        return self.abr

    def print_pots(self, name, A, B, r):
        V_r = self.V_r(A, B, r)
        dV_dr = self.dV_dr(A, B, r)
        fmt = '{}\t{}\t{:.5e}\t{:.5e}\n'
        pots_file = open(name+'.pots', 'w')
        header = '# Generated {} Potentials for {}-{}.\n#\n'+\
                '# Table Range in file (min, step, max): {}-{}-{}\n#\n'+\
                '# Atom1\tAtom2\tV_r\t-dV_dr\n'
        header = header.format(self.name(), A, B, self.r_min, self.r_step, self.r_max)
        pots_file.write(header)
        for i in range(len(V_r)):
            pots_file.write(fmt.format(A, B, V_r[i], dV_dr[i]))
        pots_file.close()