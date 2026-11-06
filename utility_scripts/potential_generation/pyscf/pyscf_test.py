from pyscf import gto, scf, dft
import numpy as np
import os

def run_scf_mol(mol):
    hf = scf.UHF(mol)
    hf.max_cycle = 5000
    hf.kernel()
    hf.analyze()
    ump2 = hf.MP2().run() # Compute second order corrections (https://pyscf.org/user/mp.html)
    return ump2.e_hf + ump2.e_corr, hf.converged # Return first order + corrections

def run_dft_mol(mol):
    hf = dft.UKS(mol)
    hf.max_cycle = 500
    hf.kernel()
    hf.analyze()
    return hf.e_tot, hf.converged

def test_scf(mol_src, dir):
    if not "./generated/" in dir:
        dir = f"./generated/{dir}"
    if not os.path.exists(dir):
        os.makedirs(dir)

    def make_mol(r):
        mol = mol_src(r)
        mol.verbose = 4
        mol.output = f"{dir}/log.txt"
        mol.build()
        return mol

    filename = f"{dir}/far.log"
    with open(filename,'w') as file:
        pass
    e_p = 0
    for r in np.arange(1.0, 26, 1):
        print(r)
        mf, converged = run_scf_mol(make_mol(r))
        if not converged:
            break
        E = mf
        with open(filename,'a') as file:
            file.write(f"{r}\t{E}\t{1.0 if converged else 0.0:.0f}\n")
        if converged and abs(E-e_p) < 0.00001:
            break
        e_p = E

    filename = f"{dir}/near.log"
    with open(filename,'w') as file:
        pass
    for r in np.arange(3.0, 0.005, -0.005):
        print(r)
        mf, converged = run_scf_mol(make_mol(r))
        E = mf
        with open(filename,'a') as file:
            file.write(f"{r}\t{E}\t{1.0 if converged else 0.0:.1f}\n")

def test_dft(mol_src, dir):
    if not "./generated/" in dir:
        dir = f"./generated/{dir}"
    if not os.path.exists(dir):
        os.makedirs(dir)

    def make_mol(r):
        mol = mol_src(r)
        mol.verbose = 4
        mol.output = f"{dir}/log.txt"
        mol.build()
        return mol

    file = open(f"{dir}/far.log", 'w')
    e_p = 0
    for r in np.arange(1.0, 26, 1):
        print(r)
        mf, converged = run_dft_mol(make_mol(r))
        E = mf
        file.write(f"{r}\t{E}\t{1.0 if converged else 0.0:.0f}\n")
        if converged and abs(E-e_p) < 0.00001:
            break
        e_p = E

    file = open(f"{dir}/near.log", 'w')
    for r in np.arange(3.0, 0.05, -0.025):
        print(r)
        mf, converged = run_dft_mol(make_mol(r))
        E = mf
        file.write(f"{r}\t{E}\t{1.0 if converged else 0.0:.0f}\n")

def Au_Li(r):
    return gto.M(
        atom=f"Au 0 0 0; Li 0 0 {r}", 
        basis={"Li":"sto3-g", "Au":"lanl2dz"},
        ecp = {'Au': 'crenbl'},
        spin=2
    )
def Xe_C(r):
    return gto.M(
        atom=f"Xe 0 0 0; C 0 0 {r}", 
        basis={"C":"sto3-g", "Xe":"lanl2dz"},
        ecp = {'Xe': 'lanl2dz'},
        spin=2
    )
def Cu_Ar(r):
    basis_Ar = "lanl2dz"#gto.basis.load("./lanl2dz ecp.0.nw", 'Ar')
    basis_Cu = "lanl2dz"#gto.basis.load("./lanl2dz ecp.0.nw", 'Cu')
    basis_Ar_ecp = gto.basis.load_ecp("./lanl2dz ecp.0.nw", 'Ar')
    basis_Cu_ecp = gto.basis.load_ecp("./lanl2dz ecp.0.nw", 'Cu')
    return gto.M(
        atom=f"Cu 0 0 0; Ar 0 0 {r}", 
        basis={"Ar":basis_Ar, "Cu":basis_Cu},
        ecp = {'Ar':basis_Ar_ecp,'Cu':basis_Cu_ecp},
        spin=1,
    )
def Cu_Na(r):
    basis_Na = "crenbl"
    basis_Cu = "crenbl"
    # basis_Na_ecp ='crenbl'# gto.basis.load_ecp("./lanl2dz ecp.0.nw", 'Na')
    # basis_Cu_ecp ='crenbl'# gto.basis.load_ecp("./lanl2dz ecp.0.nw", 'Cu')
    mol = gto.Mole()
    mol.atom = f"Cu 0 0 0; Na 0 0 {r}"
    mol.basis = {'Na':basis_Na,'Cu':basis_Cu}
    # mol.ecp = {'Cu': basis_Cu_ecp, 'Na': basis_Na_ecp}
    mol.spin = 0
    return mol
def Cu_Ga(r):
    basis_Ga = "lanl2dz"
    basis_Cu = "lanl2dz"
    basis_Ga_ecp = basis_Ga#gto.basis.load_ecp("./lanl2dz ecp.0.nw", 'Ga')
    basis_Cu_ecp = basis_Cu#gto.basis.load_ecp("./lanl2dz ecp.0.nw", 'Cu')
    mol = gto.Mole()
    mol.atom = f"Cu 0 0 0; Ga 0 0 {r}"
    mol.basis = {'Ga':basis_Ga,'Cu':basis_Cu}
    mol.ecp = {'Ga':basis_Ga_ecp,'Cu':basis_Cu_ecp}
    mol.spin = 1
    mol.charge = 1
    return mol
def Cu_K(r):
    basis_Cu = "crenbl"
    basis_K = "crenbl"
    basis_K_ecp = basis_K#gto.basis.load_ecp("./lanl2dz ecp.0.nw", 'Ga')
    basis_Cu_ecp = basis_Cu#gto.basis.load_ecp("./lanl2dz ecp.0.nw", 'Cu')
    mol = gto.Mole()
    mol.atom = f"Cu 0 0 0; K 0 0 {r}"
    mol.basis = {'K':basis_K,'Cu':basis_Cu}
    mol.ecp = {'K':basis_K_ecp,'Cu':basis_Cu_ecp}
    mol.spin = 0
    mol.charge = 0
    return mol
def Cu_K_2(r):
    basis_K = "sto3-g"
    basis_Cu = "sto3-g"
    basis_K = "lanl2dz"
    basis_Cu = "lanl2dz"
    basis_K_ecp = "lanl2dz"
    basis_Cu_ecp = "lanl2dz"
    mol = gto.Mole()
    mol.atom = f"Cu 0 0 0; K 0 0 {r}"
    mol.basis = {'K':basis_K,'Cu':basis_Cu}
    mol.ecp = {'K':basis_K_ecp,'Cu':basis_Cu_ecp}
    mol.spin = 0
    return mol

if __name__ == "__main__":
    # test_scf(Xe_C, "Xe_C_1")
    test_scf(Cu_Na, "Cu_Na")
    # test_scf(Cu_Ga, "Cu_Ga_1")
    # test_dft(Cu_Ga, "Cu_Ga_2")
    # test_scf(Cu_Ar, "Cu_Ar_2")
    # test_dft(Cu_Ar, "Cu_Ar_3")
    # test_scf(Cu_K_2, "Cu_K_6")
    # test_scf(Cu_K_2, "Cu_K_8")
    # test_dft(Cu_K, "Cu_K_4")
    print("Done")