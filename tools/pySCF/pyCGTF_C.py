from pyscf import gto
import numpy as np

mol = gto.Mole()
mol.atom = 'C 0 0 0'
#mol.atom = 'C 0 0 -1.21601; O 0 0 0.912008'
mol.basis = 'sto-3g'
mol.unit = 'Bohr'  
mol.spin = 0
mol.build()

V_nuc = mol.intor('int1e_nuc')

np.set_printoptions(precision=10, suppress=True)
print("Nuclear integrals:")
print(V_nuc)

V_nuc_tot = 0.0
for line in V_nuc:
    for value in line:
        V_nuc_tot += value

print(f"V_nuc_tot = {V_nuc_tot}")

for ib in range(mol.nbas):
    print("Shell:", ib)
    print("l:", mol.bas_angular(ib))
    print("# primitives:", mol.bas_nprim(ib))
    print("Expo:", mol.bas_exp(ib))
    print("coefficients:", mol.bas_ctr_coeff(ib))

