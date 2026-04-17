import numpy as np
from pyscf import gto, scf

def printMatrix(A):
	for row in A:
		print(" ".join(f"{x:12.6f}" for x in row))

mol = gto.Mole()
mol.atom = 'C 0 0 -1.21601; O 0 0 0.912008'
#mol.basis = 'cc-pvtz'
mol.basis = 'sto-3g'
mol.unit = 'Bohr'
mol.spin = 0
mol.build()

nc=150

mf = scf.RHF(mol)
mf.kernel()
print('-'*nc)

V_ao = mol.intor('int1e_nuc')

np.set_printoptions(precision=10, suppress=True)
print("Intégrales sur la base atomiques:")
print("="*40)
#print(V_ao)
printMatrix(V_ao)
print(f"Somme des éléments: {np.sum(V_ao):.12f}")
print('-'*nc)

print("Base:")
print("="*40)
for ib in range(mol.nbas):
    print("Shell        :", ib)
    print("l            :", mol.bas_angular(ib))
    print("# primitives :", mol.bas_nprim(ib))
    print("Expo         :", mol.bas_exp(ib))
    print("coefficients :", mol.bas_ctr_coeff(ib).reshape(-1))

print('-'*nc)


C = mf.mo_coeff

print("Coefficients de la base MO:")
print("="*40)
printMatrix(C)
print('-'*nc)

V_mo = C.T @ V_ao @ C

print("Matrice noyau-électron en base MO :")
print("="*40)
printMatrix(V_mo)
print(f"Somme des éléments: {np.sum(V_mo):.12f}")
print('-'*nc)

# calcul de E_ne par matrice densité
D = mf.make_rdm1()
E_ne  = np.einsum('ij,ji', D, V_ao)
print(f"E_ne  = {E_ne:.12f}")

nmo = C.shape[1]
nocc = mol.nelectron // 2
D_mo = np.zeros((nmo, nmo))
D_mo[:nocc, :nocc] = 2 * np.eye(nocc)
E_ne_mo = np.einsum('pq,qp', D_mo, V_mo)
print(f"E_ne (MO) = {E_ne_mo:.12f}")
print('-'*nc)
