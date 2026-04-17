import numpy as np
import pandas as pd
from pyscf import gto, scf
from pyscf.tools import molden

def printMatrix(A):
	for row in A:
		print(" ".join(f"{x:16.10f}" for x in row))

mol = gto.Mole()
mol.atom = 'C -3.31099285 0.269371293 0.00020037396 ; C -1.06359641 -0.847358644 -0.000470478814 ; H -3.47318988 2.3145633 5.67547746e-05 ; H -5.05700156 -0.797658734 0.000578319185 ; H -0.84309566 -2.88508435 0.000330765063 ; C 1.27650819 0.658637737 -0.000134107564 ; H 0.992251865 2.73610249 -0.000188909622 ; O 3.3711902, -0.231478129, 0.000206043138'
mol.basis = '631+g*'
#mol.basis = 'cc-pvtz'
#mol.basis = 'sto-3g'
mol.unit = 'Bohr'
mol.spin = 0
mol.build()

nc=150

mf = scf.RHF(mol)
mf.kernel()
print('-'*nc)

moldenfile='acrolein_HF_pySCF.molden'
print("Voir le fichier ",moldenfile)
molden.from_scf(mf, moldenfile)

V_ao = mol.intor('int1e_nuc')

np.set_printoptions(precision=10, suppress=True)
print("Intégrales sur la base atomique :")
print("="*40)
#print(V_ao)
printMatrix(V_ao)
df = pd.DataFrame(V_ao)
df.to_csv("V_ao.csv", index=False, header=False)
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

V_mo = C.T @ V_ao @ C

print("Matrice noyau-électron en base MO :")
print("="*40)
#print(V_mo)
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

print("Coef des MO pySCF:")
print("="*40)
printMatrix(C)
print('-'*nc)

"""
import pandas as pd
df = pd.read_csv('CoefGaussian.csv',  header=None, delimiter=r"\s+")
#print(df)
CGaussian=df.to_numpy().reshape(nmo,-1).T
print("Coef des MO Gaussian:")
print("="*40)
printMatrix(CGaussian)
print('-'*nc)

print("Diff Coef des MO PySCF/Gaussian:")
print("="*40)
printMatrix(C-CGaussian)
print('-'*nc)

V_moGaussian = CGaussian.T @ V_ao @ CGaussian
nocc = mol.nelectron // 2
D_mo = np.zeros((nmo, nmo))
D_mo[:nocc, :nocc] = 2 * np.eye(nocc)
E_ne_mogauss = np.einsum('pq,qp', D_mo, V_moGaussian)
print(f"nOcc               = {nocc:12d}")
print(f"E_ne (MO Gaussian) = {E_ne_mogauss:.12f}")
print(f"E_ne (MO pySCF)    = {E_ne_mo:.12f}")
print(f"Diff(meV)     )    = {(E_ne_mogauss-E_ne_mo)*27.21*1000:.12f}")
print('-'*nc)


df = pd.read_csv('Vao.csv',  header=None, delimiter=r"\s+")
#print(df)
V_ao_cdftt=df.to_numpy().T
print("Vao cdftt:")
print("="*40)
printMatrix(V_ao_cdftt)
print('-'*nc)
print("Vao pySCF:")
print("="*40)
printMatrix(V_ao)
print('-'*nc)
print("Diff vao pySCF/cdftt:")
print("="*40)
printMatrix(V_ao-V_ao_cdftt)
print('-'*nc)

V_moGaussian_cdftt = CGaussian.T @ V_ao_cdftt @ CGaussian
nocc = mol.nelectron // 2
D_mo = np.zeros((nmo, nmo))
D_mo[:nocc, :nocc] = 2 * np.eye(nocc)
E_ne_mogauss_cdftt = np.einsum('pq,qp', D_mo, V_moGaussian_cdftt)
print(f"nOcc               = {nocc:12d}")
print(f"E_ne (MO Gaussian/Vao cdftt) = {E_ne_mogauss_cdftt:.12f}")
print('-'*nc)
"""
