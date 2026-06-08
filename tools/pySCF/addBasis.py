from pyscf import gto
import numpy as np

a = 1e-6

custom_basis = {
    'C': [ [0, [1e-6, 1.0] ] ] + gto.basis.load('sto-3g', 'C') ,
    'O': gto.basis.load('sto-3g', 'O')  # base standard
}

mol = gto.M(
    atom = '''
    C 0 0 -1.21601
    O 0 0 0.912008
    ''',
    basis = custom_basis
)
mol.unit = 'Bohr'
mol.spin = 0

nc=150
print("Base:")
print("="*40)
for ib in range(mol.nbas):
    print("Shell        :", ib)
    print("l            :", mol.bas_angular(ib))
    print("# primitives :", mol.bas_nprim(ib))
    print("Expo         :", mol.bas_exp(ib))
    print("coefficients :", mol.bas_ctr_coeff(ib).reshape(-1))

print('-'*nc)

mol.build()

nc=150
print("Base:")
print("="*40)
for ib in range(mol.nbas):
    print("Shell        :", ib)
    print("l            :", mol.bas_angular(ib))
    print("# primitives :", mol.bas_nprim(ib))
    print("Expo         :", mol.bas_exp(ib))
    print("coefficients :", mol.bas_ctr_coeff(ib).reshape(-1))

print('-'*nc)

f = (2*a / np.pi)**0.75
print(f)

V_ao = mol.intor('int1e_nuc')


print("Nombre de fonctions:", mol.nao_nr())
print("integrales:\n", V_ao/f)
print("integrales:\n", V_ao)
