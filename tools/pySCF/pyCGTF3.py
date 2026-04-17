from pyscf import gto, dft
import numpy as np

def printMatrix(A):
	for row in A:
		print(" ".join(f"{x:16.10f}" for x in row))


mol = gto.Mole()
mol.atom = 'C 0 0 -1.21601177; O 0 0 0.912008824'
mol.basis = 'sto-3g'
mol.unit = 'Bohr'
mol.spin = 0
mol.build()

grid = dft.gen_grid.Grids(mol)
grid.level = 4
grid.build()

ao = dft.numint.eval_ao(mol, grid.coords)

nao = ao.shape[1]

T = np.zeros((nao, nao, nao))

weights = grid.weights
for g in range(len(weights)):
    w = weights[g]
    v = ao[g]
    T += w * np.einsum("i,j,k->ijk", v, v, v)

for i in range(T.shape[0]):
	print('i=',i)
	printMatrix(T[i])

