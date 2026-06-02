"""
===========================================================================================
Ground-State to Excited-State density matrix via TDDFT with Tamm-Dancoff Approximation (TDA)
============================================================================================
"""
import numpy as np
from pyscf import gto, scf, tdscf
from pyscf.tools import cubegen

def printMatrix(A):
	for row in A:
		print(" ".join(f"{x:16.10f}" for x in row))


MOL_ATOM = """
  H   0.000000   0.000000  0.0
  H   0.000000   0.000000  1.0
"""
#MOL_BASIS  = "sto-3g"
MOL_BASIS  = "3-21g"
MOL_CHARGE = 0
MOL_SPIN   = 0          # 2S (0 = singlet)

XC_FUNCTIONAL = "blyp"

N_STATES = 3

CUBE_NX = CUBE_NY = CUBE_NZ = 80

mol = gto.Mole()
mol.atom   = MOL_ATOM
mol.basis  = MOL_BASIS
mol.charge = MOL_CHARGE
mol.spin   = MOL_SPIN
mol.unit = "Angstrom"
mol.verbose = 3
mol.build()

mf = scf.RKS(mol)
mf.xc = XC_FUNCTIONAL
mf.conv_tol = 1e-12
mf.conv_tol_grad = 1e-10
mf.run()

mo_coeff = mf.mo_coeff
mo_occ   = mf.mo_occ
S        = mf.get_ovlp()

nao, nmo = mo_coeff.shape
nocc = np.count_nonzero(mo_occ > 0)
nvir = nmo - nocc
occ_idx  = np.where(mo_occ > 0)[0]
virt_idx = np.where(mo_occ == 0)[0]
print("mo coeff") 
printMatrix(mo_coeff)

print(f"\n{'='*62}")
print(f"  Molecule : {mol.atom.strip()[:40]}")
print(f"  Basis    : {mol.basis}   XC : {mf.xc}")
print(f"  nao={nao}  nmo={nmo}  nocc={nocc}  nvir={nvir}")
print(f"  N_elec   = {mol.nelectron}")
print(f"{'='*62}\n")

# TDDFT with TDA 
mytda = tdscf.TDA(mf)
mytda.nstates = N_STATES
mytda.run()

print(f"\n  TDA excitation energies:")
print(f"  {'State':>6}  {'E (eV)':>10}  {'f (osc.)':>10}")
print(f"  {'-'*32}")
for i, (e_ev, osc) in enumerate(
        zip(mytda.e * 27.2114, mytda.oscillator_strength())):
    print(f"  {i+1:>6}  {e_ev:>10.4f}  {osc:>10.6f}")


#  Ground-state 1-RDM (reference)
# In MO basis: diagonal with occupation numbers (0 or 2 for RKS)
dm1_gs_mo = np.diag(mo_occ)                        # (nmo, nmo)
dm1_gs_ao = mf.make_rdm1()                         # AO basis (standard)


def tda_trdm1_mo(X_raw):
    """
    Compute the excited-state 1-RDM in the MO basis from TDA amplitudes.
    The factor 2 comes from summing over alpha and beta spins.
    """
    # Normalise X so that  Tr(X^T X) = 1
    norm = np.sqrt(np.sum(X_raw ** 2))
    #X = X_raw #/ norm                           # shape (nocc, nvir)
    X = X_raw / norm                           # shape (nocc, nvir)
    print(X.shape)
    print('X=')
    printMatrix(X)

    dm1 = np.zeros((nmo, nmo))
    print(dm1.shape)
    print(nocc,nvir)
    dm1[np.ix_(occ_idx, virt_idx)] = X  # (nocc, nvirt) dans le bloc occ-virt

    return dm1


for ist in range(N_STATES):
    X_raw = np.array(mytda.xy[ist][0])  
    if X_raw.ndim == 3:
        X_raw = X_raw[0]               

    dm1_mo = tda_trdm1_mo(X_raw)
    np.save("T_RDM_state_"+str(ist+1),dm1_mo)
    np.save("X_state_"+str(ist+1),X_raw)
    print("TR-1-RDM on MO")
    printMatrix(dm1_mo)

    #  Tr-1-RDM in AO basis
    dm1_ao = mo_coeff @ dm1_mo @ mo_coeff.T
    print("Tr 1-RDM on AO")
    printMatrix(dm1_ao)

    trace = np.trace(S @ dm1_ao)

    print(f"\n  State {ist+1}  "
          f"(ΔE = {mytda.e[ist]*27.2114:.4f} eV) ")
    print(f"     Tr(S·gamma_AO) = {trace:.8f}  "
          f"[ref: {mol.nelectron}]")

    # Excited-state transition density  rho_exc(r)
    cube_exc = f"rho_tr_tda_state{ist+1}.cube"
    cubegen.density(mol, cube_exc, dm1_ao,
                    nx=CUBE_NX, ny=CUBE_NY, nz=CUBE_NZ)

    print(f"     Cube files written:")
    print(f"       {cube_exc}")

