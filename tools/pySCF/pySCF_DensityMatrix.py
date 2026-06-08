"""
Excited-State 1-RDM via TDDFT with Tamm-Dancoff Approximation (TDA)
====================================================================
"""
import numpy as np
from pyscf import gto, scf, tdscf
from pyscf.tools import cubegen

MOL_ATOM = """
  C   0.000000   0.000000  -0.655460
  O   0.000000   0.000000   0.490920
"""
#MOL_BASIS  = "6-31g*"
MOL_BASIS  = "sto-3g"
MOL_CHARGE = 0
MOL_SPIN   = 0          # 2S (0 = singlet)

XC_FUNCTIONAL = "blyp"

N_STATES = 1

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
mf.run()

mo_coeff = mf.mo_coeff      
mo_occ   = mf.mo_occ       
S        = mf.get_ovlp()  

nao, nmo = mo_coeff.shape
nocc = np.count_nonzero(mo_occ > 0)
nvir = nmo - nocc

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


# ──────────────────────────────────────────────────────────────
# 4.  Core function: TDA excited-state 1-RDM
# ──────────────────────────────────────────────────────────────
def tda_rdm1_mo(X_raw):
    """
    Compute the excited-state 1-RDM in the MO basis from TDA amplitudes.
    The factor 2 comes from summing over alpha and beta spins.
    """
    # Normalise X so that  Tr(X^T X) = 1
    norm = np.sqrt(np.sum(X_raw ** 2))
    X = X_raw / norm                           # shape (nocc, nvir)

    dm1 = np.zeros((nmo, nmo))

    # Occupied-occupied block: ground-state occ - hole depletion
    dm1[:nocc, :nocc] = 2.0 * np.eye(nocc) - 2.0 * (X @ X.T)

    # Virtual-virtual block: particle gain
    dm1[nocc:, nocc:] = 2.0 * (X.T @ X)

    # Occupied-virtual / virtual-occupied blocks are zero in TDA
    # (no coherences because Y = 0)

    return dm1


# Loop over excited states
print(f"\n{'='*62}")
print("  Excited-state 1-RDM analysis  (TDA)")
print(f"{'='*62}")

for ist in range(N_STATES):

    X_raw = np.array(mytda.xy[ist][0])  
    if X_raw.ndim == 3:
        X_raw = X_raw[0]               

    # 1-RDM in MO basis
    dm1_mo = tda_rdm1_mo(X_raw)
    print("1-RDM on MO")
    print(dm1_mo)

    #  1-RDM in AO basis
    dm1_ao = mo_coeff @ dm1_mo @ mo_coeff.T
    print("1-RDM on AO")
    print(dm1_ao)

    trace = np.trace(S @ dm1_ao)

    print(f"\n  ── State {ist+1}  "
          f"(ΔE = {mytda.e[ist]*27.2114:.4f} eV) ──")
    print(f"     Tr(S·gamma_AO) = {trace:.8f}  "
          f"[ref: {mol.nelectron}]")

    # Excited-state density  rho_exc(r)
    cube_exc = f"rho_tda_state{ist+1}.cube"
    cubegen.density(mol, cube_exc, dm1_ao,
                    nx=CUBE_NX, ny=CUBE_NY, nz=CUBE_NZ)

    print(f"     Cube files written:")
    print(f"       {cube_exc}")

# ──────────────────────────────────────────────────────────────
# 6.  Summary
# ──────────────────────────────────────────────────────────────
print(f"""
{'='*62}
  OUTPUT FILES  (one set per excited state)
{'='*62}
  dm1_mo_tda_state{{I}}.txt          1-RDM in MO basis (text)
  rho_tda_state{{I}}.cube            Excited-state density ρ_exc
  delta_rho_tda_state{{I}}.cube      Difference density Δρ
  rho_attachment_tda_state{{I}}.cube Attachment density (electron gain)
  rho_detachment_tda_state{{I}}.cube Detachment density (electron loss)
{'='*62}
  Visualise .cube files with VESTA, VMD, or Avogadro.
  In VESTA: File > Open cube, then Surface > New surface,
            set ±isovalue (~0.001 a.u.) and colour +/- lobes.
""")
