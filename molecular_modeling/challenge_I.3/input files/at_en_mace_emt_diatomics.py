from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from mace.calculators import mace_mp
import numpy as np

# Dictionary of diatomic molecules with experimental reference values
# (Bond lengths in Å, Atomization energies in eV)
diatomics = {
    "H2": {
        "symbol": "H",
        "exp_bond_length": 0.741,
        "exp_atomization_energy": 4.48,
        "init_bond": 0.8,
    },
    "N2": {
        "symbol": "N",
        "exp_bond_length": 1.098,
        "exp_atomization_energy": 9.76,
        "init_bond": 1.1,
    },
    "O2": {
        "symbol": "O",
        "exp_bond_length": 1.208,
        "exp_atomization_energy": 5.16,
        "init_bond": 1.2,
    },
    "F2": {
        "symbol": "F",
        "exp_bond_length": 1.412,
        "exp_atomization_energy": 1.66,
        "init_bond": 1.4,
    },
}

print("=" * 80)
print("DIATOMIC MOLECULES ATOMIZATION ENERGY: EMT vs MACE COMPARISON")
print("=" * 80)

# Initialize MACE calculator once outside the loop for efficiency
print("\nInitializing MACE calculator...")
calc_mace = mace_mp(model="medium", device="cpu", default_dtype="float64")

results = []

for mol_name, data in diatomics.items():
    symbol = data["symbol"]
    d_init = data["init_bond"]
    exp_bond = data["exp_bond_length"]
    exp_atom_en = data["exp_atomization_energy"]

    print("\n" + "=" * 80)
    print(f"PROCESSING MOLECULE: {mol_name}")
    print("=" * 80)

    # ============================================
    # 1. EMT Calculation
    # ============================================
    print(f"\n--- EMT Calculation for {mol_name} ---")
    try:
        # Single atom energy
        atom_emt = Atoms(symbol, calculator=EMT())
        e_atom_emt = atom_emt.get_potential_energy()

        # Molecule optimization
        molecule_emt = Atoms(f'2{symbol}', [(0., 0., 0.), (0., 0., d_init)], calculator=EMT())
        opt_emt = BFGS(molecule_emt, trajectory=f"{mol_name}_opt_emt.traj", logfile=None)
        opt_emt.run(fmax=0.01)

        d_emt = molecule_emt.get_distance(0, 1)
        e_molecule_emt = molecule_emt.get_potential_energy()
        e_atomization_emt = 2 * e_atom_emt - e_molecule_emt
        n_steps_emt = opt_emt.nsteps
        
        print(f"  Optimized bond length: {d_emt:.4f} Å")
        print(f"  Atomization energy:    {e_atomization_emt:.4f} eV")
    except Exception as e:
        print(f"  EMT calculation failed: {e}")
        d_emt, e_atomization_emt, n_steps_emt = np.nan, np.nan, 0

    # ============================================
    # 2. MACE Calculation
    # ============================================
    print(f"\n--- MACE Calculation for {mol_name} ---")
    try:
        # Single atom energy
        atom_mace = Atoms(symbol, calculator=calc_mace)
        e_atom_mace = atom_mace.get_potential_energy()

        # Molecule optimization
        molecule_mace = Atoms(f'2{symbol}', [(0., 0., 0.), (0., 0., d_init)], calculator=calc_mace)
        opt_mace = BFGS(molecule_mace, trajectory=f"{mol_name}_opt_mace.traj", logfile=None)
        opt_mace.run(fmax=0.01)

        d_mace = molecule_mace.get_distance(0, 1)
        e_molecule_mace = molecule_mace.get_potential_energy()
        e_atomization_mace = 2 * e_atom_mace - e_molecule_mace
        n_steps_mace = opt_mace.nsteps
        
        print(f"  Optimized bond length: {d_mace:.4f} Å")
        print(f"  Atomization energy:    {e_atomization_mace:.4f} eV")
    except Exception as e:
        print(f"  MACE calculation failed: {e}")
        d_mace, e_atomization_mace, n_steps_mace = np.nan, np.nan, 0

    # Store results for summary
    results.append({
        "molecule": mol_name,
        "exp_bond": exp_bond,
        "exp_atom_en": exp_atom_en,
        "d_emt": d_emt,
        "e_atom_emt": e_atomization_emt,
        "steps_emt": n_steps_emt,
        "d_mace": d_mace,
        "e_atom_mace": e_atomization_mace,
        "steps_mace": n_steps_mace,
    })

# ============================================
# 3. Comprehensive Summary Table
# ============================================
print("\n" + "=" * 80)
print("FINAL COMPREHENSIVE SUMMARY TABLE")
print("=" * 80)
print(f"{'Molecule':<10} {'Method':<8} {'Bond Length (Å)':<18} {'Atomization (eV)':<20} {'Steps':<8}")
print("-" * 70)

for res in results:
    mol = res["molecule"]
    print(f"{mol:<10} {'Exp':<8} {res['exp_bond']:<18.4f} {res['exp_atom_en']:<20.2f} {'N/A':<8}")
    print(f"{'':<10} {'EMT':<8} {res['d_emt']:<18.4f} {res['e_atom_emt']:<20.4f} {res['steps_emt']:<8}")
    print(f"{'':<10} {'MACE':<8} {res['d_mace']:<18.4f} {res['e_atom_mace']:<20.4f} {res['steps_mace']:<8}")
    print("-" * 70)

print("\nAll calculations complete! Trajectory files saved.")
