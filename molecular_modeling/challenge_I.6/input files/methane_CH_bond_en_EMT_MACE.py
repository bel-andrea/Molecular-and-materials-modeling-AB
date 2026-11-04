from ase.build import molecule
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from mace.calculators import mace_mp
import numpy as np

# Experimental reference values for methane C-H bond
EXP_CH_BOND_LENGTH = 1.087  # Å
EXP_CH_BOND_ENERGY = 4.55   # eV (approximately 439 kJ/mol)
EXP_CH_BOND_ENERGY_kJ = 439.0  # kJ/mol

print("=" * 70)
print("METHANE C-H BOND ENERGY - EMT vs MACE COMPARISON")
print("=" * 70)

def calculate_methane_ch_bond_emt():
    """Calculate methane C-H bond properties using EMT."""
    print("\n" + "-" * 35)
    print("EMT CALCULATIONS")
    print("-" * 35)
    
    try:
        # 1. Setup Methane (CH4)
        ch4 = molecule('CH4')
        ch4.calc = EMT()
        dyn_c4 = BFGS(ch4, trajectory='ch4_geom_opt_emt.traj', logfile=None)
        dyn_c4.run(fmax=0.01)
        e_methane = ch4.get_potential_energy()
        ch_bond_length_emt = ch4.get_distance(0, 1) # C-H bond length

        # 2. Setup Methyl Radical (CH3)
        ch3 = molecule('CH3')
        ch3.calc = EMT()
        dyn_m = BFGS(ch3, trajectory='ch3_geom_opt_emt.traj', logfile=None)
        dyn_m.run(fmax=0.01)
        e_methyl = ch3.get_potential_energy()

        # 3. Setup Hydrogen Atom (H)
        h_atom = molecule('H')
        h_atom.calc = EMT()
        e_h = h_atom.get_potential_energy()

        # 4. Calculate Bond Dissociation Energy (BDE)
        # BDE = E(CH3) + E(H) - E(CH4)
        bond_energy_emt = (e_methyl + e_h) - e_methane
        
        print(f"Methane Energy:       {e_methane:.4f} eV")
        print(f"Methyl Radical Energy:{e_methyl:.4f} eV")
        print(f"Hydrogen Atom Energy: {e_h:.4f} eV")
        print(f"\nEMT Results:")
        print(f"  C-H Bond Length:  {ch_bond_length_emt:.4f} Å")
        print(f"  C-H Bond Energy:  {bond_energy_emt:.4f} eV")
        print(f"  C-H Bond Energy:  {bond_energy_emt * 96.485:.2f} kJ/mol")
        
        return ch_bond_length_emt, bond_energy_emt, True
    except Exception as e:
        print(f"EMT calculation failed (likely due to missing potential parameters for H/C): {e}")
        return np.nan, np.nan, False

def calculate_methane_ch_bond_mace():
    """Calculate methane C-H bond properties using MACE."""
    print("\n" + "-" * 35)
    print("MACE CALCULATIONS")
    print("-" * 35)
    
    calc_mace = mace_mp(model="medium", device="cpu", default_dtype="float64")
    
    # 1. Setup Methane (CH4)
    ch4 = molecule('CH4')
    ch4.calc = calc_mace
    dyn_c4 = BFGS(ch4, trajectory='ch4_geom_opt_mace.traj', logfile=None)
    dyn_c4.run(fmax=0.01)
    e_methane = ch4.get_potential_energy()
    ch_bond_length_mace = ch4.get_distance(0, 1)
    print(f"Methane Energy: {e_methane:.4f} eV")

    # 2. Setup Methyl Radical (CH3)
    ch3 = molecule('CH3')
    ch3.calc = calc_mace
    dyn_m = BFGS(ch3, trajectory='ch3_geom_opt_mace.traj', logfile=None)
    dyn_m.run(fmax=0.01)
    e_methyl = ch3.get_potential_energy()
    print(f"Methyl Radical Energy: {e_methyl:.4f} eV")

    # 3. Setup Hydrogen Atom (H)
    h_atom = molecule('H')
    h_atom.calc = calc_mace
    e_h = h_atom.get_potential_energy()
    print(f"Hydrogen Atom Energy: {e_h:.4f} eV")

    # 4. Calculate Bond Dissociation Energy (BDE)
    bond_energy_mace = (e_methyl + e_h) - e_methane
    
    print(f"\nMACE Results:")
    print(f"  C-H Bond Length:  {ch_bond_length_mace:.4f} Å")
    print(f"  C-H Bond Energy:  {bond_energy_mace:.4f} eV")
    print(f"  C-H Bond Energy:  {bond_energy_mace * 96.485:.2f} kJ/mol")
    
    return ch_bond_length_mace, bond_energy_mace, True

# ============================================
# Run both calculations
# ============================================
ch_length_emt, bond_energy_emt, emt_success = calculate_methane_ch_bond_emt()
ch_length_mace, bond_energy_mace, mace_success = calculate_methane_ch_bond_mace()

# ============================================
# Comparison with Experiment
# ============================================
print("\n" + "=" * 70)
print("COMPARISON WITH EXPERIMENT")
print("=" * 70)

print("\nC-H Bond Length:")
print(f"  Experimental:  {EXP_CH_BOND_LENGTH:.4f} Å")
if emt_success:
    print(f"  EMT:           {ch_length_emt:.4f} Å  (Δ = {ch_length_emt - EXP_CH_BOND_LENGTH:+.4f} Å)")
else:
    print(f"  EMT:           N/A (Failed)")
print(f"  MACE:          {ch_length_mace:.4f} Å  (Δ = {ch_length_mace - EXP_CH_BOND_LENGTH:+.4f} Å)")

print("\nC-H Bond Energy:")
print(f"  Experimental:  {EXP_CH_BOND_ENERGY:.2f} eV ({EXP_CH_BOND_ENERGY_kJ:.1f} kJ/mol)")
if emt_success:
    print(f"  EMT:           {bond_energy_emt:.4f} eV ({bond_energy_emt * 96.485:.2f} kJ/mol)")
else:
    print(f"  EMT:           N/A (Failed)")
print(f"  MACE:          {bond_energy_mace:.4f} eV ({bond_energy_mace * 96.485:.2f} kJ/mol)")

# ============================================
# Summary Table
# ============================================
print("\n" + "=" * 70)
print("SUMMARY TABLE")
print("=" * 70)

print(f"{'Method':<12} {'C-H Length (Å)':<18} {'Bond Energy (eV)':<18} {'Bond Energy (kJ/mol)':<18}")
print("-" * 66)
if emt_success:
    print(f"{'EMT':<12} {ch_length_emt:<18.4f} {bond_energy_emt:<18.4f} {bond_energy_emt * 96.485:<18.2f}")
else:
    print(f"{'EMT':<12} {'N/A':<18} {'N/A':<18} {'N/A':<18}")
print(f"{'MACE':<12} {ch_length_mace:<18.4f} {bond_energy_mace:<18.4f} {bond_energy_mace * 96.485:<18.2f}")
print(f"{'Exp.':<12} {EXP_CH_BOND_LENGTH:<18.4f} {EXP_CH_BOND_ENERGY:<18.2f} {EXP_CH_BOND_ENERGY_kJ:<18.1f}")
print("-" * 66)

print("\nOptimization complete! Trajectory files saved.")
