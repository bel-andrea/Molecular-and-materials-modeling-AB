from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from mace.calculators import mace_mp
import numpy as np
import warnings

warnings.filterwarnings('ignore')

# Dictionary of diatomic molecules with experimental reference values
diatomics = {
    "H2": {
        "symbol": "H",
        "exp_bond_length": 0.741,
        "exp_atomization_energy": 4.48,
        "exp_frequency": 4401.2,
        "exp_zpe": 0.273,
        "init_bond": 0.8,
    },
    "N2": {
        "symbol": "N",
        "exp_bond_length": 1.098,
        "exp_atomization_energy": 9.76,
        "exp_frequency": 2358.57,
        "exp_zpe": 0.145,
        "init_bond": 1.1,
    },
    "O2": {
        "symbol": "O",
        "exp_bond_length": 1.208,
        "exp_atomization_energy": 5.16,
        "exp_frequency": 1580.0,
        "exp_zpe": 0.098,
        "init_bond": 1.2,
    },
}

print("=" * 80)
print("DIATOMIC MOLECULES THERMOCHEMISTRY: EMT vs MACE COMPARISON")
print("=" * 80)

# Initialize MACE calculator once outside the loop for efficiency
print("\nInitializing MACE calculator...")
calc_mace = mace_mp(model="medium", device="cpu", default_dtype="float64")

all_results = []

for mol_name, data in diatomics.items():
    symbol = data["symbol"]
    d_init = data["init_bond"]
    exp_bond = data["exp_bond_length"]
    exp_atom_en = data["exp_atomization_energy"]
    exp_freq = data["exp_frequency"]
    exp_zpe = data["exp_zpe"]

    print("\n" + "=" * 80)
    print(f"PROCESSING MOLECULE: {mol_name}")
    print("=" * 80)

    mol_results = {"molecule": mol_name}

    # ============================================
    # 1. EMT Calculations
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

        # Vibrational analysis
        vib_emt = Vibrations(molecule_emt, name=f"{mol_name}_vib_emt")
        vib_emt.run()
        vib_energies_emt = vib_emt.get_energies()
        vib_freqs_emt_cm1 = vib_emt.get_frequencies()
        
        positive_modes_emt = [np.real(f) for f in vib_freqs_emt_cm1 if abs(np.real(f)) > 1.0]
        stretching_freq_emt = max(positive_modes_emt) if positive_modes_emt else 0.0

        # Thermochemistry (ZPE & IdealGasThermo)
        vib_energies_filtered_emt = np.array([np.real(e) for e in vib_energies_emt if np.real(e) > 0.001])
        zpe_emt = 0.5 * np.sum(vib_energies_filtered_emt)

        thermo_emt = IdealGasThermo(
            vib_energies=vib_energies_filtered_emt,
            potentialenergy=e_molecule_emt,
            atoms=molecule_emt,
            geometry='linear',
            symmetrynumber=2,
            spin=0
        )
        temp = 298.15
        press = 101325.0
        H_emt = thermo_emt.get_enthalpy(temperature=temp, verbose=False)
        S_emt = thermo_emt.get_entropy(temperature=temp, pressure=press, verbose=False)
        G_emt = thermo_emt.get_gibbs_energy(temperature=temp, pressure=press, verbose=False)

        print(f"  Optimized bond length: {d_emt:.4f} Å")
        print(f"  Atomization energy:    {e_atomization_emt:.4f} eV")
        print(f"  Stretching frequency:  {stretching_freq_emt:.2f} cm⁻¹")
        print(f"  ZPE:                   {zpe_emt:.4f} eV")
        print(f"  Enthalpy (H):          {H_emt:.4f} eV")
        print(f"  Entropy (S):           {S_emt:.4f} eV/K")

        mol_results["emt"] = {
            "bond": d_emt, "atom_en": e_atomization_emt, "freq": stretching_freq_emt,
            "zpe": zpe_emt, "H": H_emt, "S": S_emt, "G": G_emt, "success": True
        }
    except Exception as e:
        print(f"  EMT calculation failed: {e}")
        mol_results["emt"] = {"success": False}

    # ============================================
    # 2. MACE Calculations
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

        # Vibrational analysis
        vib_mace = Vibrations(molecule_mace, name=f"{mol_name}_vib_mace")
        vib_mace.run()
        vib_energies_mace = vib_mace.get_energies()
        vib_freqs_mace_cm1 = vib_mace.get_frequencies()
        
        positive_modes_mace = [np.real(f) for f in vib_freqs_mace_cm1 if abs(np.real(f)) > 1.0]
        stretching_freq_mace = max(positive_modes_mace) if positive_modes_mace else 0.0

        # Thermochemistry (ZPE & IdealGasThermo)
        vib_energies_filtered_mace = np.array([np.real(e) for e in vib_energies_mace if np.real(e) > 0.001])
        zpe_mace = 0.5 * np.sum(vib_energies_filtered_mace)

        thermo_mace = IdealGasThermo(
            vib_energies=vib_energies_filtered_mace,
            potentialenergy=e_molecule_mace,
            atoms=molecule_mace,
            geometry='linear',
            symmetrynumber=2,
            spin=0
        )
        H_mace = thermo_mace.get_enthalpy(temperature=temp, verbose=False)
        S_mace = thermo_mace.get_entropy(temperature=temp, pressure=press, verbose=False)
        G_mace = thermo_mace.get_gibbs_energy(temperature=temp, pressure=press, verbose=False)

        print(f"  Optimized bond length: {d_mace:.4f} Å")
        print(f"  Atomization energy:    {e_atomization_mace:.4f} eV")
        print(f"  Stretching frequency:  {stretching_freq_mace:.2f} cm⁻¹")
        print(f"  ZPE:                   {zpe_mace:.4f} eV")
        print(f"  Enthalpy (H):          {H_mace:.4f} eV")
        print(f"  Entropy (S):           {S_mace:.4f} eV/K")

        mol_results["mace"] = {
            "bond": d_mace, "atom_en": e_atomization_mace, "freq": stretching_freq_mace,
            "zpe": zpe_mace, "H": H_mace, "S": S_mace, "G": G_mace, "success": True
        }
    except Exception as e:
        print(f"  MACE calculation failed: {e}")
        mol_results["mace"] = {"success": False}

    mol_results["exp"] = {
        "bond": exp_bond, "atom_en": exp_atom_en, "freq": exp_freq, "zpe": exp_zpe
    }
    all_results.append(mol_results)

# ============================================
# 3. Comprehensive Summary Tables
# ============================================
print("\n" + "=" * 90)
print("FINAL COMPREHENSIVE SUMMARY ACROSS DIATOMICS")
print("=" * 90)

for res in all_results:
    mol = res["molecule"]
    print(f"\n>>> Molecule: {mol} <<<")
    print(f"{'Method':<8} {'Bond (Å)':<12} {'Atomization (eV)':<18} {'Freq (cm⁻¹)':<16} {'ZPE (eV)':<12} {'Entropy (S)':<14}")
    print("-" * 80)
    
    exp = res["exp"]
    print(f"{'Exp.':<8} {exp['bond']:<12.4f} {exp['atom_en']:<18.2f} {exp['freq']:<16.2f} {exp['zpe']:<12.4f} {'N/A':<14}")
    
    if res.get("emt", {}).get("success"):
        e = res["emt"]
        print(f"{'EMT':<8} {e['bond']:<12.4f} {e['atom_en']:<18.4f} {e['freq']:<16.2f} {e['zpe']:<12.4f} {e['S']:<14.4f}")
    else:
        print(f"{'EMT':<8} {'Failed / Not Available':<45}")
        
    if res.get("mace", {}).get("success"):
        m = res["mace"]
        print(f"{'MACE':<8} {m['bond']:<12.4f} {m['atom_en']:<18.4f} {m['freq']:<16.2f} {m['zpe']:<12.4f} {m['S']:<14.4f}")
    else:
        print(f"{'MACE':<8} {'Failed':<45}")

print("\n" + "=" * 80)
print("All multi-molecule thermochemistry calculations complete!")
