# ASE driven Quantum Espresso geometry optimization of methane
#
import os
from ase import Atoms
from ase.build import molecule
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.optimize import BFGS
from ase.units import Ry, Bohr

os.environ['OMP_NUM_THREADS'] = '1'

# Set the path to your pseudopotential directory
PSEUDO_DIR = '/Users/andreabelendez/Desktop/molmatmodel/mol_mat_model_working/molecular-modeling/ase_mace/ase_tests/04_quantum_espresso_water/pseudo'

# 1. Create the initial methane molecule structure
ch4 = molecule('CH4')
ch4.set_cell([12, 12, 12])  # Create a vacuum box
ch4.center()               # Center the molecule in the box[cite: 9]

# 2. Configure the Quantum ESPRESSO Calculator
pseudopotentials = {
    'C': 'C.upf',
    'H': 'H.upf'
}

input_data = {
    'control': {
        'calculation': 'scf', 
        'prefix': 'ch4',
        'outdir': './outdir',
        'verbosity': 'low',
        'tstress': True,
        'tprnfor': True
    },
    'system': {
        'ecutwfc': 40.0,    # Plane wave cutoff (Ry)[cite: 9]
        'ibrav': 0,         # Use 0 for ASE to handle cell parameters[cite: 9]
        'nosym': True,      # No symmetry for molecules[cite: 9]
        'noinv': True,
        'occupations': 'smearing',  # Better convergence for molecules[cite: 9]
        'smearing': 'gaussian',
        'degauss': 0.02,
    },
    'electrons': {
        'conv_thr': 1e-8,
        'mixing_beta': 0.7,
        'diagonalization': 'david'
    }
}

command = 'mpirun -np 4 pw.x'
profile = EspressoProfile(command, pseudo_dir=PSEUDO_DIR)

calc = Espresso(profile=profile, 
                pseudopotentials=pseudopotentials,
                input_data=input_data,
                kpts=(1, 1, 1))  # Gamma-point only for molecules[cite: 9]

ch4.calc = calc

# 3. Run the Geometry Optimization
dyn = BFGS(ch4, trajectory='ch4_opt.traj', logfile='ch4_opt.log')
dyn.run(fmax=0.01)

# 4. Print optimized geometry
print("Optimized Methane Geometry:")
print(f"C-H1 Bond Length (A): {ch4.get_distance(0, 1):.4f}")
print(f"H-C-H1 Angle (deg):   {ch4.get_angle(1, 0, 2):.2f}")
print(f"Final Energy (eV):    {ch4.get_total_energy():.4f}")
