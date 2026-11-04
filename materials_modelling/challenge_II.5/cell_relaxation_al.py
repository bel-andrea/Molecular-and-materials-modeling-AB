#!/usr/bin/env python
from ase import Atoms
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.io import write
from ase.optimize import BFGS
import numpy as np
import sys, os
from ase.filters import UnitCellFilter 
sys.stdout.flush()

# ==============================================
# ASE full relaxation (unitcell + atom) - Aluminum
# ==============================================

# ==============================================
# 1. Quantum Espresso Input Parameters
# ==============================================
pseudopotentials = {'Al': 'Al.upf'} 

input_data = {
    'control': {
        'calculation': 'scf',
        'prefix': 'al',
        'outdir': './tmp',
        'verbosity': 'low',
        'tstress': True,
        'tprnfor': True
    },
    'system': {
        'ecutwfc': 60,                  # Updated based on convergence test
        'occupations': 'smearing',
        'smearing': 'gauss', # Standard smearing for metals
        'degauss': 0.01,
        'ibrav': 0,
        'nat': 4,                       # 4 atoms in conventional FCC cell
        'ntyp': 1
    },
    'electrons': {
        'conv_thr': 1.0e-10
    }
}

# ==============================================
# 2. Atomic Structure (AMCSD 001136, FCC Al)
# ==============================================
a = 4.04958  # Experimental lattice parameter in Å

atoms = Atoms(
    symbols=['Al'] * 4,
    scaled_positions=[
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5]
    ],
    cell=[
        [a, 0.0, 0.0],
        [0.0, a, 0.0],
        [0.0, 0.0, a]
    ],
    pbc=[True, True, True]
)

# ==============================================
# 3. Calculator Configuration
# ==============================================
qe_bin = "/Users/andreabelendez/anaconda3/envs/molmatmodel/bin/"
pw_command = f'mpirun -np 8 {qe_bin}/pw.x'

pw_profile = EspressoProfile(
    command=pw_command,
    pseudo_dir='../.././'
)

calc = Espresso(
    profile=pw_profile,
    pseudopotentials=pseudopotentials,
    input_data=input_data,
    kpts=(10, 10, 10)  # Balanced k-point grid for metallic Al cell relaxation
)

atoms.calc = calc

# ==============================================
# 4. Cell Relaxation
# ==============================================
print("\nStarting Aluminum cell relaxation by ASE", flush=True)

ucf = UnitCellFilter(atoms)

from ase.io.trajectory import Trajectory
traj = Trajectory('relaxation_al.traj', 'w', atoms)

opt = BFGS(ucf, trajectory=traj, logfile='relaxation_al.log')

try:
    opt.run(fmax=0.01) # Convergence criterion: max force < 0.01 eV/Å
        
    final_energy = atoms.get_potential_energy()
    forces = opt.atoms.get_forces()
    qe_style_forces = atoms.get_forces() 
    stress = atoms.get_stress()
     
    force_norms = np.linalg.norm(forces, axis=1)
    max_force = np.max(force_norms)
    max_qe_force = np.max(np.abs(qe_style_forces))
    pressure = -np.sum(stress[:3]) * 1602.1766208 / 3

    print("\nFinal results (Aluminum)", flush=True)    
    print(f"  Total Energy                  : {final_energy:>12.6f} eV", flush=True)
    print(f"  ASE-style max force (norm)    : {max_force:>8.6f} eV/Å", flush=True)
    print(f"  QE-style max force            : {max_qe_force:>8.6f} eV/Å", flush=True)
    print(f"  Pressure                      : {pressure:>8.6f} kbar", flush=True)

    write('final_relaxed_structure_al.vasp', atoms, format='vasp', direct=False)
    print("\nFinal relaxed structure saved to: final_relaxed_structure_al.vasp", flush=True)

except Exception as e:
    print(f"\nRelaxation failed: {str(e)}", flush=True)
    print("Check output files for error details", flush=True)
finally:
    traj.close()
    print("\nRelaxation complete", flush=True)
