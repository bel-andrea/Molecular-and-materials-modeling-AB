## Software checks\
  All checks were performed employing a Macbook Air M4 Chip\

### NWChem\
_Version 7.3.1.\_
  **_Primitive runs:_ ch3_zora_b3lyp_prop.out & ch3_zora_b3lyp_prop.outN2\**
NWChem 7.3.1\
Nthr 6   |   6-311G   |   Total times cpu: 7.3 s\
Nthr 2   |   6-311G   |   Total times cpu: 9.2 s\
**_Python-driven run:_  run_nwchem.py_logfileAB\**
Nproc 2   |   Wall time: 8.90 s   |   Speedup: 1.00   |   Efficiency: 50.0\
Nproc 4   |   Wall time: 7.00 s   |   Speedup: 1.27   |   Efficiency: 31.8\

### Quantum Espresso\
_PWSCF version 7.5_\
**_Primitive run:_ Tl2.so_scf.in_out_MPI8_AB\**
MPI 8   |   Wall time: 20.76 s\
**_Python-driven benchmark:_ tl_bulk_scf_04.py_logfile_AB\**
MPI 1   |   Wall time: 39.59 s\
MPI 2   |   Wall time: 24.63 s\
MPI 4   |   Wall time: 16.54 s\
MPI 8   |   Wall time: 21.14 s\
MPI 10  |   Wall time:  22.26 s\
**_AI-improved benchmark:_ tl_bulk_scf_05.py_MPI4_AB & tl_bulk_scf_05.py_MPI8_AB\**
MPI 4   |   Wall time: 16.80 s\
MPI 8   |   Wall time: 21.05 s\

### MOPAC\
_Version 23.3.5\_
** dna_pm7_threads.out\**
Nthr 6   |   Total Job Time: 28.27 s
