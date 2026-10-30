Software checks
===============

All checks were performed using a MacBook Air with an M4 chip.

NWChem
------

*Version 7.3.1*

**Primitive runs:** ``ch3_zora_b3lyp_prop.out`` and ``ch3_zora_b3lyp_prop.outN2``

=====  ========  ==============
Nthr   Basis     Total CPU time
=====  ========  ==============
6      6-311G    7.3 s
2      6-311G    9.2 s
=====  ========  ==============

**Python-driven run:** ``run_nwchem.py_logfileAB``

=====  ==========  =======  ==============
Nproc  Wall time   Speedup  Efficiency (%)
=====  ==========  =======  ==============
2      8.90 s      1.00     50.0
4      7.00 s      1.27     31.8
=====  ==========  =======  ==============

Quantum Espresso
----------------

*PWSCF version 7.5*

**Primitive run:** ``Tl2.so_scf.in_out_MPI8_AB``

=====  ==========
MPI    Wall time
=====  ==========
8      20.76 s
=====  ==========

**Python-driven benchmark:** ``tl_bulk_scf_04.py_logfile_AB``

=====  ==========
MPI    Wall time
=====  ==========
1      39.59 s
2      24.63 s
4      16.54 s
8      21.14 s
10     22.26 s
=====  ==========

**AI-improved benchmark:** ``tl_bulk_scf_05.py_MPI4_AB`` and ``tl_bulk_scf_05.py_MPI8_AB``

=====  ==========
MPI    Wall time
=====  ==========
4      16.80 s
8      21.05 s
=====  ==========

MOPAC
-----

*Version 23.3.5*

**Run:** ``dna_pm7_threads.out``

=====  ==============
Nthr   Total job time
=====  ==============
6      28.27 s
=====  ==============
