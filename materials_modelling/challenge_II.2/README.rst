Cutoff and K-points convergence test
=====================================

Task 2.2
---------
*QE calculation with tstress = .true. and tprnfor = .true. :* results found in ``Si_force_stress.out``

Task 2.3
---------
*Manually shifting the atoms:* the X position was shifted by 0.02 Å, and the convergence test was run again. 
Results are found in the ``convergence_test_shifted.out`` file inside this folder. Calculations converged with *ecutwfc* = 65 and *K-points* = (15, 15, 15).

Challenge
----------
The convergence test was run with three different pseudopotentials:

**1:** ``convergence_test_si.out`` optimized norm-conserving Vanderbilt Pseudopotential (ONCV). Inside ``class_exercise`` folder.

**2:** ``convergence_test_si_pbe_rrkjus.out`` Rappe-Rabe-Kaxiras-Joannopoulos Ultrasoft Pseudopotential (USPP).

**3:** ``convergence_test_si_pbesol_kjpaw.out`` Kresse-Joubert Projector Augmented-Wave (PAW) Pseudopotential.

======   =======   ============
Pseudo   ecutwfc   K-point
======   =======   ============
ONCV     65        (15, 15, 15)
------   -------   ------------
USPP     90        (15, 15, 15)
------   -------   ------------
PAW      110       (16, 16, 16)
======   =======   ============

