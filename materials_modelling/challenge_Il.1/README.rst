Total energy calculation
========================

Pseudopotentials
-----------------
QE calculations were run using three different pseudopotentials for Si. Files are organized in three different folders, for each pseudo used:

**1:** ``norm_conserving`` optimized norm-conserving Vanderbilt Pseudopotential (ONCV).

**2:** ``pbe_rrkjus`` Rappe-Rabe-Kaxiras-Joannopoulos Ultrasoft Pseudopotential (USPP).

**3:** ``pbesol_kjpaw`` Kresse-Joubert Projector Augmented-Wave (PAW) Pseudopotential.

The employed pseudos are also found under the ``pseudos`` folder.

Modification of *ecutwfc*
-------------------------
For each pseudopotential, QE calculations were run using values of *ecutwfc* ranging from 20 to 50, in intervals of 5.
The results were then plotted and saved as ‘‘.png‘‘ files, located inside the folders for each pseudopotential.
Overall, *ecutwfc*=20 resulted in the highest energy values.
