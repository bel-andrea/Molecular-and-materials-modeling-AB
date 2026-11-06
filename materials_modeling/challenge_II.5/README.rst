Al electronic properties
=========================

The selected experimental structure for aluminum was obtained from the American Mineralogist Crystal Structure Database (AMCSD), with the code 0011136.

https://www.rruff.net/odr/amcsd#/odr/view/724369/2187/eyI3MTk3IjoiYWx1bWludSIsImR0X2lkIjoiNzcxIiwiaWdub3JlIjoiNzcxXzczNl84MTcsNzcxXzczNl84MTgsNzcxXzczNl84MTksNzcxXzczNl83MzQifQ/1

Process & files
----------------

Firstly, the convergence test was run using data from the experimental structure (``convergence_test_al.py``& ``convergence_test_al.out``).
Calculations converged at *ecutwfc* = 60 & *K-points* = (10, 10, 10),

*Note:* for this challenge, the convergence criteria was adjusted to ΔE < 0.05 meV/atom, as calculations were attempted using the previous criteria (ΔE < 0.01 meV/atom), which did not show convergence after running for 5+ hours and up to *ecutwfc* = 120.

Then, cell relaxation was performed considering *ecutwfc* & *K-points* identified during the convergence test (``cell_relaxation_al.py``& ``cell_relaxation_al.out``). The relaxed structure was saved to ``final_relaxed_structure_al.vasp``.

Finally, electronic properties were computed (``electronic_properties_al.py``& ``electronic_properties_al.out``).
