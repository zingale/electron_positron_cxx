# tests

The following tests exercise different parts of the EOS algorithm:

* `test_bounds.cpp` : test the routines in
  `degeneracy_parameter_bounds.H` that give the bounds for the root
  finding for η.

* `test_breakpoints.cpp` : test the routines that compute the
  breakpoints for the quadrature that are defined in Gong et al. 2001.

* `test_brent.cpp` : check the implementation of Brent's method for
  root finding in `brent.H`, using the test problem from Brent's
  original paper.

* `test_difference.cpp` : test the different finite-difference methods
  implemented in `difference_utils.H` using the test problem from
  Ridders' original paper.

* `test_eos_energy.cpp` : test the first- and second-derivatives of
  specific energy (with respect to ρ and T) by comparing to
  finite-difference approximations.

* `test_eos_entropy.cpp` : test the first- and second-derivatives of
  specific entropy (with respect to ρ and T) by comparing to
  finite-difference approximations.

* `test_eos_number_density.cpp` : test the first- and
  second-derivatives of number density (with respect to ρ and T) by
  comparing to finite-difference approximations.

* `test_eos_pressure.cpp` : test the first- and second-derivatives of
  pressure (with respect to ρ and T) by comparing to finite-difference
  approximations.

* `test_eta_derivs.cpp` : test the first- and second-derivatives of η
  (with respect to ρ and T) by comparing to finite-difference
  approximations.

* `test_fermi.cpp` : test the Fermi integral first- and
  second-derivatives (with respect to η and β) computed via quadrature
  by comparing to finite-difference approximations.

* `test_maxwell.cpp` : check the Maxwell relations for thermodynamic
  consistency.

* `test_n_derivs.cpp` : check the first- and second-derivatives (with
  respect to η and β) of the number density.  These are the building
  blocks of the η derivatives that are then used for all the other
  thermodynamic derivatives.

* `test_precision.cpp` : print out the details of the floating point
  representation and test the printing tools in `util.H`.

