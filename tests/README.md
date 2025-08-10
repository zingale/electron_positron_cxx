# tests

The following tests exercise different parts of the EOS algorithm:

* `test_bounds.cpp` : This tests the routines in
  `degeneracy_parameter_bounds.H` that give the bounds for the root
  finding for η.

* `test_breakpoints.cpp` : This tests the routines that compute the
  breakpoints for the quadrature that are defined in Gong et al. 2001.

* `test_brent.cpp` : This checks the implementation of Brent's method
  for root finding in `brent.H`, using the test problem from Brent's
  original paper.

* `test_difference.cpp` : This tests the different finite-difference
  methods implemented in `difference_utils.H` using the test problem
  from Ridders' original paper.

* `test_electron_positron.cpp` : This compares the derivative computed
  in the EOS (`electron_positron.H`) to finite-difference
  approximations.

* `test_eta_derivs.cpp` : This tests the density and temperature
  first- and second-derivatives of η by comparing to finite-difference
  approximations.

* `test_fermi.cpp` : This tests the Fermi integral first- and
  second-derivatives computed via quadrature by comparing to
  finite-difference approximations.

* `test_maxwell.cpp` : This checks the Maxwell relations for
  thermodynamic consistency.



