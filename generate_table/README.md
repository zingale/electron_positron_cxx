# `generate_table`

This generates a table of Helmholtz free energy and other
thermodynamic quantities in the format used by the Timmes & Swesty
(2000) "helmeos" equation of state.  To use this table in that EOS,
you should adjust the minimum / maximum density ranges and number of
points to match what is used in the generation.

There are 4 separate tables written:

* Helmholtz free energy -- this contains 9 values at each table point:
  F, ∂F/∂ρ, ∂F/∂T, ∂²F/∂ρ², ∂²F/∂T², ∂²F/∂ρ∂T, ∂³F/∂ρ²∂T, ∂³F/∂ρ∂T², ∂⁴F/∂ρ²∂T²

* pressure derivative -- this contains 4 values at each table point:
  ∂p/∂ρ, ∂²p/∂ρ², ∂²p/∂ρ∂T, ∂³p/∂ρ²∂T

* degeneracy parameter -- this contains 4 values at each table point:
  η, ∂η/∂ρ, ∂η/∂T, ∂²η/∂ρ∂T

* electron + positron number density -- this contains 4 values at each
  table point:
  n, ∂n/∂ρ, ∂n/∂T, ∂²n/∂ρ∂T

