# Electron-Positron Equation of State

This is an electron-positron equation of state that computes the
number density, pressure, energy, and entropy of a Fermi gas (both the
electron and positron contributions), along with the first- and
second-derivatives with respect to ρ and T.

This uses the method of Gong et al. 2001 to compute the Fermi-Dirac
integrals, and follows the notation from Timmes and Arnett 1999 for
the forms of the thermodynamic quantities.

To get good results, this uses 128-bit precision, relying on GCC's
`__float128` implementation.


## Requirements

This needs a C++23 compiler.  For GCC, you need GCC >= 15.1

## Options

The precision can be changed via the `PRECISION` make variable, e.g.,

* `make PRECISION=FLOAT128` builds with 128-bit precision (using the
  `__float128` data type.

  Note: at the moment, this is only supported with GCC, and not
  CLANG because of an incompatibility with 128-bit floats and
  `std::println()` in C++ under CLANG.

* `make PRECISION=LONG_DOUBLE` builds with a `long double`, which is
  80 bits on x86 architectures.

* `make PRECISION=DOUBLE` will build with 64-bit precision, using
   `double`.

The number of quadrature points used in the integration can be set via
`QUAD_PTS`, e.g., as:

```
make QUAD_PTS=100
```

for 100 points.  Valid options are `50`, `100`, `200` (the default),
and `400`.

Be sure to do

```
make clean
```

before building with any different options.


## Driver

A basic driver that takes density, temperature, and Y_e is in `eos/`.
To build, type:

```
make
```

An example usage is:

```
$ echo 1.e-2 1.e9 0.5 | ./eos 
Enter rho, T, Ye (space-separated): 

 ρ =            0.01       T  =           1e+09      Yₑ  =             0.5

 β =      0.16863701
 η =      -5.9298938    ∂η/∂ρ  =   0.00027751945    ∂η/∂T  =   5.9298752e-09
                       ∂²η/∂ρ² =  -2.1269302e-13   ∂²η/∂T² =  -1.1859592e-17   ∂²η/∂ρ∂T =  -2.1338024e-12

n⁻ =   5.4294165e+26   ∂n⁻/∂ρ  =   1.5055394e+23   ∂n⁻/∂T  =   4.1771487e+18
                      ∂²n⁻/∂ρ² =   4.1713392e+19  ∂²n⁻/∂T² =   2.4851053e+10  ∂²n⁻/∂ρ∂T =  -3.2112223e+09
n⁺ =   5.4293864e+26   ∂n⁺/∂ρ  =   -1.505531e+23   ∂n⁺/∂T  =   4.1771487e+18
                      ∂²n⁺/∂ρ² =   4.1713392e+19  ∂²n⁺/∂T² =   2.4851053e+10  ∂²n⁺/∂ρ∂T =  -3.2112223e+09

p⁻ =   7.4991791e+19   ∂p⁻/∂ρ  =   2.0803187e+16   ∂p⁻/∂T  =   6.5212244e+11
                      ∂²p⁻/∂ρ² =   5.7685628e+12  ∂²p⁻/∂T² =       4590.1194  ∂²p⁻/∂ρ∂T =        20901188
p⁺ =   7.4991375e+19   ∂p⁺/∂ρ  =  -2.0803071e+16   ∂p⁺/∂T  =   6.5212202e+11
                      ∂²p⁺/∂ρ² =   5.7685627e+12  ∂²p⁺/∂T² =       4590.1194  ∂²p⁺/∂ρ∂T =       -20901959

e⁻ =   1.3262018e+22   ∂e⁻/∂ρ  =  -1.3261981e+24   ∂e⁻/∂T  =    1.170244e+14
                      ∂²e⁻/∂ρ² =   2.6523962e+26  ∂²e⁻/∂T² =       840753.43  ∂²e⁻/∂ρ∂T =  -1.1702436e+16
e⁺ =   1.0216387e+23   ∂e⁺/∂ρ  =  -1.0216415e+25   ∂e⁺/∂T  =   8.0099949e+14
                      ∂²e⁺/∂ρ² =    2.043283e+27  ∂²e⁺/∂T² =       4909917.4  ∂²e⁺/∂ρ∂T =  -8.0099953e+16

s⁻ =   6.5212383e+13   ∂s⁻/∂ρ  =  -6.5212223e+15   ∂s⁻/∂T  =       459011.83
                      ∂²s⁻/∂ρ² =   1.3042445e+18  ∂²s⁻/∂T² =    0.0024163237  ∂²s⁻/∂ρ∂T =       -45901194
s⁺ =   6.5212063e+13   ∂s⁺/∂ρ  =  -6.5212223e+15   ∂s⁺/∂T  =       459012.06
                      ∂²s⁺/∂ρ² =   1.3042445e+18  ∂²s⁺/∂T² =    0.0024163232  ∂²s⁺/∂ρ∂T =       -45901194
```

## Tests

There are a large number of tests that exercise different parts of the
algorithm.  In `tests/` type:

```
make
```

You can optionally set the precision and number of quadrature points
as described above.  The `tests/README.md` describes the basic
functionality of the tests.


## Derivations

The derivation of the derivatives of the different thermodynamic
quantities is provided by the notebook `degenerate_derivatives.ipynb`,
using SymPy.

These are the expressions that are coded up in `electron_positron.H`.


## Generating quadrature points / weights

The Gauss-Legendre and Gauss-Leguerre quadrature roots and weights are
generated via the notebook `generate_quadrature_weights.ipynb` using
SymPy.  This will directly write the C++ header file with the desired
number of quadrature points.


## clang-tidy

To check the codebase with `clang-tidy`, build as:

```
make USE_CLANG_TIDY=TRUE PRECISION=LONG_DOUBLE
```

At the moment, building with `PRECISION=FLOAT128` with
`clang-tidy` does not work.
