# Electron-Positron Equation of State

This is an electron-positron equation of state that computes the
number density, pressure, and energy of a Fermi gas, along with the
first derivatives.

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

for 100 points.  Valid options are `50`, `100`, `200` (the default), and `400`.

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
