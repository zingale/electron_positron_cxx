# Electron-Positron Equation of State

This is an electron-positron equation of state that computes the
number density, pressure, and energy of a Fermi gas, along with the
first derivatives.

This uses the method of Gong et al. 2001 to compute the Fermi-Dirac
integrals, and follows the notation from Timmes and Arnett 1999 for
the forms of the thermodynamic quantities.

To get good results, this uses 128-bit precision, relying on GCC's
`__float128` and the quadmath library.


## Requirements

This needs a C++23 compiler.  For GCC, you need GCC >= 15.1

## Options

The precision can be changed via the `PRECISION` make variable, e.g.,

* `make PRECISION=FLOAT128` builds with 128-bit precision (using the
  `__float128` data type.

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

In `tests/` type:

```
make
```

optionally, you can switch from `_float128` to `long double` by
editing the `GNUmakefile`.

## clang-tidy

To check the codebase with `clang-tidy`, build as:

```
make USE_CLANG_TIDY=TRUE
```
