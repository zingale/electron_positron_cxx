# Electron-Positron Equation of State

This is an electron-positron equation of state that computes the
number density, pressure, and energy of a Fermi gas, along with the
first derivatives.

This uses the method of Gong et al. 2001 to compute the Fermi-Dirac
integrals, and follows the notation from Timmes and Arnett 1999 for
the forms of the thermodynamic quantities.

To get good results, this uses 128-bit precision, relying on GCC's
`__float128` and the quadmath library.

The precision can be selected by editing the `Make.eos` file.

## Requirements

This needs a C++23 compiler.  For GCC, you need GCC >= 15.1

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
