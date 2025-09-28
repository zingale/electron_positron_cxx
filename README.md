[![DOI](https://zenodo.org/badge/1026869302.svg)](https://doi.org/10.5281/zenodo.16942327)

[![Test compilation](https://github.com/zingale/electron_positron_cxx/actions/workflows/compile_action.yml/badge.svg)](https://github.com/zingale/electron_positron_cxx/actions/workflows/compile_action.yml)
[![Test MacOS compilation](https://github.com/zingale/electron_positron_cxx/actions/workflows/compile_macos_action.yml/badge.svg)](https://github.com/zingale/electron_positron_cxx/actions/workflows/compile_macos_action.yml)
[![Run tests](https://github.com/zingale/electron_positron_cxx/actions/workflows/run_tests.yml/badge.svg)](https://github.com/zingale/electron_positron_cxx/actions/workflows/run_tests.yml)

# Electron-Positron Equation of State

This is an electron-positron equation of state that computes the
number density, pressure, energy, and entropy of a Fermi gas (both the
electron and positron contributions), along with the first-, second-,
and third-derivatives with respect to ρ and T.

This uses the method of Gong et al. 2001 to compute the Fermi-Dirac
integrals, and follows the notation from Timmes and Arnett 1999 for
the forms of the thermodynamic quantities.

To get good results, this can use 128- or 256-bit precision.

Overall, this code is slow, and it is mainly intended to be used
to tabulate EOS properties.


## Requirements

* a C++23 compiler.  For GCC, you need GCC >= 15.1.

* [libquadmath](https://gcc.gnu.org/onlinedocs/libquadmath/) for 128-bit precision

* [boost](https://www.boost.org/) for 256-bit precision

* GNU make >= 4.0.

> [!NOTE]
> On a Mac, you will need to install make via brew, since Apple's
> version is more than 10 years old.  Then use `gmake` instead of
> `make` in the building examples below.


## Root-finding for degeneracy parameter

The first step in finding the thermodynamic state is to solve for the
degeneracy parameter, η, via a root-finding process that ensures
charge neutrality (the number density of electrons from the mass
density must equal the number density of electrons minus the number
density of positrons).

We use Brent's method to do this solve, which requires a bracket that
contains the roots.  We precompute the degeneracy parameter on a grid
of ρYₑ and T, and use this to provide bounds.  The code
`generate_eta/generate_etas.cpp` will compute this grid.

## Floating point precision

Several different floating point standards are supported.  The entire
EOS is templated on a `real_t` type, which enables any real data type
to be used.  The precision can be changed via the `PRECISION`, as
described below.

### 256-bit floating point

256-bit is enabled via the Boost library, `boost::multiprecision`.
The floating point properties (as output by `tests/test_precision`)
are:

```
size of real_t is 64 bytes
machine epsilon is 9.055679078826712368e-72
minimum exponent is 10**-78862
maximum exponent is 10**78862
```

This is built as:

```
make PRECISION=BOOST256
```


### 128-bit floating point

128-bit is enabled via the GCC `__float128` data type and the quadmath
library.  The floating point properties are:

```
size of real_t is 16 bytes
machine epsilon is 1.9259299443872358530559779425849273e-34
minimum exponent is 10**-4931
maximum exponent is 10**4932
```

This is built as:

```
make PRECISION=FLOAT128
```

> [!NOTE]
> On x86 architectures, this is supported both with GCC and LLVM/clang.
> On Apple ARM processors, apple-clang does not support `__float128`,
> so you need to install the latest GCC via `brew` and then build with
> `CXX=g++-15` on the compile line.


### 80-bit floating point

80-bit precision is enabled via the x87 extensions on Intel
architectures, and uses a `long double` datatype.  The floating
point properties are:

```
size of real_t is 16 bytes
machine epsilon is 1.084202172485504434e-19
minimum exponent is 10**-4931
maximum exponent is 10**4932
```

This is built as:

```
make PRECISION=LONG_DOUBLE
```

> [!NOTE]
> On Macs with ARM processors, a `long double` is the same as
> a `double` and will give only 64-bit precision.


### 64-bit floating point

64-bit precision uses the normal `double` data type.  The
floating point properties are:

```
size of real_t is 8 bytes
machine epsilon is 2.220446049250313e-16
minimum exponent is 10**-307
maximum exponent is 10**308
```

This is build as:

```
make PRECISION=DOUBLE
```


## Quadrature

The number of quadrature points used for each subinterval in the
overall integration can be set via `QUAD_PTS`, e.g., as:

```
make QUAD_PTS=100
```

for 100 points.  Valid options are `20`, `50`, `100`, `200` (the default),
`400`, and `800`.

Be sure to do

```
make clean
```

before building with any different options.


### Generating quadrature points / weights

The Gauss-Legendre and Gauss-Leguerre quadrature roots and weights are
generated via the notebook `tools/generate_quadrature_weights.ipynb` using
SymPy.  This will directly write the C++ header file with the desired
number of quadrature points.


## Driver

A basic driver that takes density, temperature, and Y_e is in `eos/`.
To build, type:

```
make
```

An example usage is:

```
Enter rho, T, Ye (space-separated): 

 ρ =            0.01       T  =           1e+09      Yₑ  =             0.5

 β =      0.16863701

degeneracy parameter:
 η =      -5.9298938    ∂η/∂ρ  =   0.00027751945    ∂η/∂T    =   5.9298752e-09
                       ∂²η/∂ρ² =  -2.1269302e-13   ∂²η/∂ρ∂T  =  -2.1338024e-12  ∂²η/∂T²   =  -1.1859592e-17
                       ∂³η/∂ρ³ =  -2.1269302e-11   ∂³η/∂ρ²∂T =   4.9121274e-21  ∂³η/∂ρ∂T² =   2.0135641e-20  ∂³η/∂T³ =   3.5577151e-26  

number density:
  n⁻ =   5.4294165e+26   ∂n⁻/∂ρ  =   1.5055394e+23   ∂n⁻/∂T    =   4.1771487e+18
                        ∂²n⁻/∂ρ² =   4.1713392e+19  ∂²n⁻/∂ρ∂T  =  -3.2112223e+09  ∂²n⁻/∂T²   =   2.4851053e+10
                        ∂³n⁻/∂ρ³ =  -9.6064995e+10  ∂³n⁻/∂ρ²∂T =  -3.2112223e+11  ∂³n⁻/∂ρ∂T² =       30.311409  ∂³n⁻/∂T³ =       100.02177
  n⁺ =   5.4293864e+26   ∂n⁺/∂ρ  =   -1.505531e+23   ∂n⁺/∂T    =   4.1771487e+18
                        ∂²n⁺/∂ρ² =   4.1713392e+19  ∂²n⁺/∂ρ∂T  =  -3.2112223e+09  ∂²n⁺/∂T²   =   2.4851053e+10
                        ∂³n⁺/∂ρ³ =  -9.6064995e+10  ∂³n⁺/∂ρ²∂T =  -3.2112223e+11  ∂³n⁺/∂ρ∂T² =       30.311409  ∂³n⁺/∂T³ =       100.02177
  n  =   1.0858803e+27   ∂n/∂ρ   =   8.3426784e+17   ∂n/∂T     =   8.3542974e+18
                        ∂²n/∂ρ²  =   8.3426784e+19  ∂²n/∂ρ∂T   =  -6.4224446e+09  ∂²n/∂T²    =   4.9702106e+10
                        ∂³n/∂ρ³  =  -1.9212999e+11  ∂³n/∂ρ²∂T  =  -6.4224446e+11  ∂³n/∂ρ∂T²  =       60.622818  ∂³n/∂T³  =       200.04354

pressure:
  p⁻ =   7.4991791e+19   ∂p⁻/∂ρ  =   2.0803187e+16   ∂p⁻/∂T    =   6.5212244e+11
                        ∂²p⁻/∂ρ² =   5.7685628e+12  ∂²p⁻/∂ρ∂T  =        20901188  ∂²p⁻/∂T²   =       4590.1194
                        ∂³p⁻/∂ρ³ =       3896696.4  ∂³p⁻/∂ρ²∂T =      -38584.989  ∂³p⁻/∂ρ∂T² =   0.00056995441  ∂³p⁻/∂T³ =   2.4163235e-05
  p⁺ =   7.4991375e+19   ∂p⁺/∂ρ  =  -2.0803071e+16   ∂p⁺/∂T    =   6.5212202e+11
                        ∂²p⁺/∂ρ² =   5.7685627e+12  ∂²p⁺/∂ρ∂T  =       -20901959  ∂²p⁺/∂T²   =       4590.1194
                        ∂³p⁺/∂ρ³ =      -3923222.8  ∂³p⁺/∂ρ²∂T =      -38584.988  ∂³p⁺/∂ρ∂T² =   -0.0005633577  ∂³p⁺/∂T³ =   2.4163235e-05
  p  =   1.4998317e+20   ∂p/∂ρ   =   1.1537126e+11   ∂p/∂T     =   1.3042445e+12
                        ∂²p/∂ρ²  =   1.1537126e+13  ∂²p/∂ρ∂T   =      -771.69977  ∂²p/∂T²    =       9180.2389
                        ∂³p/∂ρ³  =      -26526.423  ∂³p/∂ρ²∂T  =      -77169.977  ∂³p/∂ρ∂T²  =   6.5967097e-06  ∂³p/∂T³  =   4.8326469e-05

specific internal energy:
  e⁻ =   1.3262018e+22   ∂e⁻/∂ρ  =  -1.3261981e+24   ∂e⁻/∂T    =    1.170244e+14
                        ∂²e⁻/∂ρ² =   2.6523962e+26  ∂²e⁻/∂ρ∂T  =  -1.1702436e+16  ∂²e⁻/∂T²   =       840753.43  
                        ∂³e⁻/∂ρ³ =  -7.9571885e+28  ∂³e⁻/∂ρ²∂T =   2.3404872e+18  ∂³e⁻/∂ρ∂T² =       -84075342  ∂³e⁻/∂T³ =    0.0045713038
  e⁺ =   1.0216387e+23   ∂e⁺/∂ρ  =  -1.0216415e+25   ∂e⁺/∂T    =   8.0099949e+14
                        ∂²e⁺/∂ρ² =    2.043283e+27  ∂²e⁺/∂ρ∂T  =  -8.0099953e+16  ∂²e⁺/∂T²   =       4909917.4
                        ∂³e⁺/∂ρ³ =  -6.1298489e+29  ∂³e⁺/∂ρ²∂T =   1.6019991e+19  ∂³e⁺/∂ρ∂T² =  -4.9099174e+08  ∂³e⁺/∂T³ =      0.02094908
  e  =   1.1542588e+23   ∂e/∂ρ   =  -1.1542613e+25   ∂e/∂T     =   9.1802389e+14
                        ∂²e/∂ρ²  =   2.3085226e+27  ∂²e/∂ρ∂T   =  -9.1802389e+16  ∂²e/∂T²    =       5750670.8
                        ∂³e/∂ρ³  =  -6.9255678e+29  ∂³e/∂ρ²∂T  =   1.8360478e+19  ∂³e/∂ρ∂T²  =  -5.7506708e+08  ∂³e/∂T³  =     0.025520384

specific entropy:
  s⁻ =   6.5212383e+13   ∂s⁻/∂ρ  =  -6.5212223e+15   ∂s⁻/∂T  =       459011.83
                        ∂²s⁻/∂ρ² =   1.3042445e+18  ∂²s⁻/∂T² =    0.0024163237  ∂²s⁻/∂ρ∂T =       -45901194
  s⁺ =   6.5212063e+13   ∂s⁺/∂ρ  =  -6.5212223e+15   ∂s⁺/∂T  =       459012.06
                        ∂²s⁺/∂ρ² =   1.3042445e+18  ∂²s⁺/∂T² =    0.0024163232  ∂²s⁺/∂ρ∂T =       -45901194
  s  =   1.3042445e+14   ∂s/∂ρ   =  -1.3042445e+16   ∂s/∂T   =       918023.89
                        ∂²s/∂ρ²  =   2.6084889e+18  ∂²s/∂T²  =    0.0048326469  ∂²s/∂ρ∂T  =       -91802389

```

## Generating an EOS table

The code in `generate_table/` will compute a
Helmholtz-free-energy-based table in the format used by the popular
Timmes & Swesty (2000) helmeos.  No differencing is used to compute
the high-order derivatives of the Helmholtz free energy, rather we
directly compute all derivatives via integration of the
(differentiated) Fermi-Dirac integrand.

The driver uses OpenMP to parallelize over the points in the table.

## Tests

There are a large number of tests that exercise different parts of the
algorithm.  In `tests/` type:

```
make
```

You can optionally set the precision and number of quadrature points
as described above.  The `tests/README.md` describes the basic
functionality of the tests.

Many of these will compare a difference approximation to the
derivatives to the derivative computed by the EOS by integration.


## Derivations

The derivation of the derivatives of the different thermodynamic
quantities is provided by the notebook
`notes/degenerate_derivatives.ipynb`, using SymPy.

These are the expressions that are coded up in `electron_positron.H`.



## clang-tidy

To check the codebase with `clang-tidy`, build as:

```
make USE_CLANG_TIDY=TRUE PRECISION=LONG_DOUBLE
```

At the moment, building with `PRECISION=FLOAT128` with
`clang-tidy` does not work.
