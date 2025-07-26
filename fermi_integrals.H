#ifndef FERMI_INTEGRALS_H
#define FERMI_INTEGRALS_H

// Classes that enable the computation of Fermi-Dirac integrals to
// high precision.  This uses the methods of Aparicio (1998) and Gong
// (2001), described in
// https://iopscience.iop.org/article/10.1086/313121 and
// https://www.sciencedirect.com/science/article/abs/pii/S001046550100145X?via%3Dihub.

// The method consist on breaking down the Fermi-integrals (eta,
// theta) integration domain into four subintervals and apply
// Gauss-Legendre or the Gauss laguerre quadrature methods with 200
// points over each domain to guarantee accuracy beyond
// double-precision.

// In the first subinterval, a x=z^2 change of variable is applied to
// overcome the kernel singularity near the origin.

#include <iostream>
#include <cstdlib>
#include <array>
#include <tuple>
#include <cmath>

#include "quadrature_weights.H"

template <typename T>
class BreakPoints {

    // The break points used in splitting the integral from [0, inf]
    // into separate integrals over smaller domains.  These are
    // described in Gong et al. 2001.  Different choices of the break
    // points are made depending on whether we are computing the
    // integral itself or one of its derivatives with respect to eta
    // (the degeneracy parameter).

    // this is roughly the maximum value of x such that std::exp(x) does
    // not overflow

    T max_exponent{};

    T D{};
    T sigma{};
    std::array<T, 3> a{};
    std::array<T, 3> b{};
    std::array<T, 3> c{};
    std::array<T, 3> d{};
    std::array<T, 3> e{};

public:

    explicit BreakPoints(int eta_order) :
        max_exponent{std::trunc(std::log(std::numeric_limits<T>::max()))}
    {

        if (eta_order == 0) {
            D = 3.36091L;
            sigma = 9.11856e-2L;

            a[0] = 6.77740L;
            a[1] = 3.76010L;
            a[2] = 7.56690L;

            b[0] = 1.14180L;
            b[1] = 9.37188e-2L;
            b[2] = 1.16953L;

            c[0] = 2.98255L;
            c[1] = 2.10635e-2L;
            c[2] = 7.54162e-1L;

            d[0] = 0.0;
            d[1] = 3.10836e1L;
            d[2] = 6.65585L;

            e[0] = 0.0L;
            e[1] = 1.00557L;
            e[2] = -1.28190e-1L;

        } else if (eta_order == 1) {

            D = 4.99551L;
            sigma = 9.11856e-2L;

            a[0] = 6.77740L;
            a[1] = 3.76010L;
            a[2] = 7.56690L;

            b[0] = 1.14180L;
            b[1] = 9.37188e-2L;
            b[2] = 1.16953L;

            c[0] = 2.98255L;
            c[1] = 2.10635e-2L;
            c[2] = 7.54162L;

            d[0] = 0.0L;
            d[1] = 3.95015e1L;
            d[2] = 7.64734L;

            e[0] = 0.0L;
            e[1] = 1.00557L;
            e[2] = -1.28190e-1L;

        } else if (eta_order == 2) {

            D = 3.93830L;
            sigma = 9.11856e-2L;

            a[0] = 6.77740L;
            a[1] = 3.76010L;
            a[2] = 7.56690L;

            b[0] = 1.14180L;
            b[1] = 9.37188e-2L;
            b[2] = 1.16953L;

            c[0] = 2.98255L;
            c[1] = 2.10635e-2L;
            c[2] = 7.54162L;

            d[0] = 0.0L;
            d[1] = 3.14499e1L;
            d[2] = 6.86346L;

            e[0] = 0.0L;
            e[1] = 1.00557L;
            e[2] = -1.28190e-1L;

        } else if (eta_order == 3) {

            D = 4.17444L;
            sigma = 9.11856e-2L;

            a[0] = 6.77740L;
            a[1] = 3.76010L;
            a[2] = 7.56690L;

            b[0] = 1.14180L;
            b[1] = 9.37188e-2L;
            b[2] = 1.16953L;

            c[0] = 2.98255L;
            c[1] = 2.10635e-2L;
            c[2] = 7.54162L;

            d[0] = 0.0L;
            d[1] = 3.05412e1L;
            d[2] = 7.88030L;

            e[0] = 0.0L;
            e[1] = 1.00557L;
            e[2] = -1.28190e-1L;

        } else {

            std::cerr << "invalid eta_order" << std::endl;
            abort();

        }
    }

    std::tuple<T, T, T>
    get_points(T eta) {
        // Return the 3 break points that define the 4 sub-intervals
        // of integration

        T term = sigma * (eta - D);
        T xi{};
        if (term > max_exponent) {
            xi = eta - D;
        } else {
            xi = std::log1p(std::exp(term)) / sigma;
        }

        T X_a = ((a[0] + b[0] * xi + c[0] * xi*xi) /
                 (1.0 + c[0] * xi));

        T X_b = ((a[1] + b[1] * xi + c[1] * d[1] * xi*xi) /
                 (1.0 + e[1] * xi + c[1] * xi*xi));

        T X_c = ((a[2] + b[2] * xi + c[2] * d[2] * xi*xi) /
                 (1.0 + e[2] * xi + c[2] * xi*xi));

        return {X_a - X_b, X_a, X_a + X_c};
    }
};


#endif
