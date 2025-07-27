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
#include <functional>
#include <format>
#include <print>

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

        T xsq = xi * xi;

        T X_a = ((a[0] + b[0] * xi + c[0] * xsq) /
                 (1.0 + c[0] * xi));

        T X_b = ((a[1] + b[1] * xi + c[1] * d[1] * xsq) /
                 (1.0 + e[1] * xi + c[1] * xi*xi));

        T X_c = ((a[2] + b[2] * xi + c[2] * d[2] * xsq) /
                 (1.0 + e[2] * xi + c[2] * xsq));

        return {X_a - X_b, X_a, X_a + X_c};
    }
};


template<typename T>
class FermiIntegral {

    // Construct the integral
    //
    //              ∞
    // F_k(η, β) = ∫ x^k [1 + (xβ/2)]^{1/2} / (e^{x-η} + 1) dx
    //              0
    //
    // using the method from Gong et al. 2001.  This splits the
    // integration into 4 intervals, and uses Legendre quadrature for
    // the first 3 parts and Laguerre quadrature for the last
    // interval.  For the first interval, a change of variables is
    // done to avoid singularities (effectively integrating in terms
    // of momentum instead of energy).

    // First and second derivatives with respect to $\eta$ and $\beta$
    // are supported.


    // this is roughly the maximum value of x such that std::exp(x) does
    // not overflow

    T max_exponent{};

public:

    T k{};
    T eta{};
    T beta{};

    T F{};
    T dF_deta{};
    T dF_dbeta{};
    T d2F_deta2{};
    T d2F_detadbeta{};
    T d2F_dbeta2{};

    FermiIntegral(T _k, T _eta, T _beta) :
        max_exponent{std::trunc(std::log(std::numeric_limits<T>::max()))},
        k(_k), eta(_eta), beta(_beta)
    {}

private:

    T kernel_p(T x, int eta_der, int beta_der) {

        T result{};

        // we want to work in terms of x**2
        // see Aparicio 1998 (but note they are missing a factor of 2
        // in the conversion from x to z).

        T xsq = x * x;
        T sqrt_term = std::sqrt(1.0L + 0.5L * xsq * beta);
        T num = 2.0L * std::pow(x, 2.0L * k + 1.0L) * sqrt_term;

        // this is what we are usually exponentiating
        T delta = xsq - eta;

        // this is 1.0 / (2.0 + exp(-delta) + exp(delta)
        // which is 1.0 / (2.0 * (1.0 + cosh(delta))
        T inv_cosh_term{};
        if (std::abs(delta) < max_exponent) {
            inv_cosh_term = 0.5L / (1.0L + std::cosh(delta));
        }

        // this is (exp(xsq - eta) - 1.0) / (exp(xsq - eta) + 1.0)
        T tanh_half_delta = std::tanh(0.5L * delta);

        T denomi{};
        if (delta < -max_exponent) {
            denomi = 1.0L;
        } else {
            // inv_exp_delta = np.exp(-delta)
            if (delta > max_exponent) {
                denomi = 0.0L;
            } else {
                T inv_exp_delta = std::exp(-delta);
                // 1 / (exp(x**2 - eta) + 1) rewritten
                denomi = inv_exp_delta / (1.0L + inv_exp_delta);
            }
        }

        // now construct the integrand for what we are actual computing
        if (eta_der == 0 && beta_der == 0) {
            result = num * denomi;

        } else if (eta_der == 1 && beta_der == 0) {
            // this is IB = 1 from Gong et al.
            // this corresponds to eq A.1 in terms of x**2
            result = num * inv_cosh_term;

        } else if (eta_der == 0 && beta_der == 1) {
            // this is IB = 2 from Gong et al.
            // this corresponds to eq A.2 in terms of x**2
            result = 0.5L * std::pow(x, 2.0L * k + 3.0L) / sqrt_term * denomi;

        } else if (eta_der == 2 && beta_der == 0) {
            // this is IB = 3 from Gong et al.
            // this corresponds to eq A.3 in terms of x**2
            result = num * inv_cosh_term * tanh_half_delta;

        } else if (eta_der == 1 && beta_der == 1) {
            // this is IB = 4 from Gong et al.
            // this corresponds to eq A.4 in terms of x**2
            result = 0.5L * std::pow(x, 2.0L * k + 3.0L) * inv_cosh_term / sqrt_term;

        } else if (eta_der == 0 && beta_der == 2) {
            // this is IB = 5 from Gong et al.
            // this corresponds to eq A.5 in terms of x**2
            result = -0.125L * std::pow(x, 2.0L * k + 5.0L) / std::pow(sqrt_term, 3) * denomi;
        }

        return result;

    }

    T kernel_E(T x, int eta_der, int beta_der) {

        T result{};

        // we will work in terms of x

        T sqrt_term = std::sqrt(1.0L + 0.5L * x * beta);
        T num = std::pow(x, k) * sqrt_term;

        // this is what we are usually exponentiating
        T delta = x - eta;

        // this is 1.0 / (2.0 + exp(-delta) + exp(delta)
        // which is 1.0 / (2.0 * (1.0 + cosh(delta))
        T inv_cosh_term{};
        if (std::abs(delta) < max_exponent) {
            inv_cosh_term = 0.5L / (1.0L + std::cosh(delta));
        }

        // this is (exp(x - eta) - 1.0) / (exp(x - eta) + 1.0)
        T tanh_half_delta = std::tanh(0.5L * delta);

        T denomi{};
        if (delta < -max_exponent) {
            denomi = 1.0L;
        } else {
            if (delta > max_exponent) {
                denomi = 0.0L;
            } else {
                T inv_exp_delta = std::exp(-delta);
                // 1 / (exp(x - eta) + 1) rewritten
                denomi = inv_exp_delta / (1.0L + inv_exp_delta);
            }
        }

        // now construct the integrand for what we are actual computing
        if (eta_der == 0 && beta_der == 0) {
            result = num * denomi;

        } else if (eta_der == 1 && beta_der == 0) {
            // this is IB = 1 from Gong et al.
            // this corresponds to eq A.1
            result = num * inv_cosh_term;

        } else if (eta_der == 0 && beta_der == 1) {
            // this is IB = 2 from Gong et al.
            // this corresponds to eq A.2
            result = 0.25L * std::pow(x, k + 1.0L) / sqrt_term * denomi;

        } else if (eta_der == 2 && beta_der == 0) {
            // this is IB = 3 from Gong et al.
            //  this corresponds to eq A.3
            result = num * inv_cosh_term * tanh_half_delta;

        } else if (eta_der == 1 && beta_der == 1) {
            // this is IB = 4 from Gong et al.
            // this corresponds to eq A.4
            result = 0.25L * std::pow(x, k + 1.0L) * inv_cosh_term / sqrt_term;

        } else if (eta_der == 0 && beta_der == 2) {
            // this is IB = 5 from Gong et al.
            // this corresponds to eq A.5
            result = -0.0625L * std::pow(x, k + 2.0L) / std::pow(sqrt_term, 3) * denomi;
        }

        return result;
    }

    T compute_legendre(int mode,
                       T a, T b,
                       int eta_der, int beta_der) {

        // mode = 0 means we integrate x**2
        // mode = 1 means we integrate x

        int N = static_cast<int>(x_leg.size() / 2);

        // We set the correspondence (a+b)/2 -> 0 and map the (-1,0)
        // and (0,1) intervals separately.

        T fac1 = 0.5L * (a + b);
        T fac2 = 0.5L * (b - a);

        T integral{};
        if (mode == 0) {
            for (int n = N; n < static_cast<int>(x_leg.size()); ++n) {
                integral += (kernel_p(fac1 + fac2 * x_leg[n], eta_der, beta_der) +
                             kernel_p(fac1 - fac2 * x_leg[n], eta_der, beta_der)) * w_leg[n];
            }

        } else {
            for (int n = N; n < static_cast<int>(x_leg.size()); ++n) {
                integral += (kernel_E(fac1 + fac2 * x_leg[n], eta_der, beta_der) +
                             kernel_E(fac1 - fac2 * x_leg[n], eta_der, beta_der)) * w_leg[n];
            }
        }

        integral *= 0.5L * (b - a);
        return integral;
    }

    T compute_laguerre(T a, int eta_der, int beta_der) {

        // Laguerre quadrature solves and integral of the form:
        //
        //   ∞
        // ∫  f(x) exp(-x) dx ~ ∑ f(x_i) w_i
        //  0
        //
        // where x_i are the nodes (roots of the Laguerre polynomial)
        // and w_i are the weights
        //
        // We want to integrate just f(x) from a instead of 0, then we
        // change variables x = z + a
        //
        //   ∞            ∞
        // ∫  f(x) dx = ∫  f(z + a) dz
        //  a            0
        //
        // Since Laguerre quadrature includes an exp(-x) kernel, we
        // need to scale our function by exp(x), giving
        //
        //   ∞             ∞
        // ∫  f(x) dx -> ∫  [f(x + a) exp(x)] exp(-x) dx
        //  a             0
        //
        //             ~ ∑ f(x_i) w_i exp(x) f(x + a)

        // note: the w_lag already have the exp(x) term included

        T integral{};
        for (int n = 0; n < static_cast<int>(x_lag.size()); ++n) {
            integral += kernel_E(x_lag[n] + a, eta_der, beta_der) * w_lag[n];
        }

        return integral;
    }

    T compute_fermi(int eta_der, int beta_der) {

        BreakPoints<T> bp(eta_der);
        auto [S_1, S_2, S_3] = bp.get_points(eta);

        auto I0 = compute_legendre(0, 0.0, std::sqrt(S_1),
                                   eta_der, beta_der);

        auto I1 = compute_legendre(1, S_1, S_2,
                                   eta_der, beta_der);

        auto I2 = compute_legendre(1, S_2, S_3,
                                   eta_der, beta_der);

        auto I3 = compute_laguerre(S_3, eta_der, beta_der);

        return I0 + I1 + I2 + I3;
    }

public:

    void evaluate(int mode) {
        // Perform the integration for the Fermi-Dirac function and
        // its derivatives (if desired)
        //
        // mode = 0 : F only
        // mode = 1 : F and its first derivatives
        // mode = 2 : F and its first and second derivatives


        F = compute_fermi(0, 0);

        if (mode >= 1) {
            dF_deta = compute_fermi(1, 0);
            dF_dbeta = compute_fermi(0, 1);
        }

        if (mode == 2) {
            d2F_deta2 = compute_fermi(2, 0);
            d2F_detadbeta = compute_fermi(1, 1);
            d2F_dbeta2 = compute_fermi(0, 2);
        }

    }

};

template<typename T>
std::ostream& operator<< (std::ostream& os, FermiIntegral<T>& f) {

    std::string out{};
    out += std::format("F        = {:20.10g}\n", f.F);
    out += std::format("dF/dη    = {:20.10g}\n", f.dF_deta);
    out += std::format("dF/dβ    = {:20.10g}\n", f.dF_dbeta);
    out += std::format("d²F/dη²  = {:20.10g}\n", f.d2F_deta2);
    out += std::format("d²F/dηdβ = {:20.10g}\n", f.d2F_detadbeta);
    out += std::format("d²F/dβ²  = {:20.10g}\n", f.d2F_dbeta2);

    os << out;
    return os;
}

#endif
