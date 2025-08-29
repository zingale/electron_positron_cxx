#ifndef FERMI_INTEGRALS_H
#define FERMI_INTEGRALS_H

// Compute the Fermi-Dirac integrals using the method of Gong (2001),
// https://www.sciencedirect.com/science/article/abs/pii/S001046550100145X?via%3Dihub.

// The method consist on breaking down the Fermi-integrals integration
// domain into four subintervals and apply Gauss-Legendre or the Gauss
// laguerre quadrature methods with a large number of points in each
// subdomain to guarantee accuracy beyond double-precision.

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
#include <ranges>

#include "real_type.H"
#include "mp_math.H"

#if defined(QUAD20)
#include "quadrature_weights_20.H"
#elif defined(QUAD50)
#include "quadrature_weights_50.H"
#elif defined(QUAD100)
#include "quadrature_weights_100.H"
#elif defined(QUAD400)
#include "quadrature_weights_400.H"
#else
#include "quadrature_weights_200.H"
#endif

using namespace literals;

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
        max_exponent{mp::trunc(mp::log(std::numeric_limits<T>::max()))}
    {

        if (eta_order == 0) {
            D = 3.36091_rt;
            sigma = 9.11856e-2_rt;

            a[0] = 6.77740_rt;
            a[1] = 3.76010_rt;
            a[2] = 7.56690_rt;

            b[0] = 1.14180_rt;
            b[1] = 9.37188e-2_rt;
            b[2] = 1.16953_rt;

            c[0] = 2.98255_rt;
            c[1] = 2.10635e-2_rt;
            c[2] = 7.54162e-1_rt;

            d[0] = 0.0_rt;
            d[1] = 3.10836e1_rt;
            d[2] = 6.65585_rt;

            e[0] = 0.0_rt;
            e[1] = 1.00557_rt;
            e[2] = -1.28190e-1_rt;

        } else if (eta_order == 1) {

            D = 4.99551_rt;
            sigma = 9.11856e-2_rt;

            a[0] = 6.77740_rt;
            a[1] = 3.76010_rt;
            a[2] = 7.56690_rt;

            b[0] = 1.14180_rt;
            b[1] = 9.37188e-2_rt;
            b[2] = 1.16953_rt;

            c[0] = 2.98255_rt;
            c[1] = 2.10635e-2_rt;
            c[2] = 7.54162_rt;

            d[0] = 0.0_rt;
            d[1] = 3.95015e1_rt;
            d[2] = 7.64734_rt;

            e[0] = 0.0_rt;
            e[1] = 1.00557_rt;
            e[2] = -1.28190e-1_rt;

        } else if (eta_order == 2) {

            D = 3.93830_rt;
            sigma = 9.11856e-2_rt;

            a[0] = 6.77740_rt;
            a[1] = 3.76010_rt;
            a[2] = 7.56690_rt;

            b[0] = 1.14180_rt;
            b[1] = 9.37188e-2_rt;
            b[2] = 1.16953_rt;

            c[0] = 2.98255_rt;
            c[1] = 2.10635e-2_rt;
            c[2] = 7.54162_rt;

            d[0] = 0.0_rt;
            d[1] = 3.14499e1_rt;
            d[2] = 6.86346_rt;

            e[0] = 0.0_rt;
            e[1] = 1.00557_rt;
            e[2] = -1.28190e-1_rt;

        } else if (eta_order == 3) {

            D = 4.17444_rt;
            sigma = 9.11856e-2_rt;

            a[0] = 6.77740_rt;
            a[1] = 3.76010_rt;
            a[2] = 7.56690_rt;

            b[0] = 1.14180_rt;
            b[1] = 9.37188e-2_rt;
            b[2] = 1.16953_rt;

            c[0] = 2.98255_rt;
            c[1] = 2.10635e-2_rt;
            c[2] = 7.54162_rt;

            d[0] = 0.0_rt;
            d[1] = 3.05412e1_rt;
            d[2] = 7.88030_rt;

            e[0] = 0.0_rt;
            e[1] = 1.00557_rt;
            e[2] = -1.28190e-1_rt;

        } else {

            std::cerr << "invalid eta_order" << std::endl;
            abort();

        }
    }

    auto get_points(T eta) -> std::tuple<T, T, T>
    {
        // Return the 3 break points that define the 4 sub-intervals
        // of integration

        T term = sigma * (eta - D);
        T xi{};
        if (term > max_exponent) {
            xi = eta - D;
        } else {
            xi = mp::log1p(mp::exp(term)) / sigma;
        }

        T X_a = ((a[0] + xi * (b[0] + xi * c[0])) /
                 (1.0_rt + c[0] * xi));

        T X_b = ((a[1] + xi * (b[1] + xi * c[1] * d[1])) /
                 (1.0_rt + xi * (e[1] + xi * c[1])));

        T X_c = ((a[2] + xi * (b[2] + xi * c[2] * d[2])) /
                 (1.0_rt + xi * (e[2] + xi * c[2])));

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

    // First through third derivatives with respect to $\eta$ and $\beta$
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

    T d3F_deta3{};
    T d3F_deta2dbeta{};
    T d3F_detadbeta2{};
    T d3F_dbeta3{};

    FermiIntegral(T _k, T _eta, T _beta) :
        max_exponent{mp::trunc(mp::log(std::numeric_limits<T>::max()))},
        k(_k), eta(_eta), beta(_beta)
    {}

private:

    auto kernel_p(T x, int eta_der, int beta_der) -> T
    {

        T result{};

        // we want to work in terms of x**2
        // see Aparicio 1998 (but note they are missing a factor of 2
        // in the conversion from x to z).

        T xsq = x * x;
        T sqrt_term = mp::sqrt(1.0_rt + 0.5_rt * xsq * beta);
        T num = 2.0_rt * mp::pow(x, 2.0_rt * k + 1.0_rt) * sqrt_term;

        // this is what we are usually exponentiating

        T delta = xsq - eta;

        // this is 1.0 / (2.0 + exp(-delta) + exp(delta)
        // which is 1.0 / (2.0 * (1.0 + cosh(delta))

        T inv_cosh_term{};
        if (mp::abs(delta) < max_exponent) {
            inv_cosh_term = 0.5_rt / (1.0_rt + mp::cosh(delta));
            //inv_cosh_term = 1.0_rt / (2.0_rt + std::exp(-delta) + std::exp(delta));
        }

        // this is (exp(xsq - eta) - 1.0) / (exp(xsq - eta) + 1.0)

        T tanh_half_delta = mp::tanh(0.5_rt * delta);

        T denomi{};
        if (delta < -max_exponent) {
            denomi = 1.0_rt;
        } else {
            // inv_exp_delta = np.exp(-delta)
            if (delta > max_exponent) {
                denomi = 0.0_rt;
            } else {
                //T inv_exp_delta = std::exp(-delta);
                // 1 / (exp(x**2 - eta) + 1) rewritten
                //denomi = inv_exp_delta / (1.0_rt + inv_exp_delta);
                denomi = 1.0_rt / (1.0_rt + mp::exp(delta));
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
            result = 0.5_rt * mp::pow(x, 2.0_rt * k + 3.0_rt) / sqrt_term * denomi;

        } else if (eta_der == 2 && beta_der == 0) {
            // this is IB = 3 from Gong et al.
            // this corresponds to eq A.3 in terms of x**2
            result = num * inv_cosh_term * tanh_half_delta;

        } else if (eta_der == 1 && beta_der == 1) {
            // this is IB = 4 from Gong et al.
            // this corresponds to eq A.4 in terms of x**2
            result = 0.5_rt * mp::pow(x, 2.0_rt * k + 3.0_rt) * inv_cosh_term / sqrt_term;

        } else if (eta_der == 0 && beta_der == 2) {
            // this is IB = 5 from Gong et al.
            // this corresponds to eq A.5 in terms of x**2
            result = -0.125_rt * mp::pow(x, 2.0_rt * k + 5.0_rt) / mp::pow(sqrt_term, 3) * denomi;

        } else if (eta_der == 3 && beta_der == 0) {
            // this is IB = 6 from Gong et al.
            // this corresponds to eq A.6 in terms of x**2
            result = 4.0_rt * mp::pow(x, 2.0_rt * k + 1.0_rt) * sqrt_term *
                (mp::cosh(delta) - 2.0_rt) / mp::pow(2.0_rt * (1.0_rt + mp::cosh(delta)), 2);

        } else if (eta_der == 2 && beta_der == 1) {
            // this is IB = 7 from Gong et al.
            // this corresponds to eq A.7 in terms of x**2
            result = 0.5_rt * mp::pow(x, 2.0_rt *k +3.0_rt) * inv_cosh_term *
                tanh_half_delta / sqrt_term;

        } else if (eta_der == 1 && beta_der == 2) {
            // this is IB = 8 from Gong et al.
            // this corresponds to eq A.8 in terms of x**2
            result = -0.125_rt * mp::pow(x, 2.0_rt * k +5.0_rt) * inv_cosh_term /
                mp::pow(sqrt_term, 3);

        } else if (eta_der == 0 && beta_der == 3) {
            // this is IB = 9 from Gong et al.
            // this corresponds to eq A.9 in terms of x**2
            result  = 0.03125_rt * mp::pow(x, 2.0 * k + 7.0_rt) *
                mp::pow(denomi, 3) / mp::pow(sqrt_term, 5);

        return result;

    }

    auto kernel_E(T x, int eta_der, int beta_der) -> T
    {

        T result{};

        // we will work in terms of x

        T sqrt_term = mp::sqrt(1.0_rt + 0.5_rt * x * beta);
        T num = mp::pow(x, k) * sqrt_term;

        // this is what we are usually exponentiating
        T delta = x - eta;

        // this is 1.0 / (2.0 + exp(-delta) + exp(delta)
        // which is 1.0 / (2.0 * (1.0 + cosh(delta))
        T inv_cosh_term{};
        if (mp::abs(delta) < max_exponent) {
            inv_cosh_term = 0.5_rt / (1.0_rt + mp::cosh(delta));
            //inv_cosh_term = 1.0_rt / (2.0_rt + std::exp(-delta) + std::exp(delta));
        }

        // this is (exp(x - eta) - 1.0) / (exp(x - eta) + 1.0)
        T tanh_half_delta = mp::tanh(0.5_rt * delta);

        T denomi{};
        if (delta < -max_exponent) {
            denomi = 1.0_rt;
        } else {
            if (delta > max_exponent) {
                denomi = 0.0_rt;
            } else {
                //T inv_exp_delta = std::exp(-delta);
                // 1 / (exp(x - eta) + 1) rewritten
                // denomi = inv_exp_delta / (1.0_rt + inv_exp_delta);
                denomi = 1.0_rt / (1.0_rt + mp::exp(delta));
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
            result = 0.25_rt * mp::pow(x, k + 1.0_rt) / sqrt_term * denomi;

        } else if (eta_der == 2 && beta_der == 0) {
            // this is IB = 3 from Gong et al.
            //  this corresponds to eq A.3
            result = num * inv_cosh_term * tanh_half_delta;

        } else if (eta_der == 1 && beta_der == 1) {
            // this is IB = 4 from Gong et al.
            // this corresponds to eq A.4
            result = 0.25_rt * mp::pow(x, k + 1.0_rt) * inv_cosh_term / sqrt_term;

        } else if (eta_der == 0 && beta_der == 2) {
            // this is IB = 5 from Gong et al.
            // this corresponds to eq A.5
            result = -0.0625_rt * mp::pow(x, k + 2.0_rt) / mp::pow(sqrt_term, 3) * denomi;
        }

        return result;
    }

    auto compute_legendre(int mode, T a, T b, int eta_der, int beta_der) -> T
    {

        // mode = 0 means we integrate x**2
        // mode = 1 means we integrate x

        // We set the correspondence (a+b)/2 -> 0 and map the (-1,0)
        // and (0,1) intervals separately.

        T fac1 = 0.5_rt * (a + b);
        T fac2 = 0.5_rt * (b - a);

        T integral{};
        if (mode == 0) {
            for (auto [x, w] : std::views::zip(x_leg, w_leg)) {
                integral += (kernel_p(fac1 + fac2 * x, eta_der, beta_der) +
                             kernel_p(fac1 - fac2 * x, eta_der, beta_der)) * w;
            }

        } else {
            for (auto [x, w] : std::views::zip(x_leg, w_leg)) {
                integral += (kernel_E(fac1 + fac2 * x, eta_der, beta_der) +
                             kernel_E(fac1 - fac2 * x, eta_der, beta_der)) * w;
            }
        }

        integral *= 0.5_rt * (b - a);
        return integral;
    }

    auto compute_laguerre(T a, int eta_der, int beta_der) -> T
    {

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
        for (auto [x, w] : std::views::zip(x_lag, w_lag)) {
            integral += kernel_E(x + a, eta_der, beta_der) * w;
        }

        return integral;
    }

    auto compute_fermi(int eta_der, int beta_der) -> T
    {

        // Perform the actual integration over the 4 subdomains.

        BreakPoints<T> bp(eta_der);
        auto [S_1, S_2, S_3] = bp.get_points(eta);

        auto I0 = compute_legendre(0, 0.0_rt, mp::sqrt(S_1),
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

        if (mode >= 2) {
            d2F_deta2 = compute_fermi(2, 0);
            d2F_detadbeta = compute_fermi(1, 1);
            d2F_dbeta2 = compute_fermi(0, 2);
        }

        if (mode == 3) {
            d3F_deta3 = compute_fermi(3, 0);
            d3F_deta2dbeta = compute_fermi(2, 1);
            d3F_detadbeta2 = compute_fermi(1, 2);
            d3F_dbeta3 = compute_fermi(0, 3);
        }

    }

};

template<typename T>
inline auto operator<< (std::ostream& os, FermiIntegral<T>& f) -> std::ostream&
{

    std::string out{};
    out += std::format("F         = {:20.10g}\n", f.F);
    out += std::format("∂F/∂η     = {:20.10g}\n", f.dF_deta);
    out += std::format("∂F/∂β     = {:20.10g}\n", f.dF_dbeta);

    out += std::format("∂²F/∂η²   = {:20.10g}\n", f.d2F_deta2);
    out += std::format("∂²F/∂η∂β  = {:20.10g}\n", f.d2F_detadbeta);
    out += std::format("∂²F/∂β²   = {:20.10g}\n", f.d2F_dbeta2);

    out += std::format("∂³F/∂η³   = {:20.10g}\n", f.d3F_deta3);
    out += std::format("∂³F/∂η²∂β = {:20.10g}\n", f.d3F_deta2dbeta);
    out += std::format("∂³F/∂η∂β² = {:20.10g}\n", f.d3F_detadbeta2);
    out += std::format("∂³F/∂β³   = {:20.10g}\n", f.d3F_dbeta3);

    os << out;
    return os;
}

#endif
