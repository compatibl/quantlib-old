/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2014 Klaus Spanderen
Copyright (C) 2015 CompatibL

This file is part of QuantLib, a free-software/open-source library
for financial quantitative analysts and developers - http://quantlib.org/

QuantLib is free software: you can redistribute it and/or modify it
under the terms of the QuantLib license.  You should have received a
copy of the license along with this program; if not, please email
<quantlib-dev@lists.sf.net>. The license is also available online at
<http://quantlib.org/license.shtml>.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#ifndef ad_special_functions_impl_hpp
#define ad_special_functions_impl_hpp

// Based on QuantLib\ql\math\modifiedbessel.cpp file
// Implementation of Bessel functions tuned to allow adjoint differentiation
// replacing
//           std::abs(z) > a
// with
//           cl::abs_geq(z, a)
// (works while tape recording with z being a variable)

#include "adjointtestutilities.hpp"
#include <ql/quantlib.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;

#define CL_TEST_MESSAGE(msg) BOOST_TEST_MESSAGE(msg)

namespace
{
    typedef std::complex<cl::tape_double> Complex;
    typedef std::complex<double> StdComplex;

    template <class T>  struct I {};
    template <> struct I<Real> { Real value() { return 0.0; } };
    template <> struct I<Complex > {
        Complex value() { return Complex(0.0, 1.0); }
    };
    template <class T> struct Unweighted {
        T weightSmallX(const T& x) { return 1.0; }
        T weight1LargeX(const T& x) { return std::exp(x); }
        T weight2LargeX(const T& x) { return std::exp(-x); }
    };
    template <class T> struct ExponentiallyWeighted {
        T weightSmallX(const T& x) { return std::exp(-x); }
        T weight1LargeX(const T& x) { return 1.0; }
        T weight2LargeX(const T& x) { return std::exp(-2.0*x); }
    };

    template <template <class> class W>
    Complex modifiedBesselFunction_i_impl(Real nu, const Complex& x) {
        typedef Complex T;
        //if (std::abs(x) < 13.0) {
        // this original line from QuantLib replaced by the next line
        if (!cl::abs_geq(x, 13.0)) {
            const T alpha = std::pow(0.5*x, nu)
                / GammaFunction().value(1.0 + nu);
            const T Y = 0.25*x*x;
            Size k = 1;
            T sum = alpha, B_k = alpha;

            //while (std::abs(B_k *= Y / (k*(k + nu)))>std::abs(sum)*QL_EPSILON) {
            // this original line from QuantLib replaced by the next line
            while (cl::abs_geq((B_k *= Y / (k*(k + nu))) / sum, QL_EPSILON)) {
                sum += B_k;
                QL_REQUIRE(++k < 1000, "max iterations exceeded");
            }
            return sum * W<T>().weightSmallX(x);
        }
        else {
            Real na_k = 1.0, sign = 1.0;
            T da_k = T(1.0);

            T s1 = T(1.0), s2 = T(1.0);
            for (Size k = 1; k < 30; ++k) {
                sign *= -1;
                na_k *= (4.0 * nu * nu -
                    (2.0 * static_cast<Real>(k)-1.0) *
                    (2.0 * static_cast<Real>(k)-1.0));
                da_k *= (8.0 * k) * x;
                const T a_k = na_k / da_k;

                s2 += a_k;
                s1 += sign*a_k;
            }

            const T i = I<T>().value();
            return 1.0 / std::sqrt(2 * M_PI * x) *
                (W<T>().weight1LargeX(x) * s1 +
                i * std::exp(i * nu * M_PI) * W<T>().weight2LargeX(x) * s2);
        }
    }

    template <template <class> class W>
    Complex modifiedBesselFunction_k_impl(Real nu, const Complex& x) {
        return M_PI_2 * (modifiedBesselFunction_i_impl<W>(-nu, x) -
            modifiedBesselFunction_i_impl<W>(nu, x)) /
            std::sin(M_PI * nu);
    }

    Complex modifiedBesselFunction_i(Real nu, Complex &z) {
        return modifiedBesselFunction_i_impl<Unweighted>(nu, z);
    }

    Complex modifiedBesselFunction_k(Real nu, Complex &z) {
        return modifiedBesselFunction_k_impl<Unweighted>(nu, z);
    }

    Complex modifiedBesselFunction_i_exponentiallyWeighted(Real nu, Complex &z) {
        return modifiedBesselFunction_i_impl<ExponentiallyWeighted>(nu, z);
    }

    Complex modifiedBesselFunction_k_exponentiallyWeighted(Real nu, Complex &z) {
        return modifiedBesselFunction_k_impl<ExponentiallyWeighted>(nu, z);
    }

    struct nu_standart_sample
    {
        static std::vector<Real> test_sample()
        {
            std::vector<Real> vnu = {
                -1.3,
                0.001,
                1.2,
                2.3,
                -2.3,
                -10.001
            };
            return vnu;
        }
    };
}

#endif
