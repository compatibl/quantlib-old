/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2002, 2003 Ferdinando Ametrano
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2010 Kakhkhor Abdijalilov
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

/*! \file normaldistribution.hpp
    \brief normal, cumulative and inverse cumulative distributions
*/

#ifndef quantlib_normal_distribution_t_hpp
#define quantlib_normal_distribution_t_hpp

#include <cl/tape/impl/ql/errorfunction_t.hpp>
#include <ql/errors.hpp>

namespace QuantLib {

    //! Normal distribution function
    /*! Given x, it returns its probability in a Gaussian normal distribution.
        It provides the first derivative too.

        \test the correctness of the returned value is tested by
              checking it against numerical calculations. Cross-checks
              are also performed against the
              CumulativeNormalDistribution and InverseCumulativeNormal
              classes.
    */
    template <class T>
    class NormalDistribution_T : public std::unary_function<T, T>
    {
    public:
        NormalDistribution_T(T const& average = T(0.0),
            T const& sigma = T(1.0));
        // function
        T operator()(T const& x) const;
        T derivative(T const& x) const;
    private:
        T average_, sigma_, normalizationFactor_, denominator_,
            derNormalizationFactor_;
    };

    template <class T>
    using GaussianDistribution_T = NormalDistribution_T<T>;


    //! Cumulative normal distribution function
    /*! Given x it provides an approximation to the
        integral of the gaussian normal distribution:
        formula here ...

        For this implementation see M. Abramowitz and I. Stegun,
        Handbook of Mathematical Functions,
        Dover Publications, New York (1972)
    */
    template <class T>
    class CumulativeNormalDistribution_T
        : public std::unary_function<T, T>
    {
    public:
        CumulativeNormalDistribution_T(T const& average = T(0.0),
            T const& sigma = T(1.0));
        // function
        T operator()(T const& x) const;
        T derivative(T const& x) const;
    private:
        T average_, sigma_;
        NormalDistribution_T<T> gaussian_;
        ErrorFunction_T<T> errorFunction_;
    };


    template <class T>
    inline NormalDistribution_T<T>::NormalDistribution_T(T const& average,
                                                  T const& sigma)
        : average_(average), sigma_(sigma)
    {

        QL_REQUIRE(sigma_>0.0,
                   "sigma must be greater than 0.0 ("
                   << sigma_ << " not allowed)");

        normalizationFactor_ = M_SQRT_2*M_1_SQRTPI/sigma_;
        derNormalizationFactor_ = sigma_*sigma_;
        denominator_ = 2.0*derNormalizationFactor_;
    }

    template <class T>
    inline T NormalDistribution_T<T>::operator()(T const& x) const
    {
        T deltax = x-average_;
        T exponent = -(deltax*deltax)/denominator_;
        // debian alpha had some strange problem in the very-low range
        return exponent <= -690.0 ? T(0.0) :  // exp(x) < 1.0e-300 anyway
            normalizationFactor_*std::exp(exponent);
    }

    template <class T>
    inline T NormalDistribution_T<T>::derivative(T const& x) const
    {
        return ((*this)(x) * (average_ - x)) / derNormalizationFactor_;
    }

    template <class T>
    inline CumulativeNormalDistribution_T<T>::CumulativeNormalDistribution_T(
        T const& average, T const& sigma)
        : average_(average), sigma_(sigma)
    {

        QL_REQUIRE(sigma_>0.0,
                   "sigma must be greater than 0.0 ("
                   << sigma_ << " not allowed)");
    }

    template <class T>
    inline T CumulativeNormalDistribution_T<T>::derivative(T const& x) const
    {
        T xn = (x - average_) / sigma_;
        return gaussian_(xn) / sigma_;
    }

    template <class T>
    T CumulativeNormalDistribution_T<T>::operator()(T const& z1) const {
        //QL_REQUIRE(!(z >= average_ && 2.0*average_-z > average_),
        //           "not a real number. ");
        T z = (z1 - average_) / sigma_;

        T result = 0.5 * (1.0 + errorFunction_(z*M_SQRT_2));
        if (result <= 1e-8) { //todo: investigate the threshold level
            // Asymptotic expansion for very negative z following (26.2.12)
            // on page 408 in M. Abramowitz and A. Stegun,
            // Pocketbook of Mathematical Functions, ISBN 3-87144818-4.
            T sum = 1.0, zsqr = z*z, i = 1.0, g = 1.0, x, y,
                a = QL_MAX_REAL, lasta;
            do {
                lasta = a;
                x = (4.0*i - 3.0) / zsqr;
                y = x*((4.0*i - 1) / zsqr);
                a = g*(x - y);
                sum -= a;
                g *= y;
                ++i;
                a = std::fabs(a);
            } while (lasta>a && a >= std::fabs(sum*QL_EPSILON));
            result = -gaussian_(z) / z*sum;
        }
        return result;
    }

}


#endif
