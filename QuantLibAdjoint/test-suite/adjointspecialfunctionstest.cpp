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

#include "adjointspecialfunctionstest.hpp"
#include "adjointspecialfunctionsimpl.hpp"
#include "adjointcomplexutils.hpp"
#include <ql/quantlib.hpp>


bool AdjointSpecialFunctionsTest::testModifiedBesselI()
{
    CL_TEST_MESSAGE("Testing modified Bessel Function of 1st kind...");

    bool ok = true;
    double tol = 1e-10;
    std::vector<Real> vnu = nu_standart_sample::test_sample();

    for (Real nu : vnu) {
        auto f = [nu](Complex z)
        {
            // modified Bessel function of 1st kind
            // (http://functions.wolfram.com/Bessel-TypeFunctions/BesselI/)
            return modifiedBesselFunction_i(nu, z);
        };
        auto df = [nu](Complex z)
        {
            // analytic derivative (http://functions.wolfram.com/Bessel-TypeFunctions/BesselI/)
            // dI_nu(z)/dz =  I_(nu-1)(z) - nu/z I_nu(z)
            return modifiedBesselFunction_i(nu - 1, z) - nu / z * modifiedBesselFunction_i(nu, z);
        };

        ok = ok && testFunctionDerivative<Real, small_arguments>(f, df, tol)
            && testFunctionDerivative<Complex, small_arguments>(f, df, tol);
    }
    return ok;
}

bool AdjointSpecialFunctionsTest::testModifiedBesselK()
{
    CL_TEST_MESSAGE("Testing modified Bessel Function of 2nd kind...");

    bool ok = true;
    double tol = 1e-8;
    std::vector<Real> vnu = nu_standart_sample::test_sample();

    for (Real nu : vnu) {
        auto f = [nu](Complex z)
        {
            // modified Bessel function of 2nd kind
            // (http://functions.wolfram.com/Bessel-TypeFunctions/BesselK/)
            return modifiedBesselFunction_k(nu, z);
        };
        auto df = [nu](Complex z)
        {
            // analytic derivative (http://functions.wolfram.com/Bessel-TypeFunctions/BesselK/)
            // dK_nu(z)/dz = -1/2 ( K_(nu-1)(z) + K_(nu+1)(z) )
            return -1. / 2. * (modifiedBesselFunction_k(nu - 1, z) + modifiedBesselFunction_k(nu + 1, z));
        };

        ok = ok && testFunctionDerivative<Real, small_arguments>(f, df, tol)
            && testFunctionDerivative<Complex, small_arguments>(f, df, tol);
    }
    return ok;
}

bool AdjointSpecialFunctionsTest::testModifiedBesselIWeighted()
{
    CL_TEST_MESSAGE("Testing exponentially weighted modified Bessel Function of 1st kind...");

    bool ok = true;
    double tol = 1e-8;
    std::vector<Real> vnu = nu_standart_sample::test_sample();

    for (Real nu : vnu) {
        auto f = [nu](Complex z)
        {
            // modified Bessel function of 2nd kind
            // (http://functions.wolfram.com/Bessel-TypeFunctions/BesselK/)
            // exponentially weighted: Iexp_nu(z) = I_nu(z) * exp(-z)
            return modifiedBesselFunction_i_exponentiallyWeighted(nu, z);
        };
        auto df = [nu](Complex z)
        {
            // analytic derivative (http://functions.wolfram.com/Bessel-TypeFunctions/BesselI/)
            // dIexp_nu(z)/dz =  Iexp_(nu-1)(z) - nu/z Iexp_nu(z) - Iexp_nu(z)
            return modifiedBesselFunction_i_exponentiallyWeighted(nu - 1, z)
                - nu / z * modifiedBesselFunction_i_exponentiallyWeighted(nu, z)
                - modifiedBesselFunction_i_exponentiallyWeighted(nu, z);
        };

        ok = ok && testFunctionDerivative<Real, small_arguments>(f, df, tol)
            && testFunctionDerivative<Complex, small_arguments>(f, df, tol);
    }
    return ok;
}

bool AdjointSpecialFunctionsTest::testModifiedBesselKWeighted()
{
    CL_TEST_MESSAGE("Testing exponentially weighted modified Bessel Function of 2nd kind...");

    bool ok = true;
    double tol = 1e-8;
    std::vector<Real> vnu = nu_standart_sample::test_sample();

    for (Real nu : vnu) {
        auto f = [nu](Complex z)
        {
            // modified Bessel function of 2nd kind
            // (http://functions.wolfram.com/Bessel-TypeFunctions/BesselK/)
            // exponentially weighted: Kexp_nu(z) = K_nu(z) * exp(-z)
            return modifiedBesselFunction_k_exponentiallyWeighted(nu, z);
        };
        auto df = [nu](Complex z)
        {
            // analytic derivative (http://functions.wolfram.com/Bessel-TypeFunctions/BesselK/)
            // dKexp_nu(z)/dz = -1/2 ( Kexp_(nu-1)(z) + Kexp_(nu+1)(z) ) - Kexp_nu(z)
            return -1. / 2. * (modifiedBesselFunction_k_exponentiallyWeighted(nu - 1, z)
                + modifiedBesselFunction_k_exponentiallyWeighted(nu + 1, z))
                - modifiedBesselFunction_k_exponentiallyWeighted(nu, z);
        };

        ok = ok && testFunctionDerivative<Real, small_arguments>(f, df, tol)
            && testFunctionDerivative<Complex, small_arguments>(f, df, tol);
    }
    return ok;
}

test_suite* AdjointSpecialFunctionsTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Adjoint Special Functions Differentiantion test");
    suite->add(QUANTLIB_TEST_CASE(&AdjointSpecialFunctionsTest::testModifiedBesselI));
    suite->add(QUANTLIB_TEST_CASE(&AdjointSpecialFunctionsTest::testModifiedBesselK));
    suite->add(QUANTLIB_TEST_CASE(&AdjointSpecialFunctionsTest::testModifiedBesselIWeighted));
    suite->add(QUANTLIB_TEST_CASE(&AdjointSpecialFunctionsTest::testModifiedBesselKWeighted));

    return suite;
}


#if defined CL_TAPE_COMPLEX_ENABLED

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_special_functions)

BOOST_AUTO_TEST_CASE(testModifiedBesselI)
{
    BOOST_CHECK(AdjointSpecialFunctionsTest::testModifiedBesselI());
}
BOOST_AUTO_TEST_CASE(testModifiedBesselK)
{
    BOOST_CHECK(AdjointSpecialFunctionsTest::testModifiedBesselK());
}
BOOST_AUTO_TEST_CASE(testModifiedBesselIWeighted)
{
    BOOST_CHECK(AdjointSpecialFunctionsTest::testModifiedBesselIWeighted());
}
BOOST_AUTO_TEST_CASE(testModifiedBesselKWeighted)
{
    BOOST_CHECK(AdjointSpecialFunctionsTest::testModifiedBesselKWeighted());
}

BOOST_AUTO_TEST_SUITE_END()

#endif

#endif
