/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2003 RiskMap srl
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

// Based on termstructure.cpp from Quantlib/test-suite.

#include "ql\quantlib.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtermstructureimpl.hpp"

using namespace QuantLib;
using namespace boost::unit_test_framework;

// Term structure of interest rates is yield curve which shows the relation between
// the interest rate (or cost of borrowing) and the time to maturity, known as the term.

// Method testZSpreaded()
// Testing consistency of zero-spreaded term structure.
// Zero rate is interpolated in a date point by the log-linear method using the following formula:
// zeroRate = Interpolation(testDate, yield rates).
// The zero rate is differentiated
// with respect to input swap rate vector.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same grid.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it returns false.
bool AdjointTermStructureTest::testZSpreaded()
{
    BOOST_TEST_MESSAGE("Testing consistency of zero-spreaded term structure...");

    ZSpreadedTestData testData;
    size_t n = 30;
    ZSpreadedTest test(n, &testData);

    // Tape recording.
    cl::Independent(test.rate_);
    test.calculateZeroValue();
    cl::tape_function<double> f(test.rate_, test.zeroValue_);

    // Forward mode derivatives calculation.
    test.forwardResults_.resize(n);
    std::vector<double> dX(n, 0);
    for (size_t i = 0; i < n; i++)
    {
        dX[i] = 1;
        test.forwardResults_[i] = f.Forward(1, dX)[0];
        dX[i] = 0;
    }

    // Reverse mode derivatives calulation.
    test.reverseResults_ = f.Reverse(1, std::vector<double>(1, 1));

    test.calcAnalytical();
    return test.check() && testData.makeOutput();
}

// Method testFSpreaded()
// Testing consistency of forward-spreaded term structure.
// The forward rate is the future yield of a bond. It is calculated using the yield curve.
// r_{t_1,t_2} = 1 / (d_2-d_1) * ((1 + r_2 * d_2) / (1 + r_1 * d_1) - 1)
// r_{t_1,t_2} is the forward rate between term t_1 and term t_2,
// d_1 is the time length between time 0 and term t_1 (in years),
// d_2 is the time length between time 0 and term t_2 (in years),
// r_1 is the zero-coupon yield for the time period (0, t_1),
// r_2 is the zero-coupon yield for the time period (0, t_2).
// The forward rate is differentiated
// with respect to input swap rate vectors.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same grid.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it returns false.
bool AdjointTermStructureTest::testFSpreaded()
{
    BOOST_TEST_MESSAGE("Testing consistency of forward-spreaded term structure...");

    FSpreadedTestData testData;
    size_t n = 30;
    FSpreadedTest test(n, &testData);

    // Tape recording.
    cl::Independent(test.rate_);
    test.calculateForwardValue();
    cl::tape_function<double> f(test.rate_, test.forwardValue_);

    // Forward mode derivatives calculation.
    test.forwardResults_.resize(n);
    std::vector<double> dX(n, 0);
    for (size_t i = 0; i < n; i++)
    {
        dX[i] = 1;
        test.forwardResults_[i] = f.Forward(1, dX)[0];
        dX[i] = 0;
    }

    // Reverse mode derivatives calulation.
    test.reverseResults_ = f.Reverse(1, std::vector<double>(1, 1));

    test.calcAnalytical();
    return test.check() && testData.makeOutput();
}

// Method testImplied()
// Testing consistency of yield term structure.
// Discounts are interpolated in a date point by the log-linear method using the following formula:
// discount = exp( - InterpolatedInterestRate(testDate, newSettlementDate) * time).
// The output quantities (base discount, discount and implied discount) are differentiated
// with respect to input swap rate vectors.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same grid.
// Method returns true if adjoint derivatives are close enough to finite difference derivatives; otherwise, otherwise,it returns false.
bool AdjointTermStructureTest::testImplied()
{
    BOOST_TEST_MESSAGE("Testing consistency of implied term structure...");

    ImpliedTestData testData;
    size_t n = 30;
    ImpliedTest test(n, &testData);

    // Tape recording.
    cl::Independent(test.rate_);
    test.calculateDiscountValues();
    cl::tape_function<double> f(test.rate_, test.discountValues_);

    // Adjoint derivatives calculation.
    test.adjointResults_ = f.Jacobian(test.doubleRate_);

    test.calcAnalytical();
    return test.checkAdjoint() && testData.makeOutput();
}

test_suite* AdjointTermStructureTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Term structure tests");
    suite->add(QUANTLIB_TEST_CASE(&AdjointTermStructureTest::testZSpreaded));
    suite->add(QUANTLIB_TEST_CASE(&AdjointTermStructureTest::testFSpreaded));
    suite->add(QUANTLIB_TEST_CASE(&AdjointTermStructureTest::testImplied));

    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_termstructure)

BOOST_AUTO_TEST_CASE(testConsistencyZeroSpreadedTermStructure)
{
    BOOST_CHECK(AdjointTermStructureTest::testZSpreaded());
}

BOOST_AUTO_TEST_CASE(testConsistencyForwardSpreadedTermStructure)
{
    BOOST_CHECK(AdjointTermStructureTest::testFSpreaded());
}

BOOST_AUTO_TEST_CASE(testConsistencyYieldTermStructure)
{
    BOOST_CHECK(AdjointTermStructureTest::testImplied());
}
BOOST_AUTO_TEST_SUITE_END()

#endif






