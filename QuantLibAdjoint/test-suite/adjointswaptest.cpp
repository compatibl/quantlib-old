/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2003, 2004, 2007 StatPro Italia srl
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

// based on swap.cpp file from test-suite

#include "adjointswapimpl.hpp"

bool AdjointSwapTest::testRateDependency()
{
    BOOST_MESSAGE("Testing vanilla-swap dependency on rate...");
    RateDependencyTestData testData;

    size_t n = 4;

    RateDependencyTest test(n, &testData);

    // Tape recording.
    cl::Independent(test.rate_);
    test.calculateTotalSwapNpv();
    cl::tape_function<double> f(test.rate_, test.swapNpv_);

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

bool AdjointSwapTest::testSpreadDependency()
{
    BOOST_MESSAGE("Testing vanilla-swap dependency on floating spread...");
    SpreadDependencyTestData testData;

    size_t n = 4;

    SpreadDependencyTest test(n, &testData);

    // Tape recording.
    cl::Independent(test.spread_);
    test.calculateTotalSwapNpv();
    cl::tape_function<double> f(test.spread_, test.swapNpv_);

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

bool AdjointSwapTest::testInArrears()
{
    BOOST_MESSAGE("Testing in-arrears swap calculation...");
    InArrearsTestData testData;

    size_t n = 10;

    InArrearsTest test(n, &testData);

    // Tape recording.
    cl::Independent(test.volatility_);
    test.calculateTotalSwapNpv();
    cl::tape_function<double> f(test.volatility_, test.swapNpv_);

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

test_suite* AdjointSwapTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Adjoint swap tests");
    suite->add(QUANTLIB_TEST_CASE(&AdjointSwapTest::testRateDependency));
    suite->add(QUANTLIB_TEST_CASE(&AdjointSwapTest::testSpreadDependency));
    suite->add(QUANTLIB_TEST_CASE(&AdjointSwapTest::testInArrears));
    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_swap)

BOOST_AUTO_TEST_CASE(testSwapDependencyOnRate)
{
    BOOST_CHECK(AdjointSwapTest::testRateDependency());
}

BOOST_AUTO_TEST_CASE(testSwapDependencyOnSpread)
{
    BOOST_CHECK(AdjointSwapTest::testSpreadDependency());
}

BOOST_AUTO_TEST_CASE(testInArrearsSwapCalc)
{
    BOOST_CHECK(AdjointSwapTest::testInArrears());
}

BOOST_AUTO_TEST_SUITE_END()

#endif


