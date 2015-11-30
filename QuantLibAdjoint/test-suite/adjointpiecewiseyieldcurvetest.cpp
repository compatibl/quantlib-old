/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2005, 2006, 2007, 2008, 2009 StatPro Italia srl
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

//based on piecewiseyieldcurve.cpp from test-suite

#include "adjointpiecewiseyieldcurveimpl.hpp"


using namespace QuantLib;
using namespace boost::unit_test_framework;

// Test bootstrap over Japanese Yen LIBOR swaps.
bool AdjointPiecewiseYieldCurveTest::testJpyLibor()
{
    BOOST_TEST_MESSAGE("Testing bootstrap over JPY LIBOR swaps...");

    TestData testData;

    size_t n = 20;
    PiecewiceYieldCurveTest test(n, &testData);

    // Tape recording.
    cl::Independent(test.rates_);
    test.calculateNPV();
    cl::tape_function<double> f(test.rates_, test.totalNPV_);

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


test_suite* AdjointPiecewiseYieldCurveTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Piecewise yield curve tests");
    suite->add(QUANTLIB_TEST_CASE(&AdjointPiecewiseYieldCurveTest::testJpyLibor));
    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_piecewise)

BOOST_AUTO_TEST_CASE(testJpyLiborSwapNpvs)
{
    BOOST_CHECK(AdjointPiecewiseYieldCurveTest::testJpyLibor());
}

BOOST_AUTO_TEST_SUITE_END()

#endif