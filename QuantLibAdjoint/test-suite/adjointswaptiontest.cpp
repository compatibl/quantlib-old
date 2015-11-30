/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2003 RiskMap srl
Copyright (C) 2006, 2007 Ferdinando Ametrano
Copyright (C) 2006 Marco Bianchetti
Copyright (C) 2006 Cristina Duminuco
Copyright (C) 2007, 2008 StatPro Italia srl
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

#include <ql/quantlib.hpp>
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include "adjointswaptionimpl.hpp"

// based on swaption.cpp file from test-suite
using namespace QuantLib;
using namespace boost::unit_test_framework;


// Test swaption dependency on spread.
bool AdjointSwaptionTest::testSpreadDependency()
{
    BOOST_MESSAGE("Testing swaption dependency on spread...");

    SpreadDependencyTestData testData;

    size_t n = 2;
    SpreadDependencyTest test(n, &testData);

    // Tape recording.
    cl::Independent(test.spread_);
    test.calculateTotalSwaptionNpv();
    cl::tape_function<double> f(test.spread_, test.swaptionNpv_);

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

// Test swaption correction of spread.
bool AdjointSwaptionTest::testSpreadCorrection()
{
    BOOST_MESSAGE("Testing swaption treatment of spread...");

    SpreadCorrectionTestData testData;

    size_t n = 2;
    SpreadCorrectionTest test(n, &testData);

    // Tape recording.
    cl::Independent(test.spread_);
    test.calculateTotalSwaptionNpv();
    cl::tape_function<double> f(test.spread_, test.swaptionNpv_);

    // Adjoint derivatives calculation.
    test.adjointResults_ = f.Jacobian(test.doubleSpread_);

    test.calcAnalytical();

    return test.checkAdjoint() && testData.makeOutput();
}

// Test swaption value against cached value.
bool AdjointSwaptionTest::testCachedValue()
{
    BOOST_MESSAGE("Testing swaption value against cached value...");

    CachedValueTestData testData;

    size_t n = 2;
    CachedValueTest test(n, &testData);

    // Tape recording.
    cl::Independent(test.volatility_);
    test.calculateTotalSwaptionNpv();
    cl::tape_function<double> f(test.volatility_, test.swaptionNpv_);

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


test_suite* AdjointSwaptionTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Adjoint swaption test");
    suite->add(QUANTLIB_TEST_CASE(&AdjointSwaptionTest::testSpreadDependency));
    suite->add(QUANTLIB_TEST_CASE(&AdjointSwaptionTest::testSpreadCorrection));
    suite->add(QUANTLIB_TEST_CASE(&AdjointSwaptionTest::testCachedValue));
    return suite;
}


#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_swaption)

BOOST_AUTO_TEST_CASE(testSwaptDependencyOnSpread)
{
    BOOST_CHECK(AdjointSwaptionTest::testSpreadDependency());
}

BOOST_AUTO_TEST_CASE(testSwaptCorrectionOfSpread)
{
    BOOST_CHECK(AdjointSwaptionTest::testSpreadCorrection());
}

BOOST_AUTO_TEST_CASE(testSwaptCachedValue)
{
    BOOST_CHECK(AdjointSwaptionTest::testCachedValue());
}

BOOST_AUTO_TEST_SUITE_END()

#endif
