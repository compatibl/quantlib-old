/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2008, 2009 StatPro Italia srl
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

// Based on creditdefaultswap.cpp file from test-suite.

#include "adjointcreditdefaultswapimpl.hpp"

using namespace QuantLib;
using namespace boost::unit_test_framework;


// Adjoint derivatives of NPV and fair spread of credit default swap with respect to default probabilities are tested.
// Returns true if adjoint and finite difference derivatives match with the specified accuracy, otherwise it returns false.
bool AdjointCreditDefaultSwapTest::testDefaultProbabilities()
{
    BOOST_TEST_MESSAGE("Testing adjoint derivatives of NPV and fair spread of credit default swap"
        " with respect to default probabilities...");

    DefaultProbabilitiesTestData testData;

    size_t n = 12;

    DefaultProbabilitiesTest test(n, &testData);

    // Tape recording.
    cl::Independent(test.probability_);
    test.calculateFunction();
    cl::tape_function<double> f(test.probability_, test.calculatedFunction_);

    // Adjoint derivatives calculation.
    test.adjointResults_ = f.Jacobian(test.doubleProbability_);

    test.calcAnalytical();

    return test.checkAdjoint() && testData.makeOutput();
}


// Adjoint derivatives of NPV and fair spread of credit default swap with respect to discount factors are tested.
// Returns true if adjoint and finite difference derivatives match with the specified accuracy, otherwise it returns false.
bool AdjointCreditDefaultSwapTest::testDiscountFactor()
{
    BOOST_TEST_MESSAGE("Testing adjoint derivatives of NPV and fair spread of credit default swap"
        " with respect to discount factors...");

    DiscountFactorTestData testData;

    size_t n = 12;

    DiscountFactorTest test(n, &testData);

    // Tape recording.
    cl::Independent(test.discountFactor_);
    test.calculateFunction();
    cl::tape_function<double> f(test.discountFactor_, test.calculatedFunction_);

    // Adjoint derivatives calculation.
    test.adjointResults_ = f.Jacobian(test.doubleDiscountFactor_);

    test.calcAnalytical();

    return test.checkAdjoint() && testData.makeOutput();
}


// Adjoint derivatives of NPV of credit default swap with respect to notional, spread, and recovery rate are tested.
// Returns true if adjoint and finite difference derivatives match with the specified accuracy, otherwise it returns false.
bool AdjointCreditDefaultSwapTest::testNotionalSpreadRate()
{
    BOOST_TEST_MESSAGE("Testing adjoint derivatives of NPV of credit default swap"
        " with respect to notional, spread, and recovery rate...");

    NotionalSpreadRateTest test;

    // Tape recording.
    cl::Independent(test.parameters_);
    test.calculateFunction();
    cl::tape_function<double> f(test.parameters_, test.calculatedFunction_);

    // Adjoint derivatives calculation.
    test.adjointResults_ = f.Jacobian(test.doubleParameters_);

    test.calcAnalytical();

    return test.checkAdjoint();
}


test_suite* AdjointCreditDefaultSwapTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Adjoint credit-default swap tests");
    suite->add(QUANTLIB_TEST_CASE(&testDefaultProbabilities));
    suite->add(QUANTLIB_TEST_CASE(&testDiscountFactor));
    suite->add(QUANTLIB_TEST_CASE(&testNotionalSpreadRate));
    return suite;
}


#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_credit_default_swap)

BOOST_AUTO_TEST_CASE(testCreditDefaultSwapDefaultProbabilities)
{
    BOOST_CHECK(AdjointCreditDefaultSwapTest::testDefaultProbabilities());
}
BOOST_AUTO_TEST_CASE(testCreditDefaultSwapDiscountFactor)
{
    BOOST_CHECK(AdjointCreditDefaultSwapTest::testDiscountFactor());
}
BOOST_AUTO_TEST_CASE(testCreditDefaultSwapNotionalSpreadRate)
{
    BOOST_CHECK(AdjointCreditDefaultSwapTest::testNotionalSpreadRate());
}

BOOST_AUTO_TEST_SUITE_END()

#endif


