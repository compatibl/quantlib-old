/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2007 Marco Bianchetti
Copyright (C) 2007 Giorgio Facchinetti
Copyright (C) 2006 Chiara Fornarola
Copyright (C) 2005 StatPro Italia srl
Copyright (C) 2013 Peter Caspers
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

// Based on shortratemodels.cpp file from Quantlib/test-suite.

#include "adjointshortratemodelsimpl.hpp"

// Hull and White in 1990 proposed the following model for the short rate:
// $dr(t)=(\theta(t)-a(t)*r(t))dt +\sigma(t) dW_t$, where
// a(t) is a mean reversion and \sigma(t) is a volatility
// (a and $\sigma$ are assumed to be constants in the testing  model);
// The goal of approximating algorithm is to find the values of a and $\sigma$ in
// the  Hull-White model that best fit the market prices of the swaptions
// or find the minimum: $\sum_{i=1}^N (MarketPrice_i-PriceHW_i)^2$, where:
// N is a number of swaptions, $PriceHW= PriceHW(a,\sigma))$ is a swaption prices vector,
// calculated using Hull-White model, MarketPrice is market swaption prices vector.
// Black-Sholes implied volatilities for swaptions are provided in
// std::vector<cl::tape_double> vol_ vector of independent variables.
// Calculated using the Hull-White calibration  model $\sigma$
// (see $model\righarrow \sigma()$ in the code below) is differentiated
// with respect to input volatilities vector.
// Model $\sigma$ dependence on the first element of input volatilties vector is plotted.

bool AdjointShortRateModelsTest::testCachedHullWhite()
{
    BOOST_MESSAGE("Testing Hull-White calibration against cached values using swaptions with start delay...");

    CachedHullWhiteTestData testData(1);

    size_t n = 12;

    CachedHullWhiteTest test(n, &testData);

    // Tape recording.
    cl::Independent(test.volatility_);
    test.calculateModelSigma();
    cl::tape_function<double> f(test.volatility_, test.modelSigma_);

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

    for (Size i = n; i < n; i++)
        std::cout << "ad - " << test.adjointResults_[i] << " fd - " << test.analyticalResults_[i] << std::endl;

    bool result = test.check() && testData.makeOutput();

    Settings::instance().resetEvaluationDate();

    return result;
}

bool  AdjointShortRateModelsTest::testCachedHullWhiteFixedReversion()
{
    BOOST_MESSAGE("Testing Hull-White calibration with fixed reversion against cached values...");

    CachedHullWhiteTestData testData(2);

    size_t n = 12;

    CachedHullWhiteTest test(n, &testData);

    // Tape recording.
    cl::Independent(test.volatility_);
    test.calculateModelSigma();
    cl::tape_function<double> f(test.volatility_, test.modelSigma_);

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

    bool result = test.check() && testData.makeOutput();

    Settings::instance().resetEvaluationDate();

    return result;
}

bool AdjointShortRateModelsTest::testFuturesConvexityBias()
{
    BOOST_MESSAGE("Testing Hull-White futures convexity bias...");

    FuturesConvexityBiasTest test;

    // Tape recording.
    cl::Independent(test.parameters_);
    test.calculateFunction();
    cl::tape_function<double> f(test.parameters_, test.calculatedFunction_);

    // Forward mode derivatives calculation.
    Size sizeof_indep = test.indepVarNumber();
    test.forwardResults_.resize(sizeof_indep);
    std::vector<double> dX(sizeof_indep, 0);
    for (size_t i = 0; i < sizeof_indep; i++)
    {
        dX[i] = 1;
        test.forwardResults_[i] = f.Forward(1, dX)[0];
        dX[i] = 0;
    }

    // Reverse mode derivatives calulation.
    test.reverseResults_ = f.Reverse(1, std::vector<double>(1, 1));

    test.calcAnalytical();

    return test.check();
}

test_suite*  AdjointShortRateModelsTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("CppAD Hull-White model calibration  tests");
    suite->add(QUANTLIB_TEST_CASE(&AdjointShortRateModelsTest::testCachedHullWhite));
    suite->add(QUANTLIB_TEST_CASE(&AdjointShortRateModelsTest::testCachedHullWhiteFixedReversion));
    suite->add(QUANTLIB_TEST_CASE(&AdjointShortRateModelsTest::testFuturesConvexityBias));
    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_short_rate_models)

BOOST_AUTO_TEST_CASE(testCachedHullWhite)
{
    BOOST_CHECK(AdjointShortRateModelsTest::testCachedHullWhite());
}

BOOST_AUTO_TEST_CASE(testCachedHullWhiteFixedReversion)
{
    BOOST_CHECK(AdjointShortRateModelsTest::testCachedHullWhiteFixedReversion());
}

BOOST_AUTO_TEST_CASE(testFuturesConvexityBias)
{
    BOOST_CHECK(AdjointShortRateModelsTest::testFuturesConvexityBias());
}

BOOST_AUTO_TEST_SUITE_END()

#endif