/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2008, 2009 StatPro Italia srl
Copyright (C) 2009 Ferdinando Ametrano
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

//based on defaultirobabilitycurve.cpp from test-suite

#include "adjointdefaultprobabilitycurveimpl.hpp"

bool AdjointDefaultProbabilityCurveTest::testFlatHazardRate()
{
    BOOST_TEST_MESSAGE("Testing flat hazard rate...");
    Size n = 10;
    TestFlatHazardRate test(n);

    // Tape recording.
    cl::Independent(test.hazardRate_);
    test.computeProbability();

    // End of tape recording.
    cl::tape_function<double> f(test.hazardRate_, test.computedProbability_);

    // Forward mode derivatives calculation.
    test.forwardResults_ = f.Forward(1, std::vector<double>(1, 1));

    // Reverse mode derivatives calulation.
    test.reverseResults_.resize(n );
    std::vector<double> dw(n, 0);
    for (size_t i = 0; i < n; i++)
    {
        dw[i] = 1;
        test.reverseResults_[i] = f.Reverse(1, dw)[0];
        dw[i] = 0;
    }

    // Calculate derivatives using analytical formula.
    test.calcAnalytical();

    return test.check();
}

bool AdjointDefaultProbabilityCurveTest::testFlatHazardRateTime()
{
    BOOST_TEST_MESSAGE("Testing flat hazard rate...");
    Size n = 10;
    TestFlatHazardRateTime test(n);

    // Tape recording.
    cl::Independent(test.t_);
    test.computeProbability();

    // End of tape recording.
    cl::tape_function<double> f(test.t_, test.computedProbability_);

    // Calculating derivatives using adjoint.
    test.adjointResults_ = f.Jacobian(test.doubleT_);

    // Calculate derivatives using analytical formula.
    test.calcAnalytical();

    return test.checkAdjoint();

}

bool AdjointDefaultProbabilityCurveTest::testFlatHazardConsistency()
{
    BOOST_TEST_MESSAGE("Testing piecewise-flat hazard-rate consistency...");
    Size n = 10;
    TestBootstrapFromSpread<HazardRate, BackwardFlat> testSpread(n);

    // Tape recording.
    cl::Independent(testSpread.quote_);
    testSpread.computeFairRate();

    // End of tape recording.
    cl::tape_function<double> f1(testSpread.quote_, testSpread.computedRate_);

    // Calculating derivatives using adjoint.
    testSpread.adjointResults_ = f1.Jacobian(testSpread.doubleQuote_);

    // Calculate derivatives using analytical formula.
    testSpread.calcAnalytical();

    bool result = testSpread.checkAdjoint();

    TestBootstrapFromUpfront<HazardRate, BackwardFlat> testUpfront(n);

    // Tape recording.
    cl::Independent(testUpfront.quote_);
    testUpfront.computeFairRate();

    // End of tape recording.
    cl::tape_function<double> f2(testUpfront.quote_, testUpfront.computedRate_);

    // Calculating derivatives using adjoint.
    testUpfront.adjointResults_ = f2.Jacobian(testUpfront.doubleQuote_);

    // Calculate derivatives using analytical formula.
    testUpfront.calcAnalytical();

    result &= testUpfront.checkAdjoint();

    return result;
}

bool AdjointDefaultProbabilityCurveTest::testFlatDensityConsistency()
{
    BOOST_TEST_MESSAGE("Testing piecewise-flat default-density consistency...");
    Size n = 10;
    TestBootstrapFromSpread<DefaultDensity, BackwardFlat> testSpread(n);

    // Tape recording.
    cl::Independent(testSpread.quote_);
    testSpread.computeFairRate();

    // End of tape recording.
    cl::tape_function<double> f1(testSpread.quote_, testSpread.computedRate_);

    // Calculating derivatives using adjoint.
    testSpread.adjointResults_ = f1.Jacobian(testSpread.doubleQuote_);

    // Calculate derivatives using analytical formula.
    testSpread.calcAnalytical();

    bool result = testSpread.checkAdjoint();

    TestBootstrapFromUpfront<DefaultDensity, BackwardFlat> testUpfront(n);

    // Tape recording.
    cl::Independent(testUpfront.quote_);
    testUpfront.computeFairRate();

    // End of tape recording.
    cl::tape_function<double> f2(testUpfront.quote_, testUpfront.computedRate_);

    // Calculating derivatives using adjoint.
    testUpfront.adjointResults_ = f2.Jacobian(testUpfront.doubleQuote_);

    // Calculate derivatives using analytical formula.
    testUpfront.calcAnalytical();

    result &= testUpfront.checkAdjoint();

    return result;
}

bool AdjointDefaultProbabilityCurveTest::testLinearDensityConsistency()
{
    BOOST_TEST_MESSAGE("Testing piecewise-linear default-density consistency...");
    Size n = 10;
    TestBootstrapFromSpread<DefaultDensity, Linear> testSpread(n);

    // Tape recording.
    cl::Independent(testSpread.quote_);
    testSpread.computeFairRate();

    // End of tape recording.
    cl::tape_function<double> f1(testSpread.quote_, testSpread.computedRate_);

    // Calculating derivatives using adjoint.
    testSpread.adjointResults_ = f1.Jacobian(testSpread.doubleQuote_);

    // Calculate derivatives using analytical formula.
    testSpread.calcAnalytical();

    bool result = testSpread.checkAdjoint();

    TestBootstrapFromUpfront<DefaultDensity, Linear> testUpfront(n);

    // Tape recording.
    cl::Independent(testUpfront.quote_);
    testUpfront.computeFairRate();

    // End of tape recording.
    cl::tape_function<double> f2(testUpfront.quote_, testUpfront.computedRate_);

    // Calculating derivatives using adjoint.
    testUpfront.adjointResults_ = f2.Jacobian(testUpfront.doubleQuote_);

    // Calculate derivatives using analytical formula.
    testUpfront.calcAnalytical();

    result &= testUpfront.checkAdjoint();

    return result;
}

bool AdjointDefaultProbabilityCurveTest::testLogLinearSurvivalConsistency()
{
    BOOST_TEST_MESSAGE("Testing log-linear survival-probability consistency...");
    Size n = 10;
    TestBootstrapFromSpread<SurvivalProbability, LogLinear> testSpread(n);

    // Tape recording.
    cl::Independent(testSpread.quote_);
    testSpread.computeFairRate();

    // End of tape recording.
    cl::tape_function<double> f1(testSpread.quote_, testSpread.computedRate_);

    // Calculating derivatives using adjoint.
    testSpread.adjointResults_ = f1.Jacobian(testSpread.doubleQuote_);

    // Calculate derivatives using analytical formula.
    testSpread.calcAnalytical();

    bool result = testSpread.checkAdjoint();

    TestBootstrapFromUpfront<SurvivalProbability, LogLinear> testUpfront(n);

    // Tape recording.
    cl::Independent(testUpfront.quote_);
    testUpfront.computeFairRate();

    // End of tape recording.
    cl::tape_function<double> f2(testUpfront.quote_, testUpfront.computedRate_);

    // Calculating derivatives using adjoint.
    testUpfront.adjointResults_ = f2.Jacobian(testUpfront.doubleQuote_);

    // Calculate derivatives using analytical formula.
    testUpfront.calcAnalytical();

    result &= testUpfront.checkAdjoint();

    return result;
}


test_suite* AdjointDefaultProbabilityCurveTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Default-probability curve tests");
    suite->add(QUANTLIB_TEST_CASE(
    &AdjointDefaultProbabilityCurveTest::testFlatHazardRate));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointDefaultProbabilityCurveTest::testFlatHazardRateTime));
   suite->add(QUANTLIB_TEST_CASE(
         &AdjointDefaultProbabilityCurveTest::testFlatHazardConsistency));
         suite->add(QUANTLIB_TEST_CASE(
         &AdjointDefaultProbabilityCurveTest::testFlatDensityConsistency));
         suite->add(QUANTLIB_TEST_CASE(
         &AdjointDefaultProbabilityCurveTest::testLinearDensityConsistency));
         suite->add(QUANTLIB_TEST_CASE(
         &AdjointDefaultProbabilityCurveTest::testLogLinearSurvivalConsistency));
    return suite;
}


#if defined CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_default_curves)

BOOST_AUTO_TEST_CASE(testFlatHazardRate)
{
    BOOST_CHECK(AdjointDefaultProbabilityCurveTest::testFlatHazardRate());
}
BOOST_AUTO_TEST_CASE(testFlatHazardRateTime)
{
    BOOST_CHECK(AdjointDefaultProbabilityCurveTest::testFlatHazardRateTime());
}
BOOST_AUTO_TEST_CASE(testFlatHazardConsistency)
{
    BOOST_CHECK(AdjointDefaultProbabilityCurveTest::testFlatHazardConsistency());
}
BOOST_AUTO_TEST_CASE(testFlatDensityConsistency)
{
    BOOST_CHECK(AdjointDefaultProbabilityCurveTest::testFlatDensityConsistency());
}
BOOST_AUTO_TEST_CASE(testLinearDensityConsistency)
{
    BOOST_CHECK(AdjointDefaultProbabilityCurveTest::testLinearDensityConsistency());
}
BOOST_AUTO_TEST_CASE(testLogLinearSurvivalConsistency)
{
    BOOST_CHECK(AdjointDefaultProbabilityCurveTest::testLogLinearSurvivalConsistency());
}
BOOST_AUTO_TEST_SUITE_END()

#endif