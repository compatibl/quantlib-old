/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2014 StatPro Italia srl
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

// Based on piecewisezerospreadedtermstructure.cpp from Quantlib/test-suite.

#include "adjointzerospreadedtermstructureimpl.hpp"

using namespace QuantLib;
using namespace boost::unit_test_framework;


// Term structure of interest rates is yield curve which shows the relation between
// the interest rate (or cost of borrowing) and the time to maturity, known as the term.

// Method testFlatInterpolationLeft()
// Testing flat interpolation before the first spreaded date.
// zeroRate is interpolated in a date point by the method mentioned above using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays, yield rates) + spread quote.
// The output quantities (interpolatedZeroRate and expectedRate) are differentiated
// with respect to input yield rate vector values.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same grid.
// Flat Interpolation Left returns the output value corresponding to the node value that is immediately less than the input value.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it returns false.
bool AdjointZeroSpreadedTermStructureTest::testFlatInterpolationLeft()
{

    BOOST_MESSAGE("Testing flat interpolation before the first spreaded date...");
    TestData<FlatInterpolLeft> testData;
    size_t n = 10;
    TestData<FlatInterpolLeft>::Test test(n, &testData);

    cl::Independent(test.rates_);
    test.calculateResRates(test.rates_, test.resRates_);
    cl::tape_function<double> f(test.rates_, test.resRates_);

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

// Method testFlatInterpolationRight()
// Testing flat interpolation after the last spreaded date.
// zeroRate is interpolated in a date point by the method mentioned above using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays, yield rates) + spread quote.
// The output quantities (interpolatedZeroRate and expectedRate) are differentiated
// with respect to input yield rate vector values.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same grid.
// Flat Interpolation Right returns the output value corresponding to the node value that is immediately greater than the input value.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it returns false.
bool AdjointZeroSpreadedTermStructureTest::testFlatInterpolationRight()
{

    BOOST_MESSAGE("Testing flat interpolation after the last spreaded date...");

    TestData<FlatInterpolRight> testData;
    size_t n = 10;
    TestData<FlatInterpolRight>::Test test(n, &testData);

    cl::Independent(test.rates_);
    test.calculateResRates(test.rates_, test.resRates_);
    cl::tape_function<double> f(test.rates_, test.resRates_);

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

// Method testLinearInterpolationMultipleSpreads()
// Testing linear interpolation with more than two spreaded dates.
// zeroRate is interpolated in a date point by the method mentioned above using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays, yield rates) + spread quote.
// The output quantities (interpolatedZeroRate and expectedRate) are differentiated
// with respect to input yield rate vector values.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same grid.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it returns false.
bool AdjointZeroSpreadedTermStructureTest::testLinearInterpolationMultipleSpreads()
{

    BOOST_MESSAGE("Testing linear interpolation with more than two spreaded dates...");

    TestData<LinearInterpolMultipleSpreads> testData;
    size_t n = 10;
    TestData<LinearInterpolMultipleSpreads>::Test test(n, &testData);

    cl::Independent(test.rates_);
    test.calculateResRates(test.rates_, test.resRates_);
    cl::tape_function<double> f(test.rates_, test.resRates_);

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

// Method testLinearInterpolation()
// Testing linear interpolation between two dates.
// zeroRate is interpolated in a date point by the method mentioned above using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays, yield rates) + spread quote.
// The output quantities (interpolatedZeroRate and expectedRate) are differentiated
// with respect to input yield rate vector values.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same grid.
// Linear Interpolation fits a line between the adjacent nodes, and returns the point on that line corresponding to the input x-value.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it returns false.
bool AdjointZeroSpreadedTermStructureTest::testLinearInterpolation()
{

    BOOST_MESSAGE("Testing linear interpolation between two dates...");

    TestData<LinearInterpol> testData;
    size_t n = 10;
    TestData<LinearInterpol>::Test test(n, &testData);

    cl::Independent(test.rates_);
    test.calculateResRates(test.rates_, test.resRates_);
    cl::tape_function<double> f(test.rates_, test.resRates_);

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

// Method testForwardFlatInterpolation()
// Testing forward flat interpolation between two dates.
// zeroRate is interpolated in a date point by the method mentioned above using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays, yield rates) + spread quote.
// The output quantities (interpolatedZeroRate and expectedRate) are differentiated
// with respect to input yield rate vector values.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same grid.
// Forward Flat Interpolation returns the value of the left node.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it returns false.
bool AdjointZeroSpreadedTermStructureTest::testForwardFlatInterpolation()
{

    BOOST_MESSAGE("Testing forward flat interpolation between two dates...");

    TestData<ForwardFlatInterpol> testData;
    size_t n = 10;
    TestData<ForwardFlatInterpol>::Test test(n, &testData);

    cl::Independent(test.rates_);
    test.calculateResRates(test.rates_, test.resRates_);
    cl::tape_function<double> f(test.rates_, test.resRates_);

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

// Method testBackwardFlatInterpolation()
// Testing backward flat interpolation between two dates.
// zeroRate is interpolated in a date point by the method mentioned above using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays, yield rates) + spread quote.
// The output quantities (interpolatedZeroRate and expectedRate) are differentiated
// with respect to input yield rate vector values.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same grid.
// Backward Flat Interpolation returns the value of the right node.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it returns false.
bool AdjointZeroSpreadedTermStructureTest::testBackwardFlatInterpolation()
{

    BOOST_MESSAGE("Testing backward flat interpolation between two dates...");

    TestData<BackwardFlatInterpol> testData;
    size_t n = 10;
    TestData<BackwardFlatInterpol>::Test test(n, &testData);

    cl::Independent(test.rates_);
    test.calculateResRates(test.rates_, test.resRates_);
    cl::tape_function<double> f(test.rates_, test.resRates_);

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

// Method testDefaultInterpolation()
// Testing default interpolation between two dates.
// Linear Interpolation is used as a Default Interpolation.
// zeroRate is interpolated in a date point by the method mentioned above using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays, yield rates) + spread quote.
// The output quantities (interpolatedZeroRate and expectedRate) are differentiated
// with respect to input yield rate vector values.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same grid.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it returns false.
bool AdjointZeroSpreadedTermStructureTest::testDefaultInterpolation()
{

    BOOST_MESSAGE("Testing default interpolation between two dates...");

    TestData<DefaultInterpol> testData;
    size_t n = 10;
    TestData<DefaultInterpol>::Test test(n, &testData);

    cl::Independent(test.rates_);
    test.calculateResRates(test.rates_, test.resRates_);
    cl::tape_function<double> f(test.rates_, test.resRates_);

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

// Method testSetInterpolationFactory()
// Testing factory constructor with additional parameters.
// As a parameter for interpolation factory used Cubic Spline Interpolation -
// interpolation with special piecewise cubic polynomial.
// zeroRate is interpolated in a date point by the method mentioned above using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays, yield rates) + spread quote.
// The output quantities (interpolatedZeroRate and expectedRate) are differentiated
// with respect to input yield rate vector values.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same grid.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it returns false.
bool AdjointZeroSpreadedTermStructureTest::testSetInterpolationFactory()
{

    BOOST_MESSAGE("Testing factory constructor with additional parameters...");

    TestData<SetInterpolFactory> testData;
    size_t n = 10;
    TestData<SetInterpolFactory>::Test test(n, &testData);

    cl::Independent(test.rates_);
    test.calculateResRates(test.rates_, test.resRates_);
    cl::tape_function<double> f(test.rates_, test.resRates_);

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

// Method testQuoteChanging()
// Testing quote update.
// Backward Flat Interpolation with quote changing.
// zeroRate is interpolated in a date point by the method mentioned above using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays, yield rates) + spread quote.
// The output quantities (interpolatedZeroRate and expectedRate) are differentiated
// with respect to input yield rate vector values.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same grid.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it returns false.
bool AdjointZeroSpreadedTermStructureTest::testQuoteChanging()
{

    BOOST_MESSAGE("Testing quote update...");

    TestData<QuoteChanging> testData;
    size_t n = 10;
    TestData<QuoteChanging>::Test test(n, &testData);

    cl::Independent(test.rates_);
    test.calculateResRates(test.rates_, test.resRates_);
    cl::tape_function<double> f(test.rates_, test.resRates_);

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

    cl::Independent(test.rates_);
    test.resetSpreads(test.rates_, test.resRates_);
    cl::tape_function<double> f2(test.rates_, test.resRates_);

    // Forward mode derivatives calculation.
    test.forwardResults_.resize(n);
    dX.resize(n, 0);
    std::vector<double> dy;
    // Reverse mode derivatives calulation.
    dy = f2.Reverse(1, std::vector<double>(1, 1));

    for (size_t i = 0; i < n; i++)
    {
        if (std::fabs(test.reverseResults_[i] - dy[i]) > 1e-10)
            BOOST_ERROR("Derivatives mismatch after changing parameter spread");
    }

    return testData.makeOutput();
}

test_suite* AdjointZeroSpreadedTermStructureTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Interpolated piecewise zero spreaded yield curve tests");
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointZeroSpreadedTermStructureTest::testFlatInterpolationLeft));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointZeroSpreadedTermStructureTest::testFlatInterpolationRight));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointZeroSpreadedTermStructureTest::testLinearInterpolationMultipleSpreads));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointZeroSpreadedTermStructureTest::testLinearInterpolation));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointZeroSpreadedTermStructureTest::testBackwardFlatInterpolation));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointZeroSpreadedTermStructureTest::testForwardFlatInterpolation));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointZeroSpreadedTermStructureTest::testDefaultInterpolation));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointZeroSpreadedTermStructureTest::testSetInterpolationFactory));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointZeroSpreadedTermStructureTest::testQuoteChanging));

    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_zero_spreaded_term_structure)

BOOST_AUTO_TEST_CASE(testFlatInterpolationLeft)
{
    BOOST_CHECK(AdjointZeroSpreadedTermStructureTest::testFlatInterpolationLeft());
}

BOOST_AUTO_TEST_CASE(testFlatInterpolationRight)
{
    BOOST_CHECK(AdjointZeroSpreadedTermStructureTest::testFlatInterpolationRight());
}

BOOST_AUTO_TEST_CASE(testLinearInterpolationMultipleSpreads)
{
    BOOST_CHECK(AdjointZeroSpreadedTermStructureTest::testLinearInterpolationMultipleSpreads());
}

BOOST_AUTO_TEST_CASE(testLinearInterpolation)
{
    BOOST_CHECK(AdjointZeroSpreadedTermStructureTest::testLinearInterpolation());
}

BOOST_AUTO_TEST_CASE(testForwardFlatInterpolation)
{
    BOOST_CHECK(AdjointZeroSpreadedTermStructureTest::testForwardFlatInterpolation());
}

BOOST_AUTO_TEST_CASE(testBackwardFlatInterpolation)
{
    BOOST_CHECK(AdjointZeroSpreadedTermStructureTest::testBackwardFlatInterpolation());
}

BOOST_AUTO_TEST_CASE(testDefaultLinearInterpolation)
{
    BOOST_CHECK(AdjointZeroSpreadedTermStructureTest::testDefaultInterpolation());
}

BOOST_AUTO_TEST_CASE(testSetInterpolationFactory)
{
    BOOST_CHECK(AdjointZeroSpreadedTermStructureTest::testSetInterpolationFactory());
}

BOOST_AUTO_TEST_CASE(testBackwardFlatInterpQuoteDependence)
{
    BOOST_CHECK(AdjointZeroSpreadedTermStructureTest::testQuoteChanging());
}

BOOST_AUTO_TEST_SUITE_END()

#endif


