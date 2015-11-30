/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2006, 2008 Ferdinando Ametrano
Copyright (C) 2006 François du Vignaud
Copyright (C) 2007 Cristina Duminuco
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

// Based on swaptionvolatilitymatrix.cpp from Quantlib/test - suite.

#include "adjointswaptionvolatilitymatrixtest.hpp"
#include "adjointswaptionvolatilitymatriximpl.hpp"

using namespace QuantLib;
using namespace boost::unit_test_framework;


// This test adjoint differentiation of implied volatilities on expected volatility
// of set of swaptions, which are represented in matrix form.
// Swaption volatility matrix provides the at-the-money volatility for a given
// swaption by interpolating a volatility matrix whose elements
// are the market volatilities of a set of swaptions given
// option reference date and swap lengths.
// The volatility matrix M must be defined so that:
//  -the number of rows equals the number of option dates;
//  -the number of columns equals the number of swap tenors;
//  -M[i][j] contains the volatility corresponding
// to the i-th option and j -th tenor.

// TestSwaptionVolMatrixFloatDateFloatData uses floating reference date and floating market data.
// The output implied volatilities are differentiated with respect to input expected
// volatility vector values. Adjoint differentiation results are checked using central finite
// difference method.
bool AdjointSwaptionVolatilityMatrixTest::testSwaptionVolMatrixFloatDateFloatData()
{
    BOOST_TEST_MESSAGE("Testing swaption volatility matrix...");

    size_t optionSize = 5;
    size_t n = optionSize*maxSwapSize;
    TestData testData(DataType::floatDateFloatData);

    SwaptionMatrixTest test(optionSize, &testData);

    // Tape recording.
    cl::Independent(test.expVol_);
    test.calculateImplVol();
    cl::tape_function<double> f(test.expVol_, test.implVol_);

    // Forward mode derivatives calculation.
    test.forwardResults_.resize(n*n);
    std::vector<double> dX(n, 0), dy(n);
    for (size_t i = 0; i < n; i++)
    {
        dX[i] = 1;
        dy = f.Forward(1, dX);
        dX[i] = 0;
        for (size_t j = 0; j < n; j++)
            test.forwardResults_[i*n + j] = dy[j];
    }

    // Reverse mode derivatives calulation.
    test.reverseResults_.resize(n*n);
    std::vector<double> dw(n, 0);

    for (size_t i = 0; i < n; i++)
    {
        dw[i] = 1;
        dy = f.Reverse(1, dw);
        dw[i] = 0;
        for (size_t j = 0; j < n; j++)
            test.reverseResults_[i*n + j] = dy[j];
    }

    test.calcAnalytical();

    return test.check() && testData.makeOutput();

}

// TestSwaptionVolMatrixFixedDateFloatData uses fixed reference date and floating market data.
// The output implied volatilities are differentiated with respect to input expected
// volatility vector values. Adjoint differentiation results are checked using central finite
// difference method.
bool AdjointSwaptionVolatilityMatrixTest::testSwaptionVolMatrixFixedDateFloatData()
{
    BOOST_TEST_MESSAGE("Testing swaption volatility matrix...");

    size_t optionSize = 5;
    size_t n = optionSize*maxSwapSize;
    TestData testData(DataType::fixedDateFloatData);

    SwaptionMatrixTest test(optionSize, &testData);

    // Tape recording.
    cl::Independent(test.expVol_);
    test.calculateImplVol();
    cl::tape_function<double> f(test.expVol_, test.implVol_);

    // Forward mode derivatives calculation.
    test.forwardResults_.resize(n*n);
    std::vector<double> dX(n, 0), dy(n);
    for (size_t i = 0; i < n; i++)
    {
        dX[i] = 1;
        dy = f.Forward(1, dX);
        dX[i] = 0;
        for (size_t j = 0; j < n; j++)
            test.forwardResults_[i*n + j] = dy[j];
    }

    // Reverse mode derivatives calulation.
    test.reverseResults_.resize(n*n);
    std::vector<double> dw(n, 0);

    for (size_t i = 0; i < n; i++)
    {
        dw[i] = 1;
        dy = f.Reverse(1, dw);
        dw[i] = 0;
        for (size_t j = 0; j < n; j++)
            test.reverseResults_[i*n + j] = dy[j];
    }
    test.calcAnalytical();

    return test.check() && testData.makeOutput();

}

// TestSwaptionVolMatrixFloatDateFixedData uses floating reference date and fixed market data.
// The output implied volatilities are differentiated with respect to input expected
// volatility vector values. Adjoint differentiation results are checked using central finite
// difference method.
bool AdjointSwaptionVolatilityMatrixTest::testSwaptionVolMatrixFloatDateFixedData()
{
    BOOST_TEST_MESSAGE("Testing swaption volatility matrix...");

    size_t optionSize = 5;
    size_t n = optionSize*maxSwapSize;
    TestData testData(DataType::floatDateFixedData);

    SwaptionMatrixTest test(optionSize, &testData);

    // Tape recording.
    cl::Independent(test.expVol_);
    test.calculateImplVol();
    cl::tape_function<double> f(test.expVol_, test.implVol_);

    // Forward mode derivatives calculation.
    test.forwardResults_.resize(n*n);
    std::vector<double> dX(n, 0), dy(n);
    for (size_t i = 0; i < n; i++)
    {
        dX[i] = 1;
        dy = f.Forward(1, dX);
        dX[i] = 0;
        for (size_t j = 0; j < n; j++)
            test.forwardResults_[i*n + j] = dy[j];
    }

    // Reverse mode derivatives calulation.
    test.reverseResults_.resize(n*n);
    std::vector<double> dw(n, 0);

    for (size_t i = 0; i < n; i++)
    {
        dw[i] = 1;
        dy = f.Reverse(1, dw);
        dw[i] = 0;
        for (size_t j = 0; j < n; j++)
            test.reverseResults_[i*n + j] = dy[j];
    }
    test.calcAnalytical();

    return test.check() && testData.makeOutput();

}

// TestSwaptionVolMatrixFixedDateFloatData uses fixed reference date and fixed market data.
// The output implied volatilities are differentiated with respect to input expected
// volatility vector values. Adjoint differentiation results are checked using central finite
// difference method.
bool AdjointSwaptionVolatilityMatrixTest::testSwaptionVolMatrixFixedDateFixedData()
{
    BOOST_TEST_MESSAGE("Testing swaption volatility matrix...");

    size_t optionSize = 5;
    size_t n = optionSize*maxSwapSize;
    TestData testData(DataType::fixedDateFixedData);

    SwaptionMatrixTest test(optionSize, &testData);

    // Tape recording.
    cl::Independent(test.expVol_);
    test.calculateImplVol();
    cl::tape_function<double> f(test.expVol_, test.implVol_);

    // Forward mode derivatives calculation.
    test.forwardResults_.resize(n*n);
    std::vector<double> dX(n, 0), dy(n);
    for (size_t i = 0; i < n; i++)
    {
        dX[i] = 1;
        dy = f.Forward(1, dX);
        dX[i] = 0;
        for (size_t j = 0; j < n; j++)
            test.forwardResults_[i*n + j] = dy[j];
    }

    // Reverse mode derivatives calulation.
    test.reverseResults_.resize(n*n);
    std::vector<double> dw(n, 0);

    for (size_t i = 0; i < n; i++)
    {
        dw[i] = 1;
        dy = f.Reverse(1, dw);
        dw[i] = 0;
        for (size_t j = 0; j < n; j++)
            test.reverseResults_[i*n + j] = dy[j];
    }

    test.calcAnalytical();

    return test.check() && testData.makeOutput();
}


test_suite* AdjointSwaptionVolatilityMatrixTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Swaption Volatility Matrix tests");

   suite->add(QUANTLIB_TEST_CASE(
        &AdjointSwaptionVolatilityMatrixTest::testSwaptionVolMatrixFloatDateFloatData));
     suite->add(QUANTLIB_TEST_CASE(
        &AdjointSwaptionVolatilityMatrixTest::testSwaptionVolMatrixFixedDateFloatData));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointSwaptionVolatilityMatrixTest::testSwaptionVolMatrixFloatDateFixedData));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointSwaptionVolatilityMatrixTest::testSwaptionVolMatrixFixedDateFixedData));

    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_swaption_matrix)

BOOST_AUTO_TEST_CASE(testSwaptionMatrixFloatDateFloatData)
{
    BOOST_CHECK(AdjointSwaptionVolatilityMatrixTest::testSwaptionVolMatrixFloatDateFloatData());
}
BOOST_AUTO_TEST_CASE(testSwaptionMatrixFixedDateFloatData)
{
    BOOST_CHECK(AdjointSwaptionVolatilityMatrixTest::testSwaptionVolMatrixFixedDateFloatData());
}
BOOST_AUTO_TEST_CASE(testSwaptionMatrixFloatDateFixedData)
{
    BOOST_CHECK(AdjointSwaptionVolatilityMatrixTest::testSwaptionVolMatrixFloatDateFixedData());
}
BOOST_AUTO_TEST_CASE(testSwaptionMatrixFixedDateFixedData)
{
    BOOST_CHECK(AdjointSwaptionVolatilityMatrixTest::testSwaptionVolMatrixFixedDateFixedData());
}

BOOST_AUTO_TEST_SUITE_END()

#endif