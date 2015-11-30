/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2005 StatPro Italia srl
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

// Based on pathgenerator.cpp from Quantlib/test-suite.

#include "adjointpathgeneratorimpl.hpp"

using namespace QuantLib;
using namespace std;
using namespace boost::unit_test_framework;


// Test multiple path generation for Black-Scholes-Merton process.
bool AdjointPathGeneratorTest::testBlackSholesPathGenerator()
{

    TestData<BlackScholes> testData;

    size_t sizeOfInd = 10;
    size_t sizeOfDep = 2 * sizeOfInd;
    TestData<BlackScholes>::Test test(sizeOfInd, &testData);

    // Tape recording.
    cl::Independent(test.sigma_);
    test.generateAssets();
    cl::tape_function<double> f(test.sigma_, test.assets_);

    // Forward mode derivatives calculation.
    test.forwardResults_.resize(sizeOfInd*sizeOfDep);
    std::vector<double> dX(sizeOfInd, 0), dy;
    for (size_t i = 0; i < sizeOfInd; i++)
    {
        dX[i] = 1;
        dy = f.Forward(1, dX);
        for (size_t j = 0; j < 2; j++)
            test.forwardResults_[(i + sizeOfInd * j)*sizeOfInd + i] = dy[i + sizeOfInd * j];
        dX[i] = 0;
    }

    // Reverse mode derivatives calulation.
    test.reverseResults_.resize(sizeOfInd*sizeOfDep);
    dX.resize(sizeOfDep, 0);
    for (size_t i = 0; i < sizeOfDep; i++)
    {
        dX[i] = 1;
        dy = f.Reverse(1, dX);
        for (size_t j = 0; j < sizeOfInd; j++)
            test.reverseResults_[i*sizeOfInd + j] = dy[j];

        dX[i] = 0;
    }

    test.calcAnalytical();

    return test.check() && testData.makeOutput();
}

bool AdjointPathGeneratorTest::testOrnsteinUhlenbeckPathGenerator()
{
    TestData<OrnstUhlenbeck> testData;

    size_t sizeOfInd = 10;
    size_t sizeOfDep = 2 * sizeOfInd;
    TestData<OrnstUhlenbeck>::Test test(sizeOfInd, &testData);

    // Tape recording.
    cl::Independent(test.sigma_);
    test.generateAssets();
    cl::tape_function<double> f(test.sigma_, test.assets_);

    // Forward mode derivatives calculation.
    test.forwardResults_.resize(sizeOfInd*sizeOfDep);
    std::vector<double> dX(sizeOfInd, 0), dy;
    for (size_t i = 0; i < sizeOfInd; i++)
    {
        dX[i] = 1;
        dy = f.Forward(1, dX);
        for (size_t j = 0; j < 2; j++)
            test.forwardResults_[(i + sizeOfInd * j)*sizeOfInd + i] = dy[i + sizeOfInd * j];
        dX[i] = 0;
    }

    // Reverse mode derivatives calulation.
    test.reverseResults_.resize(sizeOfInd*sizeOfDep);
    dX.resize(sizeOfDep, 0);
    for (size_t i = 0; i < sizeOfDep; i++)
    {
        dX[i] = 1;
        dy = f.Reverse(1, dX);
        for (size_t j = 0; j < sizeOfInd; j++)
            test.reverseResults_[i*sizeOfInd + j] = dy[j];

        dX[i] = 0;
    }

    test.calcAnalytical();

    return test.check() && testData.makeOutput();
}

bool AdjointPathGeneratorTest::testSquareRootPathGenerator()
{

    TestData<SqRoot> testData;

    size_t sizeOfInd = 10;
    size_t sizeOfDep = 2 * sizeOfInd;
    TestData<SqRoot>::Test test(sizeOfInd, &testData);

    // Tape recording.
    cl::Independent(test.sigma_);
    test.generateAssets();
    cl::tape_function<double> f(test.sigma_, test.assets_);

    // Forward mode derivatives calculation.
    test.forwardResults_.resize(sizeOfInd*sizeOfDep);
    std::vector<double> dX(sizeOfInd, 0), dy;
    for (size_t i = 0; i < sizeOfInd; i++)
    {
        dX[i] = 1;
        dy = f.Forward(1, dX);
        for (size_t j = 0; j < 2; j++)
            test.forwardResults_[(i + sizeOfInd * j)*sizeOfInd + i] = dy[i + sizeOfInd * j];
        dX[i] = 0;
    }

    // Reverse mode derivatives calulation.
    test.reverseResults_.resize(sizeOfInd*sizeOfDep);
    dX.resize(sizeOfDep, 0);
    for (size_t i = 0; i < sizeOfDep; i++)
    {
        dX[i] = 1;
        dy = f.Reverse(1, dX);
        for (size_t j = 0; j < sizeOfInd; j++)
            test.reverseResults_[i*sizeOfInd + j] = dy[j];

        dX[i] = 0;
    }

    test.calcAnalytical();

    return test.check() && testData.makeOutput();
}

test_suite* AdjointPathGeneratorTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Path generation tests");
    suite->add(QUANTLIB_TEST_CASE(&AdjointPathGeneratorTest::testBlackSholesPathGenerator));
    suite->add(QUANTLIB_TEST_CASE(&AdjointPathGeneratorTest::testOrnsteinUhlenbeckPathGenerator));
    suite->add(QUANTLIB_TEST_CASE(&AdjointPathGeneratorTest::testSquareRootPathGenerator));
    return suite;
}
#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_path_generation)

BOOST_AUTO_TEST_CASE(testBlackSholesPathGenerator)
{
    BOOST_CHECK(AdjointPathGeneratorTest::testBlackSholesPathGenerator());
}
BOOST_AUTO_TEST_CASE(testOrnsteinUhlenbeckPathGenerator)
{
    BOOST_CHECK(AdjointPathGeneratorTest::testOrnsteinUhlenbeckPathGenerator());
}
BOOST_AUTO_TEST_CASE(testSquareRootPathGenerator)
{
    BOOST_CHECK(AdjointPathGeneratorTest::testSquareRootPathGenerator());
}

BOOST_AUTO_TEST_SUITE_END()

#endif


