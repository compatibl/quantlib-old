/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
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

// Based on bonds.cpp file in testsuite

#define CL_OPEN_GENERAL_OUPUT

#include "adjointbondportfoliotest.hpp"
#include "adjointbondportfolioimpl.hpp"

using namespace QuantLib;
using namespace std;
using namespace boost::unit_test_framework;

// A bond portfolio consists of a sum of bonds with ranging yield rates: r_1,...,r_n and fixed face amount.
// The portfolio value is BondPortfolioValue= \sum_{i=1}^n e^{-r_i*t_i}.
// In this test we calculate derivatives of BondPortfolioValue with respect to r_i, i=1,...,n; n -- number of  bonds in the portfolio.
bool AdjointBondPortfolioTest::testBondPortfolio()
{
    BOOST_TEST_MESSAGE("Testing the adjoint differentiation of bond portfolio price to a change in yield forward rates...");

    TestData testData;

    size_t n = 100;
    BondPortfolioTest test(n, &testData);

    // Tape recording.
    cl::Independent(test.rate_);
    test.calculatePortfolioPrice();
    cl::tape_function<double> f(test.rate_, test.portfolioPrice_);

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

test_suite*  AdjointBondPortfolioTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("AD Bond Portfolio  test");
    suite->add(QUANTLIB_TEST_CASE(&AdjointBondPortfolioTest::testBondPortfolio));
    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_bond_portfolio)

BOOST_AUTO_TEST_CASE(testBondPortfolio)
{
    BOOST_CHECK(AdjointBondPortfolioTest::testBondPortfolio());
}

BOOST_AUTO_TEST_SUITE_END()

#endif