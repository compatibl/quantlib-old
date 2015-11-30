/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2006 Allen Kuo
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

#include "adjointfraportfolioimpl.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include <ql/quantlib.hpp>

using namespace QuantLib;
using namespace std;
using namespace boost::unit_test_framework;


// This example tests adjoint differentiation of forward value
// of Forward Rate Agreement portfolio on initial forward rates.
// Method returns true if adjoint derivatives equals finite difference derivative; otherwise, it throws exception and returns false.

bool AdjointFRAPortfolioTest::testFRAPortfolio()
{
    BOOST_TEST_MESSAGE("Testing the adjoint differentiation of contract value of FRA portfolio to a change in initial forward rates...");

    TestData testData;

    size_t n = 20;
    FRAPortfolioTest test(n, &testData);

    // Tape recording.
    cl::Independent(test.fraRate_);
    test.calculateTotalContractValue();
    cl::tape_function<double> f(test.fraRate_, test.totalContractValue_);

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

test_suite* AdjointFRAPortfolioTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("AD FRA Portfolio test");
    suite->add(QUANTLIB_TEST_CASE(&AdjointFRAPortfolioTest::testFRAPortfolio));
    return suite;
}


#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_fra_portfolio)

BOOST_AUTO_TEST_CASE(testFRAPortfolio)
{
    BOOST_CHECK(AdjointFRAPortfolioTest::testFRAPortfolio());
}

BOOST_AUTO_TEST_SUITE_END()

#endif
