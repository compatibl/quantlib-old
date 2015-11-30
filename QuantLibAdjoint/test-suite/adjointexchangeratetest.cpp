/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2004 StatPro Italia srl
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

// Based on exchangerate.cpp from QuantLib \\ test-suite

#include "adjointexchangerateimpl.hpp"

using namespace QuantLib;
using namespace std;
using namespace boost::unit_test_framework;

bool AdjointExchangeRateTest::testDerived() {

    BOOST_TEST_MESSAGE("Testing derived exchange rates...");

    size_t n = 2;
    ExchangeTest test(n);

    // Tape recording.
    cl::Independent(test.rates_);
    test.calculateAmount();
    cl::tape_function<double> f(test.rates_, test.derivedAmount_);

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

    return test.check();
}

test_suite*  AdjointExchangeRateTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Exchange rate test");
    suite->add(QUANTLIB_TEST_CASE(&AdjointExchangeRateTest::testDerived));
    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_exchane_rate)

BOOST_AUTO_TEST_CASE(testDerivedExchangeRate)
{
    BOOST_CHECK(AdjointExchangeRateTest::testDerived());
}

BOOST_AUTO_TEST_SUITE_END()

#endif