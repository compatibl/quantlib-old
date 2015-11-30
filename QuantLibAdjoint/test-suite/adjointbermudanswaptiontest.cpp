/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2005, 2007 StatPro Italia srl
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


#include "adjointbermudanswaptionimpl.hpp"


bool AdjointBermudanSwaptionTest::testTreeSwaptionEngine()
{

    BOOST_TEST_MESSAGE("Testing spread sensitivity of Bermudan swaption (with TreeSwaptionEngine) ...");

    TestData testData(tree);

    size_t n = 3;
    SwaptionTest test(n, &testData);

    // Tape recording.
    cl::Independent(test.spread_);
    test.calculatePrices();
    cl::tape_function<double> f(test.spread_, test.totalPrice_);

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


bool AdjointBermudanSwaptionTest::testFdHullWhiteSwaptionEngine()
{

    BOOST_TEST_MESSAGE("Testing spread sensitivity of Bermudan swaption (with FdHullWhiteSwaptionEngine) ...");

    TestData testData(fdHullWhite);

    size_t n = 3;
    SwaptionTest test(n, &testData);

    // Tape recording.
    cl::Independent(test.spread_);
    test.calculatePrices();
    cl::tape_function<double> f(test.spread_, test.totalPrice_);

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
    bool ok = test.check();

    ok &= testData.makeOutput();

    return ok;
}


test_suite* AdjointBermudanSwaptionTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Bermudan swaption tests");
    suite->add(QUANTLIB_TEST_CASE(&AdjointBermudanSwaptionTest::testTreeSwaptionEngine));
    suite->add(QUANTLIB_TEST_CASE(&AdjointBermudanSwaptionTest::testFdHullWhiteSwaptionEngine));
    return suite;
}


#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_bermudan_swaption)

BOOST_AUTO_TEST_CASE(testBermudanSwaptionTreeSwaptionEngine)
{
    BOOST_CHECK(AdjointBermudanSwaptionTest::testTreeSwaptionEngine());
}
BOOST_AUTO_TEST_CASE(testBermudanSwaptionFdHullWhiteSwaptionEngine)
{
    BOOST_CHECK(AdjointBermudanSwaptionTest::testFdHullWhiteSwaptionEngine());
}

BOOST_AUTO_TEST_SUITE_END()

#endif