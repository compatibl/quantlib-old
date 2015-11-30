/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2013 Gary Kennedy
Copyright (C) 2015 Peter Caspers
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

// Based on blackformula.cpp file from Quantlib//test-suite.

#include "adjointblackformulaimpl.hpp"

using namespace QuantLib;
using namespace boost::unit_test_framework;

bool AdjointBlackFormulaTest::testBachelierImpliedVolMoneyness()
{
    BOOST_TEST_MESSAGE("Testing Bachelier implied vol dependence on moneyness...");

    Size n = 20;
    Test test(n, DataType::moneyness);

    cl::Independent(test.data_[DataType::moneyness]);
    test.calculateImplBpVol();
    cl::tape_function<double> f(test.data_[DataType::moneyness], test.implBpVol_);

    // Calculating derivatives using adjoint.
    test.adjointResults_ = f.Jacobian(test.indepDouble_);

    // Calculate derivatives using analytical formula.
    test.calcAnalytical();

    return test.checkAdjoint();
}


bool AdjointBlackFormulaTest::testBachelierImpliedVolMaturityTime()
{
    BOOST_TEST_MESSAGE("Testing Bachelier implied vol dependence on maturity time...");

    Size n = 20;
    Test test(n, DataType::tte);

    cl::Independent(test.data_[DataType::tte]);
    test.calculateImplBpVol();
    cl::tape_function<double> f(test.data_[DataType::tte], test.implBpVol_);

    // Calculating derivatives using adjoint.
    test.adjointResults_ = f.Jacobian(test.indepDouble_);

    // Calculate derivatives using analytical formula.
    test.calcAnalytical();

    return test.checkAdjoint();
}

bool AdjointBlackFormulaTest::testBachelierImpliedVolStdDev()
{
    BOOST_TEST_MESSAGE("Testing Bachelier implied vol dependence on stdDev...");

    Size n = 20;
    Test test(n, DataType::stdDev);

    cl::Independent(test.data_[DataType::stdDev]);
    test.calculateImplBpVol();
    cl::tape_function<double> f(test.data_[DataType::stdDev], test.implBpVol_);

    // Calculating derivatives using adjoint.
    test.adjointResults_ = f.Jacobian(test.indepDouble_);

    // Calculate derivatives using analytical formula.
    test.calcAnalytical();

    return test.checkAdjoint();
}

test_suite* AdjointBlackFormulaTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Adjoint Bachelier Black formula tests");

    suite->add(QUANTLIB_TEST_CASE(&AdjointBlackFormulaTest::testBachelierImpliedVolMaturityTime));
    suite->add(QUANTLIB_TEST_CASE(&AdjointBlackFormulaTest::testBachelierImpliedVolStdDev));
    suite->add(QUANTLIB_TEST_CASE(&AdjointBlackFormulaTest::testBachelierImpliedVolMoneyness));

    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_blackformula)

BOOST_AUTO_TEST_CASE(testBachelierImpliedVolMaturityTime)
{
    BOOST_CHECK(AdjointBlackFormulaTest::testBachelierImpliedVolMaturityTime());
}

BOOST_AUTO_TEST_CASE(testBachelierImpliedVolStdDev)
{
    BOOST_CHECK(AdjointBlackFormulaTest::testBachelierImpliedVolStdDev());
}

BOOST_AUTO_TEST_CASE(testBachelierImpliedVolMoneyness)
{
    BOOST_CHECK(AdjointBlackFormulaTest::testBachelierImpliedVolMoneyness());
}

BOOST_AUTO_TEST_SUITE_END()

#endif