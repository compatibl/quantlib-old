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

#include "adjointgreeksimpl.hpp"

bool AdjointGreeksTest::test()
{
    BOOST_TEST_MESSAGE("Testing Greeks calculation with tape compression...");

    // Initializations...
    size_t n = 3;
    GreeksTestData td;
    GreeksTest test(n, &td.outPerform_);

    // Use only one input variable for tape recording.
    std::vector<cl::tobject> X = test.X_;

    // Declare an independent variable.
    cl::Independent(X);

    // Calculate the target function for only one input value (during tape recording).
    std::vector<cl::tobject> Y = { black_scholes_price(X) };

    // Declare tape function and stop the tape recording.
    cl::tape_function<cl::tvalue> f(X, Y);

    f.forward(0, test.x_val_);

    // Reverse mode calculation.
    std::vector<cl::tvalue> rev = f.reverse(1, std::vector<cl::tvalue>{ 1.0 });

    // Result checking and output generation...
    test.setReverseResults(rev);
    test.calcAnalytical();
    bool ok = test.checkAdjoint();
    ok &= td.makeOutput();
    return ok;
}


test_suite* AdjointGreeksTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Adjoint greeks test (with tape compression).");

    suite->add(QUANTLIB_TEST_CASE(&AdjointGreeksTest::test));

    return suite;
}


#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_greeks)

BOOST_AUTO_TEST_CASE(testGreeks)
{
    BOOST_CHECK(AdjointGreeksTest::test());
}

BOOST_AUTO_TEST_SUITE_END()

#endif
