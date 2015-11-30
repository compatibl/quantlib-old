/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2005, 2008 Klaus Spanderen
Copyright (C) 2007 StatPro Italia srl
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

#include "adjointbatesmodelimpl.hpp"

// These examples test adjoint differentiation of net present values (NPV) 
// on divident rates (testBatesModelOnDividentRate) 
// calculated with the Bates model (uses complex numbers).
//
// This particular set-up of parameters (see CommonVars in adjointbatesmodelimpl.hpp)
// is tuned to match the Black formula; the corresponding cross check is performed in the test.
// Returns true if
// 1) results are identical to the Black formula, and
// 2) both tangent and adjoint derivatives equals to finite difference derivative;
// otherwise returns false.
//
// Returns true if 
// 1) results are identical to the Black formula, and
// 2) both tangent and adjoint derivatives equals to finite difference derivative; 
// otherwise returns false.

bool AdjointBatesModelTest::testBatesModelOnDividentRate()
{
    BOOST_TEST_MESSAGE("Testing Bates model (like Black) on divident rates ...");

    bool ok = true;
    TestDataDividentRate testData;
    Size n = 3;
    BatesTestDividentRate test(n, &testData);

    // Check with BS formula
    ok &= test.CheckBatesVsBlack();

    // Tape recording
    cl::Independent(test.dividents_);
    test.calculateTotalNPV();
    cl::tape_function<double> f(test.dividents_, test.totalNPV_);

    // Forward mode derivatives calculation
    test.forwardResults_.resize(n);
    std::vector<double> dX(n, 0);
    for (size_t i = 0; i < n; i++)
    {
        dX[i] = 1;
        test.forwardResults_[i] = f.Forward(1, dX)[0];
        dX[i] = 0;
    }

    // Reverse mode derivatives calulation
    test.reverseResults_ = f.Reverse(1, std::vector<double>(1, 1));

    // Finite difference derivatives calulation
    test.calcAnalytical();

    // check all derivatives
    ok &= test.check();
    // ... and make output
    ok &= testData.makeOutput();

    Settings::instance().resetEvaluationDate();
    return ok;
}

test_suite* AdjointBatesModelTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Bates model tests");
    suite->add(QUANTLIB_TEST_CASE(&AdjointBatesModelTest::testBatesModelOnDividentRate));
    return suite;
}


#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_bates_model)

BOOST_AUTO_TEST_CASE(testBatesModelDividentRate)
{
    BOOST_CHECK(AdjointBatesModelTest::testBatesModelOnDividentRate());
}

BOOST_AUTO_TEST_SUITE_END()

#endif
