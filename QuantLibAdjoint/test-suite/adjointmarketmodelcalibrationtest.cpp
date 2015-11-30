/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2007 Ferdinando Ametrano
Copyright (C) 2007 Marco Bianchetti
Copyright (C) 2007 Cristina Duminuco
Copyright (C) 2007 Mark Joshi
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

// Based on marketmodel_smmcapletalphacalibration.cpp from test-suite.

// Algorithm of calibration can be found here:
// http://fbe.unimelb.edu.au/__data/assets/pdf_file/0019/806302/172.pdf

#include "adjointmarketmodelcalibrationimpl.hpp"

bool AdjointMarketModelCalibrationTest::testCapletVolatilities()
{
    BOOST_MESSAGE("Testing dependency of caplet volatilities from CTSMMCapletAlphaFormCalibration on input caplet volatilities...");
    Size n = 12;
    TestData testData;
    MarketModelTest test(n, &testData);

    // Start tape recording.
    cl::Independent(test.capletVols_);
    test.calibrate();

    // End tape recording.
    cl::tape_function<double> f(test.capletVols_, test.volsAndErrs_);

    // Calculating derivatives with adjoint.
    test.adjointResults_ = f.Jacobian(test.doubleCapletVols_);
    // Calculating derivatives using finite difference.
    test.calcAnalytical();

    return test.checkAdjoint() && testData.makeOutput();
}


test_suite* AdjointMarketModelCalibrationTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("AdjointMarketModelCalibration test");
    suite->add(QUANTLIB_TEST_CASE(&AdjointMarketModelCalibrationTest::testCapletVolatilities));
    return suite;
}


#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_market_model_calibration)

BOOST_AUTO_TEST_CASE(testMarketModelCalibration)
{
    BOOST_CHECK(AdjointMarketModelCalibrationTest::testCapletVolatilities());
}
BOOST_AUTO_TEST_SUITE_END()

#endif
