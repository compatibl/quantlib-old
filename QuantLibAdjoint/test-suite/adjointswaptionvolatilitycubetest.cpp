/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2006 Cristina Duminuco
Copyright (C) 2006, 2008 Ferdinando Ametrano
Copyright (C) 2006 Katiuscia Manzoni
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
#include "adjointswaptionvolatilitycubeimpl.hpp"


using namespace QuantLib;
using namespace boost::unit_test_framework;


bool AdjointSwaptionVolatilityCubeTest::testImpliedVolatilitiesDiffVols() {

    BOOST_TEST_MESSAGE("Testing swaption volatility cube Volatility (implied vols differentiation)...");

    SwaptionVolatilityCubeTestData testData(TestType::VolAtmVol);


    SwaptionVolatilityCubeTestDiffVols test(&testData);

    cl::Independent(test.vols_[test.optionTenorNumber_]);
    // Calculate volatilites.
    test.calculateBlackVar();
    cl::tape_function<double> f(test.vols_[test.optionTenorNumber_], test.BlackVariance_);

    // Calculating derivatives using adjoint.
    test.adjointResults_ = f.Jacobian(test.indepDouble_);

    // Calculate derivatives using the analytical formula.
    test.calcAnalytical();

    return test.checkAdjoint();
}

bool AdjointSwaptionVolatilityCubeTest::testBlackVariancesDiffVols() {

    BOOST_TEST_MESSAGE("Testing swaption volatility cube blackVariance (implied vols differentiation)...");

    SwaptionVolatilityCubeTestData testData(TestType::VolBlackVariance);

    SwaptionVolatilityCubeTestDiffVols test(&testData);

    cl::Independent(test.vols_[test.optionTenorNumber_]);
    // Calculate black variances.
    test.calculateBlackVar();
    cl::tape_function<double> f(test.vols_[test.optionTenorNumber_], test.BlackVariance_);

    // Calculating derivatives using adjoint.
    test.adjointResults_ = f.Jacobian(test.indepDouble_);

    // Calculate derivatives using analytical formula.
    test.calcAnalytical();

    return test.checkAdjoint();
}



test_suite* AdjointSwaptionVolatilityCubeTest::suite() {
    test_suite* suite = BOOST_TEST_SUITE("Swaption Volatility Cube tests");

    suite->add(QUANTLIB_TEST_CASE(&AdjointSwaptionVolatilityCubeTest::testImpliedVolatilitiesDiffVols));
    suite->add(QUANTLIB_TEST_CASE(&AdjointSwaptionVolatilityCubeTest::testBlackVariancesDiffVols));

    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_swaption_volatility_cube)

BOOST_AUTO_TEST_CASE(testSwaptionVolatilityCubeImpliedVolatilitiesDiffVols)
{
    BOOST_CHECK(AdjointSwaptionVolatilityCubeTest::testImpliedVolatilitiesDiffVols());
}

BOOST_AUTO_TEST_CASE(testSwaptionVolatilityCubeBlackVariancesDiffVols)
{
    BOOST_CHECK(AdjointSwaptionVolatilityCubeTest::testBlackVariancesDiffVols());
}

BOOST_AUTO_TEST_SUITE_END()

#endif