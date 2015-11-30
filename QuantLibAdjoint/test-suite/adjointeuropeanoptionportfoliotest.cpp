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


#include "adjointeuropeanoptionportfolioimpl.hpp"


bool AdjointEuropeanOptionPortfolioTest::testDeltaCallPortfolio()
{
    BOOST_TEST_MESSAGE("Testing Delta calculations for Europeam Option Portfolio...");

    size_t n = testPortfolioSize;
    Test<Delta> test(n);

    // Tape recording.
    cl::Independent(test.data_[stock]);
    test.calculatePrices();
    cl::tape_function<double> f(test.data_[stock], test.totalPrice_);

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

    TestData<Delta> testData;

    return ok && testData.makeOutput();
}


bool AdjointEuropeanOptionPortfolioTest::testVegaCallPortfolio()
{
    BOOST_TEST_MESSAGE("Testing Vega calculations for Europeam Option Portfolio...");

    Test<Vega> test(testPortfolioSize);

    bool ok = test.test();

    TestData<Vega> testData;

    return ok && testData.makeOutput();
}


bool AdjointEuropeanOptionPortfolioTest::testThetaCallPortfolio()
{
    BOOST_TEST_MESSAGE("Testing Theta calculations for Europeam Option Portfolio...");

    Test<Theta> test(testPortfolioSize);

    bool ok = test.test();

    TestData<Theta> testData;

    return ok && testData.makeOutput();
}


bool AdjointEuropeanOptionPortfolioTest::testRhoCallPortfolio()
{
    BOOST_TEST_MESSAGE("Testing Rho calculations for Europeam Option Portfolio...");

    Test<Rho> test(testPortfolioSize);

    bool ok = test.test();

    TestData<Rho> testData;

    return ok && testData.makeOutput();

}


bool AdjointEuropeanOptionPortfolioTest::testGammaCallPortfolio()
{
    BOOST_TEST_MESSAGE("Testing Gamma calculations for Europeam Option Portfolio...");

    Test<Gamma> test(testPortfolioSize);

    bool ok = test.test();

    TestData<Gamma> testData;

    return ok && testData.makeOutput();

}


bool AdjointEuropeanOptionPortfolioTest::testVommaCallPortfolio()
{
    BOOST_TEST_MESSAGE("Testing Vomma calculations for Europeam Option Portfolio...");

    Test<Vomma> test(testPortfolioSize);

    bool ok = test.test();

    TestData<Vomma> testData;

    return ok && testData.makeOutput();

}


bool AdjointEuropeanOptionPortfolioTest::testVannaCallPortfolio()
{
    BOOST_TEST_MESSAGE("Testing Vanna calculations for Europeam Option Portfolio...");

    Test<Vanna> test(testPortfolioSize);

    bool ok = test.test();

    TestData<Vanna> testData;

    return ok && testData.makeOutput();

}


bool AdjointEuropeanOptionPortfolioTest::testCharmCallPortfolio()
{
    BOOST_TEST_MESSAGE("Testing Charm calculations for Europeam Option Portfolio...");

    Test<Charm> test(testPortfolioSize);

    bool ok = test.test();

    TestData<Charm> testData;

    return ok && testData.makeOutput();
}


bool AdjointEuropeanOptionPortfolioTest::testVetaCallPortfolio()
{
    BOOST_TEST_MESSAGE("Testing Veta calculations for Europeam Option Portfolio...");

    Test<Veta> test(testPortfolioSize);

    bool ok = test.test();

    TestData<Veta> testData;

    return ok && testData.makeOutput();
}


bool AdjointEuropeanOptionPortfolioTest::testVeraCallPortfolio()
{
    BOOST_TEST_MESSAGE("Testing Vera calculations for Europeam Option Portfolio...");

    Test<Vera> test(testPortfolioSize);

    bool ok = test.test();

    TestData<Vera> testData;

    return ok && testData.makeOutput();
}


bool AdjointEuropeanOptionPortfolioTest::testSpeedCallPortfolio()
{
    BOOST_TEST_MESSAGE("Testing Speed calculations for Europeam Option Portfolio...");

    Test<Speed> test(testPortfolioSize);

    bool ok = test.test();

    TestData<Speed> testData;

    return ok && testData.makeOutput();
}


bool AdjointEuropeanOptionPortfolioTest::testUltimaCallPortfolio()
{
    BOOST_TEST_MESSAGE("Testing Ultima calculations for Europeam Option Portfolio...");

    Test<Ultima> test(testPortfolioSize);

    bool ok = test.test();

    TestData<Ultima> testData;

    return ok && testData.makeOutput();
}


bool AdjointEuropeanOptionPortfolioTest::testZommaCallPortfolio()
{
    BOOST_TEST_MESSAGE("Testing Zomma calculations for Europeam Option Portfolio...");


    Test<Zomma> test(testPortfolioSize);

    bool ok = test.test();

    TestData<Zomma> testData;

    return ok && testData.makeOutput();
}


bool AdjointEuropeanOptionPortfolioTest::testColorCallPortfolio()
{
    BOOST_TEST_MESSAGE("Testing Color calculations for Europeam Option Portfolio...");

    Test<Color> test(testPortfolioSize);

    bool ok = test.test();

    TestData<Color> testData;

    return ok && testData.makeOutput();
}


test_suite* AdjointEuropeanOptionPortfolioTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Adjoint with European option portfolio tests");
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testDeltaCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testVegaCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testThetaCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testRhoCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testGammaCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testVommaCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testVannaCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testCharmCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testVetaCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testVeraCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testSpeedCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testUltimaCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testZommaCallPortfolio));
    suite->add(QUANTLIB_TEST_CASE(&AdjointEuropeanOptionPortfolioTest::testColorCallPortfolio));
    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_european_option_portfolio)

BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioDelta)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testDeltaCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioVega)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testVegaCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioTheta)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testThetaCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioRho)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testRhoCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioGamma)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testGammaCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioVomma)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testVommaCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioVanna)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testVannaCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioCharm)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testCharmCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioVeta)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testVetaCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioVera)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testVeraCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioSpeed)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testSpeedCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioUltima)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testUltimaCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioZomma)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testZommaCallPortfolio());
}
BOOST_AUTO_TEST_CASE(testEuropeanOptionPortfolioColor)
{
    BOOST_CHECK(AdjointEuropeanOptionPortfolioTest::testColorCallPortfolio());
}

BOOST_AUTO_TEST_SUITE_END()

#endif