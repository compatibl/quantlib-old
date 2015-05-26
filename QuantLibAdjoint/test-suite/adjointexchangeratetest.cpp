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

//based on exchangerate.cpp from test-suite

#include "adjointexchangeratetest.hpp"
#include "utilities.hpp"
#include <ql/exchangerate.hpp>
#include <ql/currencies/europe.hpp>
#include <ql/currencies/america.hpp>
#include <ql/currencies/asia.hpp>
#include <ql/currencies/exchangeratemanager.hpp>

using namespace QuantLib;
using namespace std;
using namespace boost::unit_test_framework;

bool AdjointExchangeRateTest::testDerived() {

    BOOST_TEST_MESSAGE("Testing derived exchange rates...");

    Currency EUR = EURCurrency(), USD = USDCurrency(), GBP = GBPCurrency();

    ExchangeRate eur_usd = ExchangeRate(EUR, USD, 1.2042);
    ExchangeRate eur_gbp = ExchangeRate(EUR, GBP, 0.6612);

    Money m1 = 50000.0 * GBP;
    Money m2 = 100000.0 * USD;

    Money::conversionType = Money::NoConversion;
    std::vector<cl::TapeDouble> rates = { eur_usd.rate(), eur_gbp.rate()};

    //Beginning of tape recording

    Independent(rates);
    ExchangeRate derived = ExchangeRate::chain(ExchangeRate(EUR, USD, rates[0]), ExchangeRate(EUR, GBP, rates[1]));
    Real calculated = derived.exchange(m1).value();
    //Real expected = Money(m1.value()* eur_usd / ur_gbp, USD).value();
    std::vector<cl::TapeDouble> output(1);
    output[0] = calculated;
    cl::TapeFunction<double > f(rates, output);

    //Start differentiation in Reverse mode
    vector<double> sy(1,1);
    vector<double> sf;
    sf = f.Reverse(1, sy);

    //Check with explicit analytical formula

    Real tolerance = 1e-2;
    bool result = true;
    if (std::abs(m1.value() / rates[1] - sf[0]) > tolerance * (std::abs(sf[0]) + 0.01)){
        BOOST_FAIL("Wrong result: \n"
            << "    expected:   " << m1.value() / rates[1] << "\n"
            << "    calculated: " << sf[0]);
        return false;
    }
     if(std::abs(-m1.value()*rates[0] / std::pow(rates[1], 2) - sf[1]) > tolerance * (std::abs(sf[1]) + 0.01)){
        BOOST_FAIL("Wrong result: \n"
            << "    expected:   " << -m1.value()*rates[0] / std::pow(rates[1], 2) << "\n"
            << "    calculated: " << sf[1]);
        result = false;
    }
    return result;
}

test_suite*  AdjointExchangeRateTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("CppAD Hull-White model calibration  tests");
    suite->add(QUANTLIB_TEST_CASE(&AdjointExchangeRateTest::testDerived));
    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_ExchaneRateTest)

BOOST_AUTO_TEST_CASE(testDerivedExchangeRate)
{
    BOOST_CHECK(AdjointExchangeRateTest::testDerived());
}

BOOST_AUTO_TEST_SUITE_END()

#endif