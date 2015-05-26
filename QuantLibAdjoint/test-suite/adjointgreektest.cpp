/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
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

/* Copyright (C) 2015 CompatibL */

#include "adjointgreektest.hpp"
#include "utilities.hpp"
#include <ql/quantlib.hpp>
#include <boost/timer.hpp>
#include <iostream>

#ifdef CL_TAPE_CPPAD
#include <cppad/local/ad_fun.hpp>
#endif

using namespace QuantLib;
using namespace boost::unit_test_framework;

bool AdjointGreekTest::testDelta()
{
    BOOST_TEST_MESSAGE("Testing the CppAD differentiation for Delta from the Black–Scholes formula...");
#ifdef CL_TAPE_CPPAD
       size_t n = 1;
    size_t m = 1;
    std::vector<Real::value_type> stock(n);
    stock[0] = 100.0;
    Real strike = 104.0;
    Real rate = 0.05;
    Real sigma = 0.20;
    Real time = 1;
    boost::timer timer;
    //Declare stock as an independent variable
    CppAD::Independent(stock);
    Real forward = stock[0] * std::exp(rate*time);
    Real stdDev = sigma*std::sqrt(time);
    BlackCalculator BC(Option::Call, strike, forward, stdDev);
    std::vector<Real::value_type> Price(1);
    //Dependent variables
    Price[0] = BC.value().value();
    // create f: x -> y and stop tape recording
    CppAD::ADFun<double> f(stock, Price);
    std::vector<double> sx(n);
    std::vector<double> sy(m);
    std::vector<double> sf;
    sx[0] = 1.;
    //Forward mode
    sf = f.Forward(1, sx);
    double seconds_f = timer.elapsed();
    std::cout << "\nTime in forward mode (s) " << seconds_f << std::endl;
    std::cout << "delta AD Forward  mode = " << sf[0] << std::endl;
    //Reverse mode
    sf.clear();
    sy[0] = 1.;
    sf = f.Reverse(1, sy);
    double seconds_r = timer.elapsed() - seconds_f;
    std::cout << "Time in reverse mode (s) " << seconds_r << std::endl;
    std::cout << "delta AD Reverse  mode = " << sf[0] << std::endl;
    Real deltaBS = BC.delta(stock[0]);
    double seconds_theor = timer.elapsed() - seconds_f - seconds_r;
    std::cout << "Time by theoretical formula (s) " << seconds_theor << std::endl;
    std::cout << "delta theoretic = " << deltaBS << std::endl;

    return std::fabs(sf[0] - deltaBS) < 0.01;
#endif
}

bool AdjointGreekTest::testGamma()
{
    BOOST_TEST_MESSAGE("Testing the CppAD differentiation for Gamma from the Black–Scholes formula...");
#ifdef CL_TAPE_CPPAD
    size_t n = 1;
    size_t m = 1;
    std::vector<Real::value_type> stock(n);
    stock[0] = 100.0;
    Real strike = 104.0;
    Real rate = 0.05;
    Real sigma = 0.20;
    Real time = 1;
    boost::timer timer;
    //Declare stock as an independent variable
    CppAD::Independent(stock);
    Real forward = stock[0] * std::exp(rate*time);
    Real stdDev = sigma*std::sqrt(time);
    BlackCalculator BC(Option::Call, strike, forward, stdDev);
    std::vector<Real::value_type> Price(1);
    //Dependent variables
    Price[0] = BC.value().value();
    // create f: x -> y and stop tape recording
    CppAD::ADFun<double> f(stock, Price);
    std::vector<double> sx(n);
    std::vector<double> sy(m);
    std::vector<double> sf;
    sx[0] = 1.;
    //Forward mode
    sf = f.Forward(1, sx);
    sx[0] = 0;
    sf = f.Forward(2, sx);
    double seconds_f = timer.elapsed();
    std::cout << "\nTime in forward mode (s) " << seconds_f << std::endl;
    std::cout << "gamma AD Forward  mode = " << 2*sf[0] << std::endl;
    //Reverse mode
    sf.clear();
    sy[0] = 1.;
    sf = f.Reverse(2, sy);
    double seconds_r = timer.elapsed() - seconds_f;
    std::cout << "Time in reverse mode (s) " << seconds_r << std::endl;
    std::cout << "gamma AD Reverse  mode = " << sf[1] << std::endl;
    Real deltaBS = BC.gamma(stock[0]);
    double seconds_theor = timer.elapsed() - seconds_f - seconds_r;
    std::cout << "Time by theoretical formula (s) " << seconds_theor << std::endl;
    std::cout << "gamma theoretic = " << deltaBS << std::endl;

    return std::fabs(sf[1] - deltaBS) < 0.01;
#endif
}

bool AdjointGreekTest::testVega()
{
    BOOST_TEST_MESSAGE("Testing the CppAD differentiation for Vega from the Black–Scholes formula...");
#ifdef CL_TAPE_CPPAD
    size_t n = 1;
    size_t m = 1;
    std::vector<Real::value_type> sigma(n);
    Real stock = 100.0;
    Real strike = 104.0;
    Real rate = 0.05;
    sigma[0] = 0.20;
    Real time = 1;
    boost::timer timer;
    //Declare stock as an independent variable
    CppAD::Independent(sigma);
    Real forward = stock * std::exp(rate*time);
    Real stdDev = sigma[0] *std::sqrt(time);
    BlackCalculator BC(Option::Call, strike, forward, stdDev);
    std::vector<Real::value_type> Price(1);
    //Dependent variables
    Price[0] = BC.value().value();
    // create f: x -> y and stop tape recording
    CppAD::ADFun<double> f(sigma, Price);
    std::vector<double> sx(n);
    std::vector<double> sy(m);
    std::vector<double> sf;
    sx[0] = 1.;
    //Forward mode
    sf = f.Forward(1, sx);
    double seconds_f = timer.elapsed();
    std::cout << "\nTime in forward mode (s) " << seconds_f << std::endl;
    std::cout << "vega AD Forward  mode = " << sf[0] << std::endl;
    //Reverse mode
    sf.clear();
    sy[0] = 1.;
    sf = f.Reverse(1, sy);
    double seconds_r = timer.elapsed() - seconds_f;
    std::cout << "Time in reverse mode (s) " << seconds_r << std::endl;
    std::cout << "vega AD Reverse  mode = " << sf[0] << std::endl;
    Real vegaBS = BC.vega(time);
    double seconds_theor = timer.elapsed() - seconds_f - seconds_r;
    std::cout << "Time by theoretical formula (s) " << seconds_theor << std::endl;
    std::cout << "vega theoretic = " << vegaBS << std::endl;

    return std::fabs(sf[0] - vegaBS) < 0.01;
#endif
}

bool AdjointGreekTest::testTheta()
{
    BOOST_TEST_MESSAGE("Testing the CppAD differentiation for Theta from the Black–Scholes formula...");
#ifdef CL_TAPE_CPPAD
    size_t n = 1;
    size_t m = 1;
    std::vector<Real::value_type> time(n);
    Real stock = 100.0;
    Real strike = 104.0;
    Real rate = 0.05;
    Real sigma = 0.20;
    time[0] = 1;
    boost::timer timer;
    //Declare stock as an independent variable
    CppAD::Independent(time);
    Real forward = stock * std::exp(rate*time[0]);
    Real stdDev = sigma * std::sqrt(time[0]);
    Real discount = std::exp(-rate * time[0]);
    BlackCalculator BC(Option::Call, strike, forward, stdDev, discount);
    std::vector<Real::value_type> Price(1);
    //Dependent variables
    Price[0] = BC.value().value();
    // create f: x -> y and stop tape recording
    CppAD::ADFun<double> f(time, Price);
    std::vector<double> sx(n);
    std::vector<double> sy(m);
    std::vector<double> sf;
    sx[0] = 1.;
    //Forward mode
    sf = f.Forward(1, sx);
    double seconds_f = timer.elapsed();
    std::cout << "\nTime in forward mode (s) " << seconds_f << std::endl;
    std::cout << "theta AD Forward  mode = " << sf[0] << std::endl;
    //Reverse mode
    sf.clear();
    sy[0] = 1.;
    sf = f.Reverse(1, sy);
    double seconds_r = timer.elapsed() - seconds_f;
    std::cout << "Time in reverse mode (s) " << seconds_r << std::endl;
    std::cout << "theta AD Reverse  mode = " << sf[0] << std::endl;
    Real thetaBS = BC.theta(stock, time[0]);
    double seconds_theor = timer.elapsed() - seconds_f - seconds_r;
    std::cout << "Time by theoretical formula (s) " << seconds_theor << std::endl;
    std::cout << "theta theoretic = " << thetaBS << std::endl;

    return std::fabs(-sf[0] - thetaBS) < 0.01;
#endif
}

bool AdjointGreekTest::testRho()
{
    BOOST_TEST_MESSAGE("Testing the CppAD differentiation for Rho from the Black–Scholes formula...");
#ifdef CL_TAPE_CPPAD
    size_t n = 1;
    size_t m = 1;
    std::vector<Real::value_type> rate(n);
    Real stock = 100.0;
    Real strike = 104.0;
    rate[0] = 0.05;
    Real sigma = 0.20;
    Real time = 1;
    boost::timer timer;
    //Declare stock as an independent variable
    CppAD::Independent(rate);
    Real forward = stock * std::exp(rate[0]*time);
    Real stdDev = sigma * std::sqrt(time);
    Real discount = std::exp(-rate[0] * time);
    BlackCalculator BC(Option::Call, strike, forward, stdDev, discount);
    std::vector<Real::value_type> Price(1);
    //Dependent variables
    Price[0] = BC.value().value();
    // create f: x -> y and stop tape recording
    CppAD::ADFun<double> f(rate, Price);
    std::vector<double> sx(n);
    std::vector<double> sy(m);
    std::vector<double> sf;
    sx[0] = 1.;
    //Forward mode
    sf = f.Forward(1, sx);
    double seconds_f = timer.elapsed();
    std::cout << "\nTime in forward mode (s) " << seconds_f << std::endl;
    std::cout << "rho AD Forward  mode = " << sf[0] << std::endl;
    //Reverse mode
    sf.clear();
    sy[0] = 1.;
    sf = f.Reverse(1, sy);
    double seconds_r = timer.elapsed() - seconds_f;
    std::cout << "Time in reverse mode (s) " << seconds_r << std::endl;
    std::cout << "rho AD Reverse  mode = " << sf[0] << std::endl;
    Real rhoBS = BC.rho(time);
    double seconds_theor = timer.elapsed() - seconds_f - seconds_r;
    std::cout << "Time by theoretical formula (s) " << seconds_theor << std::endl;
    std::cout << "rho theoretic = " << rhoBS << std::endl;

    return std::fabs(sf[0] - rhoBS) < 0.01;
#endif
}

test_suite* AdjointGreekTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("CppAD Greek tests");
    suite->add(QUANTLIB_TEST_CASE(&AdjointGreekTest::testDelta));
    suite->add(QUANTLIB_TEST_CASE(&AdjointGreekTest::testGamma));
    suite->add(QUANTLIB_TEST_CASE(&AdjointGreekTest::testRho));
    suite->add(QUANTLIB_TEST_CASE(&AdjointGreekTest::testTheta));
    suite->add(QUANTLIB_TEST_CASE(&AdjointGreekTest::testVega));
    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_greek)

    BOOST_AUTO_TEST_CASE(testDelta)
    {
        BOOST_CHECK(AdjointGreekTest::testDelta());
    }

    BOOST_AUTO_TEST_CASE(testGamma)
    {
        BOOST_CHECK(AdjointGreekTest::testGamma());
    }

    BOOST_AUTO_TEST_CASE(testRho)
    {
        BOOST_CHECK(AdjointGreekTest::testRho());
    }

    BOOST_AUTO_TEST_CASE(testTheta)
    {
        BOOST_CHECK(AdjointGreekTest::testTheta());
    }

    BOOST_AUTO_TEST_CASE(testVega)
    {
        BOOST_CHECK(AdjointGreekTest::testVega());
    }

BOOST_AUTO_TEST_SUITE_END()

#endif