/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2005, 2007, 2009, 2014 Klaus Spanderen
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

#include "adjointhestonprocessimpl.hpp"


bool AdjointHestonProcessTest::testHestonProcessPhi()
{
    BOOST_TEST_MESSAGE("Testing haracteristic functhion for the exact sampling of the Heston process ...");

    CommonVars data;
    std::vector<double> X0 = { 0.001, 0.02 };
    std::vector<Real> X = { X0[0], X0[1] };

    cl::Independent(X);

    std::complex<Real> var(X[0], X[1]);

    Real nu_0 = 0.5;
    Real nu_t = 0.75;
    Time dt = 0.3;

    std::complex<Real> haract = Phi1(data.process_, var, nu_0, nu_t, dt);

    std::vector<Real> Y(2);
    Y[0] = haract.real();
    Y[1] = haract.imag();

    cl::tape_function<double> f(X, Y);

    std::vector<double> adjDer = f.Jacobian(X0);

    // Check corectness of the differentiation without using the fact of existance of complex derivative.
    double shift = 0.0000001;
    double tol = 0.0001;
    std::vector<std::complex<Real>> comp;
    comp.push_back(std::complex<Real>(1, 0));
    comp.push_back(std::complex<Real>(0, 1));
    for (int i = 0; i < 2; i++)
    {
        std::vector<double> XshiftedUp = X0;
        std::vector<double> XshiftedDown = X0;
        XshiftedUp[i] += shift;
        XshiftedDown[i] -= shift;
        std::complex<Real> varUp(XshiftedUp[0], XshiftedUp[1]);
        std::complex<Real> varDown(XshiftedDown[0], XshiftedDown[1]);
        std::complex<Real> haractUp = Phi1(data.process_, varUp, nu_0, nu_t, dt);
        std::complex<Real> haractDown = Phi1(data.process_, varDown, nu_0, nu_t, dt);

        std::complex<Real> Dhar = (haractUp - haractDown) / 2 / shift / comp[i];

        if (std::abs(Dhar.real() - adjDer[3 * i]) > tol)
        {
            return false;
        }

        if (i)
        {
            if (std::abs(Dhar.imag() + adjDer[1]) > tol)
            {
                return false;
            }
        }
        else
        {
            if (std::abs(Dhar.imag() - adjDer[2]) > tol)
            {
                return false;
            }
        }
    }

    return true;
}

bool AdjointHestonProcessTest::testHestonProcessExpectation()
{
    BOOST_TEST_MESSAGE("Testing differentiation of expectation of the Heston process ...");

    TestDataRespectDeltaTime data(VariableName::expectation);
    int n = 5;
    HestonProcessRespectToDeltaTimeTest test(&data, n);

    cl::Independent(test.deltaTimes_);
    // Calculate volatilites.
    test.calculateVariables();
    cl::tape_function<double> f(test.deltaTimes_, test.SumVariables_);

    // Reverse mode derivatives calulation.
    test.reverseResults_ = f.Reverse(1, std::vector<double>(1, 1));

    // Forward mode derivatives calulation.
    test.forwardResults_.resize(n);
    std::vector<double> dX(n, 0);
    for (size_t i = 0; i < n; i++)
    {
        dX[i] = 1;
        test.forwardResults_[i] = f.Forward(1, dX)[0];
        dX[i] = 0;
    }

    // Calculate derivatives using analytical formula.
    test.calcAnalytical();

    return test.check() && data.makeOutput();
}

bool AdjointHestonProcessTest::testHestonProcessEvolve()
{
    BOOST_TEST_MESSAGE("Testing differentiation of evolve of the Heston process ...");

    TestDataRespectDeltaTime data(VariableName::evolve);
    int n = 5;
    HestonProcessRespectToDeltaTimeTest test(&data, n);

    cl::Independent(test.deltaTimes_);
    // Calculate volatilites.
    test.calculateVariables();
    cl::tape_function<double> f(test.deltaTimes_, test.SumVariables_);

    // Reverse mode derivatives calulation.
    test.reverseResults_ = f.Reverse(1, std::vector<double>(1, 1));

    // Forward mode derivatives calulation.
    test.forwardResults_.resize(n);
    std::vector<double> dX(n, 0);
    for (size_t i = 0; i < n; i++)
    {
        dX[i] = 1;
        test.forwardResults_[i] = f.Forward(1, dX)[0];
        dX[i] = 0;
    }

    // Calculate derivatives using analytical formula.
    test.calcAnalytical();

    return test.check() && data.makeOutput();
}

bool AdjointHestonProcessTest::testHestonProcessDrift()
{
    BOOST_TEST_MESSAGE("Testing differentiation of drift of the Heston process ...");

    CommonVars data;
    HestonProcessDriftTest test(&data);

    // Tape recording.
    cl::Independent(test.values_);
    test.calculateDrift();
    cl::tape_function<double> f(test.values_, test.drift_);

    // Calculates derivatives using adjoint.
    test.adjointResults_ = f.Jacobian(test.indepDouble_);

    test.calcAnalitical();
    return test.checkAdjoint();
}

bool AdjointHestonProcessTest::testHestonProcessStdDeviation()
{
    BOOST_TEST_MESSAGE("Testing differentiation of deviation of the Heston process ...");

    HestonProcessTestData data(VariableName::stdDeviation);
    HestonProcessMatrixFunctionsTest test(&data);

    // Tape recording.
    Independent(test.deltaTimes_);
    test.calculateMatrixRes();
    cl::tape_function<double> f(test.deltaTimes_, test.Results_);

    // Calculates derivatives using adjoint.
    int num = test.Results_.size();
    std::vector<double> w(1, 1);
    test.adjointResults_ = f.Forward(1, w);

    test.calcAnalitical();

    return test.checkAdjoint();
}

bool AdjointHestonProcessTest::testHestonProcessCovariance()
{
    BOOST_TEST_MESSAGE("Testing differentiation of covariance of the Heston process ...");

    HestonProcessTestData data(VariableName::covariance);
    HestonProcessMatrixFunctionsTest test(&data);

    // Tape recording.
    Independent(test.deltaTimes_);
    test.calculateMatrixRes();
    cl::tape_function<double> f(test.deltaTimes_, test.Results_);

    // Calculates derivatives using adjoint.
    int num = test.Results_.size();
    std::vector<double> w(1, 1);
    test.adjointResults_ = f.Forward(1, w);

    test.calcAnalitical();

    return test.checkAdjoint();
}



test_suite* AdjointHestonProcessTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Heston process tests");
    suite->add(QUANTLIB_TEST_CASE(&AdjointHestonProcessTest::testHestonProcessPhi));
    suite->add(QUANTLIB_TEST_CASE(&AdjointHestonProcessTest::testHestonProcessExpectation));
    suite->add(QUANTLIB_TEST_CASE(&AdjointHestonProcessTest::testHestonProcessDrift));
    suite->add(QUANTLIB_TEST_CASE(&AdjointHestonProcessTest::testHestonProcessEvolve));
    suite->add(QUANTLIB_TEST_CASE(&AdjointHestonProcessTest::testHestonProcessStdDeviation));
    suite->add(QUANTLIB_TEST_CASE(&AdjointHestonProcessTest::testHestonProcessCovariance));

    return suite;
}


#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_heston_process)

BOOST_AUTO_TEST_CASE(testHestonProcessPhi)
{
    BOOST_CHECK(AdjointHestonProcessTest::testHestonProcessPhi());
}

BOOST_AUTO_TEST_CASE(testHestonProcessExpectation)
{
    BOOST_CHECK(AdjointHestonProcessTest::testHestonProcessExpectation());
}

BOOST_AUTO_TEST_CASE(testHestonProcessDrift)
{
    BOOST_CHECK(AdjointHestonProcessTest::testHestonProcessDrift());
}

BOOST_AUTO_TEST_CASE(testHestonProcessEvolve)
{
    BOOST_CHECK(AdjointHestonProcessTest::testHestonProcessEvolve());
}

BOOST_AUTO_TEST_CASE(testHestonProcessStdDeviation)
{
    BOOST_CHECK(AdjointHestonProcessTest::testHestonProcessStdDeviation());
}

BOOST_AUTO_TEST_CASE(testHestonProcessCovariance)
{
    BOOST_CHECK(AdjointHestonProcessTest::testHestonProcessCovariance());
}


BOOST_AUTO_TEST_SUITE_END()

#endif