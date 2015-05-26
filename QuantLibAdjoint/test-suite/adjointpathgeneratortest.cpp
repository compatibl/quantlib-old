/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2005 StatPro Italia srl
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

//based on pathgenerator.cpp from test-suite

#include "adjointpathgeneratortest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include <ql/methods/montecarlo/mctraits.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/processes/ornsteinuhlenbeckprocess.hpp>
#include <ql/processes/squarerootprocess.hpp>
#include <ql/processes/stochasticprocessarray.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/utilities/dataformatters.hpp>

using namespace QuantLib;
using namespace std;
using namespace boost::unit_test_framework;

namespace
{
    std::vector<Real> testSingle(const boost::shared_ptr<StochasticProcess1D>& process,
                                 const std::string& tag, bool brownianBridge,
                                 Real expected, Real antithetic)
    {
        typedef PseudoRandom::rsg_type rsg_type;
        typedef PathGenerator<rsg_type>::sample_type sample_type;
        std::vector<Real> calculated(2);
        BigNatural seed = 42;
        Time length = 10;
        Size timeSteps = 12;
        rsg_type rsg = PseudoRandom::make_sequence_generator(timeSteps, seed);
        PathGenerator<rsg_type> generator(process, length, timeSteps,
                                          rsg, brownianBridge);
        Size i;
        for (i = 0; i<100; i++)
            generator.next();

        sample_type sample = generator.next();
        calculated[0] = sample.value.back();
        Real error = std::fabs(calculated[0] - expected);
        Real tolerance = 2.0e-6;
        if (error > tolerance)
        {
            BOOST_ERROR("using " << tag << " process "
                        << (brownianBridge ? "with " : "without ")
                        << "brownian bridge:\n"
                        << std::setprecision(13)
                        << "    calculated: " << calculated[0] << "\n"
                        << "    expected:   " << expected << "\n"
                        << "    error:      " << error << "\n"
                        << "    tolerance:  " << tolerance);
        }

        sample = generator.antithetic();
        calculated[1] = sample.value.back();
        error = std::fabs(calculated[1] - antithetic);
        tolerance = 2.0e-6;
        if (error > tolerance)
        {
            BOOST_ERROR("using " << tag << " process "
                        << (brownianBridge ? "with " : "without ")
                        << "brownian bridge:\n"
                        << "antithetic sample:\n"
                        << std::setprecision(13)
                        << "    calculated: " << calculated[1] << "\n"
                        << "    expected:   " << antithetic << "\n"
                        << "    error:      " << error << "\n"
                        << "    tolerance:  " << tolerance);
        }
        return calculated;

    }

    std::vector<Real> testMultiple(const boost::shared_ptr<StochasticProcess>& process,
                                   const std::string& tag,
                                   Real expected[], Real antithetic[])
    {
        typedef PseudoRandom::rsg_type rsg_type;
        typedef MultiPathGenerator<rsg_type>::sample_type sample_type;

        BigNatural seed = 42;
        Time length = 10;
        Size timeSteps = 12;
        Size assets = process->size();
        std::vector<Real> calculated(2 * assets);
        rsg_type rsg = PseudoRandom::make_sequence_generator(timeSteps*assets,
                                                             seed);
        MultiPathGenerator<rsg_type> generator(process,
                                               TimeGrid(length, timeSteps),
                                               rsg, false);
        Size i, j;
        for (i = 0; i < 100; i++)
            generator.next();

        sample_type sample = generator.next();

        Real error, tolerance = 2.0e-7;
        for (j = 0; j < assets; j++)
            calculated[j] = sample.value[j].back();
        for (j = 0; j<assets; j++)
        {
            error = std::fabs(calculated[j] - expected[j]);
            if (error > tolerance)
            {
                BOOST_ERROR("using " << tag << " process "
                            << "(" << io::ordinal(j + 1) << " asset:)\n"
                            << std::setprecision(13)
                            << "    calculated: " << calculated[j] << "\n"
                            << "    expected:   " << expected[j] << "\n"
                            << "    error:      " << error << "\n"
                            << "    tolerance:  " << tolerance);
            }
        }

        sample = generator.antithetic();
        for (j = 0; j < assets; j++)
            calculated[assets + j] = sample.value[j].back();
        for (j = 0; j<assets; j++)
        {
            error = std::fabs(calculated[assets + j] - antithetic[j]);
            if (error > tolerance)
            {
                BOOST_ERROR("using " << tag << " process "
                            << "(" << io::ordinal(j + 1) << " asset:)\n"
                            << "antithetic sample:\n"
                            << std::setprecision(13)
                            << "    calculated: " << calculated[assets + j] << "\n"
                            << "    expected:   " << antithetic[j] << "\n"
                            << "    error:      " << error << "\n"
                            << "    tolerance:  " << tolerance);
            }
        }
        return calculated;
    }

    vector<vector<Real>> findTestSingleResult(Real sigmaValue)
    {
        Handle<Quote> x0(boost::shared_ptr<Quote>(new SimpleQuote(100.0)));
        Handle<YieldTermStructure> r(flatRate(0.05, Actual360()));
        Handle<YieldTermStructure> q(flatRate(0.02, Actual360()));
        Handle<BlackVolTermStructure> sigma(flatVol(sigmaValue, Actual360()));

        vector<vector<Real>> calculated_result;

        calculated_result.push_back(testSingle(boost::shared_ptr<StochasticProcess1D>(
            new BlackScholesMertonProcess(x0, q, r, sigma)),
            "Black-Scholes", false, 26.13784357783, 467.2928561411));

        calculated_result.push_back(testSingle(boost::shared_ptr<StochasticProcess1D>(
            new BlackScholesMertonProcess(x0, q, r, sigma)),
            "Black-Scholes", true, 60.28215549393, 202.6143139999));

        calculated_result.push_back(testSingle(boost::shared_ptr<StochasticProcess1D>(
            new OrnsteinUhlenbeckProcess(0.1, sigmaValue)),
            "Ornstein-Uhlenbeck", false, -0.8372003433557, 0.8372003433557));

        calculated_result.push_back(testSingle(boost::shared_ptr<StochasticProcess1D>(
            new SquareRootProcess(0.1, 0.1, sigmaValue, 10.0)),
            "square-root", false, 1.70608664108, 6.024200546031));
        return calculated_result;
    }
}




bool AdjointPathGeneratorTest::testPathGenerator()
{
    BOOST_TEST_MESSAGE("Testing 1-D path generation against cached values...");
    bool result = false;

#ifdef CL_TAPE_CPPAD

    SavedSettings backup;
    boost::timer timer;
    Settings::instance().evaluationDate() = Date(26, April, 2005);

    //Running Adjoint Differentiation for Sigma
    std::vector<cl::TapeDouble> sigma_Vector(1);
    sigma_Vector[0] = 0.20;

    //Beginning of tape recording
    Independent(sigma_Vector);

    vector<vector<Real>> calculated_result = findTestSingleResult(sigma_Vector[0]);
    Size n = calculated_result.size() * 2;
    std::vector<cl::TapeDouble> calculated(n);

    for (Size i = 0; i < calculated_result.size(); i++)
    {
        calculated[i * 2] = calculated_result[i][0];
        calculated[i * 2 + 1] = calculated_result[i][1];
    }

    cl::TapeFunction<double> f(sigma_Vector, calculated);
    double timeTapeRecording = timer.elapsed();
    //End of tape recording. Differentiaion wiil be held with respect to  the independent variables vector
    vector<double> sf_Forward, sf_Reverse;

    //Start differentiation in Forward mode
    double timeForward = gradForward(f, sf_Forward, false, false);

    //Start differentiation in Reverse mode
    double timeReverse = gradReverse(f, sf_Reverse, false, false);

    //Finite differences
    double h = 1.0e-10;
    timer.restart();
    std::vector<Real> sf_Finite(n);
    sigma_Vector[0] += h;

    vector<vector<Real>> calculated_result_h = findTestSingleResult(sigma_Vector[0]);
    for (Size i = 0; i < calculated_result.size(); i++)
    {
        sf_Finite[i * 2] = (calculated_result_h[i][0] - calculated_result[i][0]) / h;
        sf_Finite[i * 2 + 1] = (calculated_result_h[i][1] - calculated_result[i][1]) / h;
    }
    double timeAnalytical = timer.elapsed();

    //Check results
    double tol = 1e-4;
    result = checkWithFiniteDiff(sf_Forward, sf_Reverse, sf_Finite, tol, tol);
#endif
    return result;
}


vector<vector<Real>> findTestMultipleResult(Real sigmaValue)
{
    Handle<Quote> x0(boost::shared_ptr<Quote>(new SimpleQuote(100.0)));
    Handle<YieldTermStructure> r(flatRate(0.05, Actual360()));
    Handle<YieldTermStructure> q(flatRate(0.02, Actual360()));
    Handle<BlackVolTermStructure> sigma(flatVol(sigmaValue, Actual360()));

    Matrix correlation(3, 3);
    correlation[0][0] = 1.0; correlation[0][1] = 0.9; correlation[0][2] = 0.7;
    correlation[1][0] = 0.9; correlation[1][1] = 1.0; correlation[1][2] = 0.4;
    correlation[2][0] = 0.7; correlation[2][1] = 0.4; correlation[2][2] = 1.0;
    std::vector<boost::shared_ptr<StochasticProcess1D> > processes(3);
    boost::shared_ptr<StochasticProcess> process;
    vector<vector<Real>> calculated;
    processes[0] = boost::shared_ptr<StochasticProcess1D>(
        new BlackScholesMertonProcess(x0, q, r, sigma));
    processes[1] = boost::shared_ptr<StochasticProcess1D>(
        new BlackScholesMertonProcess(x0, q, r, sigma));
    processes[2] = boost::shared_ptr<StochasticProcess1D>(
        new BlackScholesMertonProcess(x0, q, r, sigma));
    process = boost::shared_ptr<StochasticProcess>(
        new StochasticProcessArray(processes, correlation));
    // commented values must be used when Halley's correction is enabled

    Real result1[] = {
        188.2235868185,
        270.6713069569,
        113.0431145652 };

    Real result1a[] = {
        64.89105742957,
        45.12494404804,
        108.0475146914 };

    calculated.push_back(testMultiple(process, "Black-Scholes", result1, result1a));

    processes[0] = boost::shared_ptr<StochasticProcess1D>(
        new OrnsteinUhlenbeckProcess(0.1, sigmaValue));
    processes[1] = boost::shared_ptr<StochasticProcess1D>(
        new OrnsteinUhlenbeckProcess(0.1, sigmaValue));
    processes[2] = boost::shared_ptr<StochasticProcess1D>(
        new OrnsteinUhlenbeckProcess(0.1, sigmaValue));
    process = boost::shared_ptr<StochasticProcess>(
        new StochasticProcessArray(processes, correlation));
    Real result2[] = {
        0.2942058437284,
        0.5525006418386,
        0.02650931054575 };
    Real result2a[] = {
        -0.2942058437284,
        -0.5525006418386,
        -0.02650931054575 };
    calculated.push_back(testMultiple(process, "Ornstein-Uhlenbeck", result2, result2a));

    processes[0] = boost::shared_ptr<StochasticProcess1D>(
        new SquareRootProcess(0.1, 0.1, sigmaValue, 10.0));
    processes[1] = boost::shared_ptr<StochasticProcess1D>(
        new SquareRootProcess(0.1, 0.1, sigmaValue, 10.0));
    processes[2] = boost::shared_ptr<StochasticProcess1D>(
        new SquareRootProcess(0.1, 0.1, sigmaValue, 10.0));
    process = boost::shared_ptr<StochasticProcess>(
        new StochasticProcessArray(processes, correlation));
    Real result3[] = {
        4.279510844897,
        4.943783503533,
        3.590930385958 };
    Real result3a[] = {
        2.763967737724,
        2.226487196647,
        3.503859264341 };
    calculated.push_back(testMultiple(process, "square-root", result3, result3a));
    return calculated;
}

bool AdjointPathGeneratorTest::testMultiPathGenerator()
{

    BOOST_TEST_MESSAGE("Testing n-D path generation against cached values...");
    bool result = false;

#ifdef CL_TAPE_CPPAD

    SavedSettings backup;
    boost::timer timer;
    Settings::instance().evaluationDate() = Date(26, April, 2005);
    //Running Adjoint Differentiation for Sigma
    std::vector<cl::TapeDouble> sigma_Vector(1);
    sigma_Vector[0] = 0.20;

    //Beginning of tape recording
    Independent(sigma_Vector);

    vector<vector<Real>> calculated_result = findTestMultipleResult(sigma_Vector[0]);
    Size n = 6 * calculated_result.size();
    std::vector<cl::TapeDouble> calculated(n);
    for (Size i = 0; i < calculated_result.size(); i++)
    for (Size j = 0; j < 6; j++)
        calculated[i * 6 + j] = calculated_result[i][j];


    cl::TapeFunction<double> f(sigma_Vector, calculated);
    double timeTapeRecording = timer.elapsed();
    //End of tape recording. Differentiaion wiil be held with respect to  the independent variables vector
    vector<double> sf_Forward, sf_Reverse;

    //Start differentiation in Forward mode
    double timeForward = gradForward(f, sf_Forward, false, false);

    //Start differentiation in Reverse mode
    double timeReverse = gradReverse(f, sf_Reverse, false, false);

    //Finite differences
    double h = 1.0e-10;
    timer.restart();
    std::vector<Real> sf_Finite(n);
    sigma_Vector[0] += h;
    vector<vector<Real>> calculated_result_h = findTestMultipleResult(sigma_Vector[0]);
    for (Size i = 0; i < calculated_result_h.size(); i++)
    for (Size j = 0; j < 6; j++)
        sf_Finite[i * 6 + j] = (calculated_result_h[i][j] - calculated_result[i][j]) / h;
    double timeAnalytical = timer.elapsed();

    //Check results
    double tol = 1e-4;
    result = checkWithFiniteDiff(sf_Forward, sf_Reverse, sf_Finite, tol, tol);
#endif
    return result;
}


test_suite* AdjointPathGeneratorTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Path generation tests");
    suite->add(QUANTLIB_TEST_CASE(&AdjointPathGeneratorTest::testPathGenerator));
    suite->add(QUANTLIB_TEST_CASE(&AdjointPathGeneratorTest::testMultiPathGenerator));
    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_path)

BOOST_AUTO_TEST_CASE(testPathGenerator)
{
    BOOST_CHECK(AdjointPathGeneratorTest::testPathGenerator());
}
BOOST_AUTO_TEST_CASE(testMultiPathGenerator)
{
    BOOST_CHECK(AdjointPathGeneratorTest::testMultiPathGenerator());
}

BOOST_AUTO_TEST_SUITE_END()

#endif