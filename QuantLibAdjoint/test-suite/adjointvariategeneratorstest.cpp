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

#include "adjointvariategeneratorstest.hpp"

#include <boost/random.hpp>
#include "adjointtestutilities.hpp"

using namespace boost::unit_test_framework;
using namespace QuantLib;

namespace
{
    // Function getProperties(Distribution distribution)
    // Returns estimate of properties of variate generator with specefied Distribution.
    // Parameter: specified distribution.
    // Returns:
    //  properties[0] - Estimate of mean.
    //  properties[1] - Estimate of variance.
    template<typename Distribution>
    std::vector<Real> getProperties(Distribution distribution)
    {
        boost::mt19937 engine(42);
        boost::variate_generator<boost::mt19937&, Distribution> variateGenerator(engine, distribution);

        std::vector<Real> moments(2, 0.0);
        Size sampleSize = 100000;
        for (Size i = 0; i < sampleSize; i++)
        {
            // Variate.
            Real X = variateGenerator();
            // Estimate of first initial moment.
            moments[0] += X;
            // Estimate of second initial moment.
            moments[1] += std::pow(X, 2.0);
        }
        std::vector<Real> properties(2, 0.0);
        // Estimate of mean.
        properties[0] = moments[0] / sampleSize;
        // Estimate of variance.
        properties[1] = (moments[1] - std::pow(moments[0], 2.0) / sampleSize) / (sampleSize - 1);

        return properties;
    }

    // Checking consistency of vectors elements.
    bool checkWithAnalyticalDeriv(
        const std::vector<double>& adjointDeriv,
        const std::vector<Real>& analyticalDeriv,
        Real relativeTol,
        Real absoluteTol)
    {
        bool result = true;
        if (adjointDeriv.size() != analyticalDeriv.size())
        {
            result = false;
            BOOST_ERROR("\nAn adjoint derivatives vector and a analytical "
                        "derivatives vector have different sizes.");
        }
        for (size_t i = 0, n = adjointDeriv.size(); i < n; ++i)
        {
            Real err = std::abs(adjointDeriv[i] - analyticalDeriv[i]);
            Real toler = std::max(absoluteTol,
                                  relativeTol * std::max(std::abs(adjointDeriv[i]), std::abs(analyticalDeriv[i])));
            if (err > toler)
            {
                result = false;
                BOOST_ERROR("\nAdjoint derivative and analytical derivative at position " << i << " mismatch."
                            << "\n  adjoint: " << adjointDeriv[i]
                            << "\n  analytical: " << analyticalDeriv[i]
                            << "\n  tolerance: " << toler);
            }
        }
        return result;
    }
}

// Testing boost variate generator with different distributions with RealType = QuantLib::Real
// Properties of disributions (mean estimate and variance estimate) calculated by sample and then
// differentiated with respect to input value(s).
// Each test method returns true if adjoint derivatives are close enough to analytical derivatives,
// calculated using finite differences method; otherwise, it returns false.

// Method testExponentialVariate()
bool AdjointVariateGeneratorsTest::testExponentialVariate()
{
    BOOST_MESSAGE("Testing boost exponential distribution with Real...");

    std::vector<Real> lambda = { 4.0 };

    // Start of tape recording.
    // Mark \lambda as independent variable.
    cl::Independent(lambda);

    boost::random::exponential_distribution<Real> distribution(lambda[0]);

    std::vector<Real> properties = getProperties(distribution);

    // End of tape recording.
    // Differentiaion will be held with respect to the independent variables vector.
    cl::tape_function<double> f(lambda, properties);

    int propertiesNumber = properties.size();

    // Adjoint differentiation in Forward mode.
    std::vector<double> adjoint_derivatives(propertiesNumber);
    gradForward(f, adjoint_derivatives, false, false);

    // Calculation of derivatives using finite differences method.
    std::vector<Real> analytic_derivatives(properties.size());

    double h = 1e-8;

    lambda[0] += h;
    boost::random::exponential_distribution<Real> distributionfd(lambda[0]);
    std::vector<Real> propertiesfd = getProperties(distributionfd);
    lambda[0] -= h;

    for (int j = 0; j < properties.size(); j++)
    {
        analytic_derivatives[j] = (propertiesfd[j] - properties[j]) / h;
    }

    return checkWithAnalyticalDeriv(adjoint_derivatives, analytic_derivatives, 1e-4, 1e-10);
}

// Method testGammaVariate()
bool AdjointVariateGeneratorsTest::testGammaVariate()
{
    BOOST_MESSAGE("Testing boost gamma distribution with Real...");

    Real k = 10.0;
    Real theta = 3.0;
    std::vector<Real> parameters = { k, theta };

    // Start of tape recording.
    // Mark parameters as independent variables.
    cl::Independent(parameters);

    boost::random::gamma_distribution<Real> distribution(parameters[0], parameters[1]);

    std::vector<Real> properties = getProperties(distribution);

    // End of tape recording.
    // Differentiaion will be held with respect to the independent variables vector.
    cl::tape_function<double> f(parameters, properties);

    // Adjoint differentiation in Reverse mode.
    std::vector<double> adjoint_derivatives(properties.size() * parameters.size());
    gradReverse(f, adjoint_derivatives, false, false);

    // Calculation of derivatives using finite differences method.
    std::vector<Real> analytic_derivatives(properties.size() * parameters.size());

    double h = 1e-8;

    for (int i = 0; i < parameters.size(); i++)
    {
        parameters[i] += h;
        boost::random::gamma_distribution<Real> distributionfd(parameters[0], parameters[1]);
        std::vector<Real> propertiesfd = getProperties(distributionfd);
        parameters[i] -= h;

        for (int j = 0; j < properties.size(); j++)
        {
            analytic_derivatives[2 * j + i] = (propertiesfd[j] - properties[j]) / h;
        }
    }

    return checkWithAnalyticalDeriv(adjoint_derivatives, analytic_derivatives, 1e-3, 1e-10);
}

// Method testNormalVariate()
bool AdjointVariateGeneratorsTest::testNormalVariate()
{
    BOOST_MESSAGE("Testing boost normal distribution with Real...");

    Real mu = 5.0;
    Real sigma = 3.0;
    std::vector<Real> parameters = { mu, sigma };

    // Start of tape recording.
    // Mark parameters as independent variables.
    cl::Independent(parameters);

    boost::random::normal_distribution<Real> distribution(parameters[0], parameters[1]);

    std::vector<Real> properties = getProperties(distribution);

    // End of tape recording.
    // Differentiaion will be held with respect to the independent variables vector.
    cl::tape_function<double> f(parameters, properties);

    // Adjoint differentiation in Reverse mode.
    std::vector<double> adjoint_derivatives(properties.size() * parameters.size());
    gradReverse(f, adjoint_derivatives, false, false);

    // Calculation of derivatives using finite differences method.
    std::vector<Real> analytic_derivatives(properties.size() * parameters.size());

    double h = 1e-8;

    for (int i = 0; i < parameters.size(); i++)
    {
        parameters[i] += h;
        boost::random::normal_distribution<Real> distributionfd(parameters[0], parameters[1]);
        std::vector<Real> propertiesfd = getProperties(distributionfd);
        parameters[i] -= h;

        for (int j = 0; j < properties.size(); j++)
        {
            analytic_derivatives[2 * j + i] = (propertiesfd[j] - properties[j]) / h;
        }
    }

    return checkWithAnalyticalDeriv(adjoint_derivatives, analytic_derivatives, 1e-2, 1e-2);
}

// Method testLogNormalVariate()
bool AdjointVariateGeneratorsTest::testLogNormalVariate()
{
    BOOST_MESSAGE("Testing boost log-normal distribution with Real...");

    Real mu = 5.0;
    Real sigma = 3.0;
    std::vector<Real> parameters = { mu, sigma };

    // Start of tape recording.
    // Mark parameters as independent variables.
    cl::Independent(parameters);

    boost::random::lognormal_distribution<Real> distribution(parameters[0], parameters[1]);

    std::vector<Real> properties = getProperties(distribution);

    // End of tape recording.
    // Differentiaion will be held with respect to the independent variables vector.
    cl::tape_function<double> f(parameters, properties);

    // Adjoint differentiation in Reverse mode.
    std::vector<double> adjoint_derivatives(properties.size() * parameters.size());
    gradReverse(f, adjoint_derivatives, false, false);

    // Calculation of derivatives using finite differences method.
    std::vector<Real> analytic_derivatives(properties.size() * parameters.size());

    double h = 1e-8;

    for (int i = 0; i < parameters.size(); i++)
    {
        parameters[i] += h;
        boost::random::lognormal_distribution<Real> distributionfd(parameters[0], parameters[1]);
        std::vector<Real> propertiesfd = getProperties(distributionfd);
        parameters[i] -= h;

        for (int j = 0; j < properties.size(); j++)
        {
            analytic_derivatives[2 * j + i] = (propertiesfd[j] - properties[j]) / h;
        }
    }

    return checkWithAnalyticalDeriv(adjoint_derivatives, analytic_derivatives, 1e-4, 1e-10);
}

test_suite* AdjointVariateGeneratorsTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("AdjointVariateGenerators test");
    suite->add(QUANTLIB_TEST_CASE(&AdjointVariateGeneratorsTest::testExponentialVariate));
    suite->add(QUANTLIB_TEST_CASE(&AdjointVariateGeneratorsTest::testGammaVariate));
    suite->add(QUANTLIB_TEST_CASE(&AdjointVariateGeneratorsTest::testNormalVariate));
    suite->add(QUANTLIB_TEST_CASE(&AdjointVariateGeneratorsTest::testLogNormalVariate));
    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_variate_generator)

BOOST_AUTO_TEST_CASE(testExponentialVariate)
{
    BOOST_CHECK(AdjointVariateGeneratorsTest::testExponentialVariate());
}

BOOST_AUTO_TEST_CASE(testNormalVariate)
{
    BOOST_CHECK(AdjointVariateGeneratorsTest::testNormalVariate());
}

BOOST_AUTO_TEST_CASE(testGammaVariate)
{
    BOOST_CHECK(AdjointVariateGeneratorsTest::testGammaVariate());
}

BOOST_AUTO_TEST_CASE(testLogNormalVariate)
{
    BOOST_CHECK(AdjointVariateGeneratorsTest::testLogNormalVariate());
}

BOOST_AUTO_TEST_SUITE_END()

#endif
