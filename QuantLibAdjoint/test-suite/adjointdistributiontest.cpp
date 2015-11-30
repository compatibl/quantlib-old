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

#include <ql/quantlib.hpp>
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include "adjointdistributiontest.hpp"
#include "boost\math\distributions.hpp"

using namespace QuantLib;
using namespace boost::unit_test_framework;

namespace
{
    bool checkWithAnalyticalFormula(std::vector<double> adjoint, std::vector<Real> analytical)
    {
        Size N;
        if (adjoint.size() != analytical.size())
            return false;
        else
            N = adjoint.size();
        bool result = true;
        for (Size i = 0; i < N; i++)
        {
            Real tollerance = 1e-4 * std::max(std::abs(adjoint[i]), std::abs(analytical[i]));
            if (std::abs(adjoint[i] - analytical[i]) > tollerance)
            {
                result = false;
                BOOST_FAIL("Wrong result in [" << i << "]:"
                           << "\n    expected(analytical):   " << analytical[i]
                           << "\n    calculated(adjoint): " << adjoint[i]);
            }
        }
        return result;
    }
}

// Beta distribution.
bool AdjointDistributionTest::testBetaDistribution()
{
    BOOST_TEST_MESSAGE("Testing adjoint differentiation for beta distribution...");

    bool result = false;

    // Set distribution parameters:
    // \alpha > 0
    // \beta > 0
    std::vector<cl::tape_double> alpha(1);
    alpha[0] = 3.0;

    Real beta = 2.5;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> meanvariance(2);

    // Start of tape recording.
    // Mark \alpha as independent variable.
    Independent(alpha);

    boost::math::beta_distribution<Real> betaDistribution(alpha[0], beta);

    // Analytical formula for Mean:
    // \frac{\alpha}{\alpha+\beta}
    meanvariance[0] = boost::math::mean(betaDistribution);

    // Analytical formula for Variance:
    // \frac{\alpha\beta}{(\alpha+\beta)^2(\alpha+\beta+1)}
    meanvariance[1] = boost::math::variance(betaDistribution);

    // End of tape recording.
    cl::tape_function<double> f(alpha, meanvariance);

    std::vector<double> sfForward;
    // Start of  differentiation in Forward mode.
    gradForward(f, sfForward, false, false);

    // Create vector of independent variables:
    // 0 < x < 1
    std::vector<cl::tape_double> X(1);
    X[0] = 0.7;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> cdf(1);

    // Start of tape recording.
    // Mark X as independent variable.
    Independent(X);

    // Analytical formula for CDF:
    // I_x(\alpha,\beta) is the regularized incomplete beta function.
    cdf[0] = boost::math::cdf(betaDistribution, X[0]);

    // End of tape recording.
    cl::tape_function<double> f2(X, cdf);

    std::vector<double> sfReverse;
    // Start of differentiation in Reverse mode.
    gradReverse(f2, sfReverse, false, false);

    // Create vector of adjoint derivatives.
    std::vector<double> sfAdjoint;
    sfAdjoint.push_back(sfForward[0]);
    sfAdjoint.push_back(sfForward[1]);
    sfAdjoint.push_back(sfReverse[0]);

    // Create vector of analytical derivatives.
    std::vector<Real> sfAnalytical;

    // Analytical formula for Mean derivative with respect to \alpha:
    // \frac{\beta}{(\alpha+\beta)^2}
    sfAnalytical.push_back(beta / ((alpha[0] + beta) * (alpha[0] + beta)));

    // Analytical formula for Variance derivative with respect to \alpha:
    // \frac{\beta(\beta^2-\beta\alpha+\beta-\alpha-2\alpha^2)}{(\alpha+\beta)^3(\alpha+\beta+1)^2}
    sfAnalytical.push_back((beta * (beta*beta - beta*alpha[0] + beta - alpha[0] - 2 * alpha[0] * alpha[0]))
                           / ((alpha[0] + beta) * (alpha[0] + beta) * (alpha[0] + beta) * (alpha[0] + beta + 1) * (alpha[0] + beta + 1)));

    // Analytical formula for PDF:
    // \frac{x^{\alpha-1}(1-x)^{\beta-1}}{\mathrm{B}(\alpha,\beta)}
    // \mathrm{B}(\alpha,\beta) is the beta function.
    sfAnalytical.push_back(boost::math::pdf(betaDistribution, X[0]));

    // Check with explicit analytical formulas.
    result = checkWithAnalyticalFormula(sfAdjoint, sfAnalytical);
    return result;
}

// Cauchy distribution.
bool AdjointDistributionTest::testCauchyDistribution()
{
    BOOST_TEST_MESSAGE("Testing adjoint differentiation for cauchy distribution...");

    bool result = false;

    // Set distribution parameters:
    // Location: x0
    // scale: \gamma > 0
    Real x0 = -2.0;
    Real gamma = 1.5;

    boost::math::cauchy_distribution<Real> cauchyDistribution(x0, gamma);

    // Create vector of independent variables:
    std::vector<cl::tape_double> X(1);
    X[0] = 12;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> cdf(1);

    // Start of tape recording.
    // Mark X as independent variable.
    Independent(X);

    // Analytical formula for CDF:
    // \frac{1}{\pi}\arctan\left(\frac{x - x_0}{\gamma}\right)+\frac{1}{2}
    cdf[0] = boost::math::cdf(cauchyDistribution, X[0]);

    // End of tape recording.
    cl::tape_function<double> f(X, cdf);

    std::vector<double> sfReverse;
    // Start of differentiation in Reverse mode.
    gradReverse(f, sfReverse, false, false);

    // Create vector of adjoint derivatives.
    std::vector<double> sfAdjoint;
    sfAdjoint.push_back(sfReverse[0]);


    // Create vector of analytical derivatives.
    std::vector<Real> sfAnalytical;

    // Analytical formula for PDF:
    // \frac{1}{\pi \gamma \left [1 + \left(\frac{x-x_0}{\gamma}\right)^2 \right ]}
    sfAnalytical.push_back(boost::math::pdf(cauchyDistribution, X[0]));

    // Check with explicit analytical formulas.
    result = checkWithAnalyticalFormula(sfAdjoint, sfAnalytical);
    return result;
}

// Chi-squared distribution.
bool AdjointDistributionTest::testChiSquaredDistribution()
{
    BOOST_TEST_MESSAGE("Testing adjoint differentiation for chi-squared distribution...");

    bool result = false;

    // Set distribution parameter:
    // Degrees of freedom: k \in N
    Real k = 2.0;

    boost::math::chi_squared_distribution<Real> chiSquaredDistribution(k);

    // Create vector of independent variables:
    // 0 \leq x < INF
    std::vector<cl::tape_double> X(1);
    X[0] = 8;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> cdf(1);

    // Start of tape recording.
    // Mark X as independent variable.
    Independent(X);

    // Analytical formula for CDF:
    // \frac{1}{\Gamma \left(\frac{k}{2} \right)}\gamma\left( \frac{k}{2},\frac{x}{2}\right)
    // \Gamma(k) is the gamma function.
    // \gamma(k, x) is the lower incomplete gamma function.
    cdf[0] = boost::math::cdf(chiSquaredDistribution, X[0]);

    // End of tape recording.
    cl::tape_function<double> f(X, cdf);

    std::vector<double> sfReverse;
    // Start of differentiation in Reverse mode.
    gradReverse(f, sfReverse, false, false);

    // Create vector of adjoint derivatives.
    std::vector<double> sfAdjoint;
    sfAdjoint.push_back(sfReverse[0]);


    // Create vector of analytical derivatives.
    std::vector<Real> sfAnalytical;

    // Analytical formula for PDF:
    // \frac{1}{2^\frac{k}{2}\Gamma \left(\frac{k}{2} \right)}x^{\frac{k}{2}-1}e^{-\frac{x}{2}}
    // \Gamma (k) is gamma function
    sfAnalytical.push_back(boost::math::pdf(chiSquaredDistribution, X[0]));

    // Check with explicit analytical formulas.
    result = checkWithAnalyticalFormula(sfAdjoint, sfAnalytical);
    return result;
}

// Exponential distribution.
bool AdjointDistributionTest::testExponentialDistribution()
{
    BOOST_TEST_MESSAGE("Testing adjoint differentiation for exponential distribution...");

    bool result = false;

    // Set distribution parameter:
    // \lambda > 0
    std::vector<cl::tape_double> lambda(1);
    lambda[0] = 2.0;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> meanvariance(2);

    // Start of tape recording.
    // Mark \lambda as independent variable.
    Independent(lambda);

    boost::math::exponential_distribution<Real> expDistribution(lambda[0]);

    // Analytical formula for Mean:
    // \lambda^{-1}
    meanvariance[0] = boost::math::mean(expDistribution);

    // Analytical formula for Variance:
    // \lambda^{-2}
    meanvariance[1] = boost::math::variance(expDistribution);

    // End of tape recording.
    cl::tape_function<double> f(lambda, meanvariance);

    std::vector<double> sfForward;
    // Start of  differentiation in Forward mode.
    gradForward(f, sfForward, false, false);

    // Create vector of independent variables:
    // x \geq 0
    std::vector<cl::tape_double> X(1);
    X[0] = 2.0;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> cdf(1);

    // Start of tape recording.
    // Mark X as independent variable.
    Independent(X);

    // Analytical formula for CDF:
    // 1 - e^{-\lambda x }
    cdf[0] = boost::math::cdf(expDistribution, X[0]);

    // End of tape recording.
    cl::tape_function<double> f2(X, cdf);

    std::vector<double> sfReverse;
    // Start of differentiation in Reverse mode.
    gradReverse(f2, sfReverse, false, false);

    // Create vector of adjoint derivatives.
    std::vector<double> sfAdjoint;
    sfAdjoint.push_back(sfForward[0]);
    sfAdjoint.push_back(sfForward[1]);
    sfAdjoint.push_back(sfReverse[0]);

    // Create vector of analytical derivatives.
    std::vector<Real> sfAnalytical;

    // Analytical formula for Mean derivative with respect to \lambda:
    // \frac{-1}{\lambda^2}
    sfAnalytical.push_back(-1 / (lambda[0] * lambda[0]));

    // Analytical formula for Variance derivative with respect to \lambda:
    // \frac{-2}{\lambda^3}
    sfAnalytical.push_back(-2 / (lambda[0] * lambda[0] * lambda[0]));

    // Analytical formula for PDF:
    // \lambda e^{-\lambda x}
    sfAnalytical.push_back(boost::math::pdf(expDistribution, X[0]));

    // Check with explicit analytical formulas.
    result = checkWithAnalyticalFormula(sfAdjoint, sfAnalytical);
    return result;
}

// Fisher-F distribution.
bool AdjointDistributionTest::testFisherFDistribution()
{
    BOOST_TEST_MESSAGE("Testing adjoint differentiation for fisher-f distribution...");

    bool result = false;

    // Set distribution parameters:
    // Degrees of freedom:
    // d_1 > 0
    // d_2 > 0 (For Mean: d_2 > 2, For Variance: d_2 > 4)
    Real d1 = 2.0;

    std::vector<cl::tape_double> d2(1);
    d2[0] = 5.0;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> meanvariance(2);

    // Start of tape recording.
    // Mark d2 as independent variable.
    Independent(d2);

    boost::math::fisher_f_distribution<Real> fisherfDistribution(d1, d2[0]);

    // Analytical formula for Mean:
    // \frac{d_2}{d_2-2}
    // d_2 > 2
    meanvariance[0] = boost::math::mean(fisherfDistribution);

    // Analytical formula for Variance:
    // \frac{2d_2^2(d_1 + d_2 - 2)}{d_1(d_2 - 2)^2(d_2 - 4)}
    // d_2 > 4
    meanvariance[1] = boost::math::variance(fisherfDistribution);

    // End of tape recording.
    cl::tape_function<double> f(d2, meanvariance);

    std::vector<double> sfForward;
    // Start of  differentiation in Forward mode.
    gradForward(f, sfForward, false, false);

    // Create vector of independent variables:
    // x \geq 0
    std::vector<cl::tape_double> X(1);
    X[0] = 3.0;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> cdf(1);

    // Start of tape recording.
    // Mark X as independent variable.
    Independent(X);


    // Analytical formula for CDF:
    // I_{\frac{d_1 x}{d_1 x + d_2}} \left(\tfrac{d_1}{2}, \tfrac{d_2}{2} \right)
    cdf[0] = boost::math::cdf(fisherfDistribution, X[0]);

    // End of tape recording.
    cl::tape_function<double> f2(X, cdf);

    std::vector<double> sfReverse;
    // Start of differentiation in Reverse mode.
    gradReverse(f2, sfReverse, false, false);

    // Create vector of adjoint derivatives.
    std::vector<double> sfAdjoint;
    sfAdjoint.push_back(sfForward[0]);
    sfAdjoint.push_back(sfForward[1]);
    sfAdjoint.push_back(sfReverse[0]);

    // Create vector of analytical derivatives.
    std::vector<Real> sfAnalytical;

    // Analytical formula for Mean derivative with respect to d_2:
    // \frac {-2}{(d_2-2)^2}
    sfAnalytical.push_back(-2 / ((d2[0] - 2)*(d2[0] - 2)));

    // Analytical formula for Variance derivative with respect to d_2:
    // -\frac {2 d_2(d_1(d_2^2+2d_2-16)+6d_2^2-28d_2+32)}{d_1(d_2-4)^2(d_2-2)^3}
    sfAnalytical.push_back(-(2 * d2[0] * (d1*(d2[0] * d2[0] + 2 * d2[0] - 16) + 6 * d2[0] * d2[0] - 28 * d2[0] + 32))
                           / (d1*(d2[0] - 4)*(d2[0] - 4)*(d2[0] - 2)*(d2[0] - 2)*(d2[0] - 2)));

    // Analytical formula for PDF:
    // \frac{\sqrt{\frac{(d_1\,x)^{d_1}\,\,d_2^{d_2}}{(d_1\,x+d_2)^{d_1+d_2}}}}
    //      {x\mathrm {B}\left(\frac { d_1 }{2}, \frac { d_2 }{2}\right)}
    // \mathrm{B}(\alpha,\beta) is the beta function.
    sfAnalytical.push_back(boost::math::pdf(fisherfDistribution, X[0]));

    // Check with explicit analytical formulas.
    result = checkWithAnalyticalFormula(sfAdjoint, sfAnalytical);
    return result;
}

// Gamma distribution.
bool AdjointDistributionTest::testGammaDistribution()
{
    BOOST_TEST_MESSAGE("Testing adjoint differentiation for gamma distribution...");

    bool result = false;

    // Set distribution parameters:
    // Shape : k > 0
    // Scale : \Theta > 0
    Real k = 2.0;

    std::vector<cl::tape_double> theta(1);
    theta[0] = 2.0;



    // Create vector of dependent variables.
    std::vector<cl::tape_double> meanvariance(2);

    // Start of tape recording.
    // Mark \Theta as independent variable.
    Independent(theta);

    boost::math::gamma_distribution<Real> gammaDistribution(k, theta[0]);

    // Analytical formula for Mean:
    // k \Theta
    meanvariance[0] = boost::math::mean(gammaDistribution);

    // Analytical formula for Variance:
    // k \Theta^{-2}
    meanvariance[1] = boost::math::variance(gammaDistribution);

    // End of tape recording.
    cl::tape_function<double> f(theta, meanvariance);

    std::vector<double> sfForward;
    // Start of  differentiation in Forward mode.
    gradForward(f, sfForward, false, false);

    // Create vector of independent variables:
    // x > 0
    std::vector<cl::tape_double> X(1);
    X[0] = 2.0;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> cdf(1);

    // Start of tape recording.
    // Mark X as independent variable.
    Independent(X);


    // Analytical formula for CDF:
    // \frac{1}{\Gamma(k)} \gamma\left( k, \frac{x}{\Theta}\right)
    // \Gamma(k) is the gamma function.
    // \gamma\left( k, \frac{x}{\Theta}\right) is the lower incomplete gamma function.
    cdf[0] = boost::math::cdf(gammaDistribution, X[0]);

    // End of tape recording.
    cl::tape_function<double> f2(X, cdf);

    std::vector<double> sfReverse;
    // Start of differentiation in Reverse mode.
    gradReverse(f2, sfReverse, false, false);

    // Create vector of adjoint derivatives.
    std::vector<double> sfAdjoint;
    sfAdjoint.push_back(sfForward[0]);
    sfAdjoint.push_back(sfForward[1]);
    sfAdjoint.push_back(sfReverse[0]);

    // Create vector of analytical derivatives.
    std::vector<Real> sfAnalytical;

    // Analytical formula for Mean derivative with respect to \Theta:
    // k
    sfAnalytical.push_back(k);

    // Analytical formula for Variance derivative with respect to \Theta:
    // 2k\Theta
    sfAnalytical.push_back(2 * k * theta[0]);

    // Analytical formula for PDF:
    // \frac{1}{\Gamma(k)\Theta^k}x^{k-1}e^{\frac{-x}{\Theta}}
    // \Gamma(k) is the gamma function.
    sfAnalytical.push_back(boost::math::pdf(gammaDistribution, X[0]));

    // Check with explicit analytical formulas.
    result = checkWithAnalyticalFormula(sfAdjoint, sfAnalytical);
    return result;
}

// Laplace distribution.
bool AdjointDistributionTest::testLaplaceDistribution()
{
    BOOST_TEST_MESSAGE("Testing adjoint differentiation for laplace distribution...");

    bool result = false;

    // Set distribution parameters:
    // Location : \mu
    // Scale : b > 0
    Real mu = 1.0;
    Real b = 2.0;

    boost::math::laplace_distribution<Real> laplaceDistribution(mu, b);

    // Create vector of independent variables:
    std::vector<cl::tape_double> X(1);
    X[0] = -2.0;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> cdf(1);

    // Start of tape recording.
    // Mark X as independent variable.
    Independent(X);

    // Analytical formula for CDF:
    //  \begin { cases }
    //  \frac12 \exp \left(\frac { x - \mu }{b} \right) & \mbox { if }x < \mu \\
        //  1 - \frac12 \exp \left(-\frac { x - \mu }{b} \right) & \mbox { if }x \geq \mu
    //  \end { cases }
    cdf[0] = boost::math::cdf(laplaceDistribution, X[0]);

    // End of tape recording.
    cl::tape_function<double> f2(X, cdf);

    std::vector<double> sfReverse;
    // Start of differentiation in Reverse mode.
    gradReverse(f2, sfReverse, false, false);

    // Create vector of adjoint derivatives.
    std::vector<double> sfAdjoint;
    sfAdjoint.push_back(sfReverse[0]);

    // Create vector of analytical derivatives.
    std::vector<Real> sfAnalytical;

    // Analytical formula for PDF:
    // \frac{1}{2b} \exp \left(-\frac{|x-\mu|}b \right)
    sfAnalytical.push_back(boost::math::pdf(laplaceDistribution, X[0]));

    // Check with explicit analytical formulas.
    result = checkWithAnalyticalFormula(sfAdjoint, sfAnalytical);
    return result;
}

// Log-normal distribution.
bool AdjointDistributionTest::testLognormalDistribution()
{
    BOOST_TEST_MESSAGE("Testing adjoint differentiation for log-normal distribution...");

    bool result = false;

    // Set distribution parameters:
    // Location: \mu
    // Scale: \sigma > 0
    Real mu = 2.0;

    std::vector<cl::tape_double> sigma(1);
    sigma[0] = 2.0;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> meanvariance(2);

    // Start of tape recording.
    // Mark \sigma as independent variable.
    Independent(sigma);

    boost::math::lognormal_distribution<Real> lognormalDistribution(mu, sigma[0]);

    // Analytical formula for Mean:
    // e^{\mu+\sigma^2/2}
    meanvariance[0] = boost::math::mean(lognormalDistribution);

    // Analytical formula for Variance:
    // (e^{\sigma^2}-1) e^{2\mu+\sigma^2}
    meanvariance[1] = boost::math::variance(lognormalDistribution);

    // End of tape recording.
    cl::tape_function<double> f(sigma, meanvariance);

    std::vector<double> sfForward;
    // Start of  differentiation in Forward mode.
    gradForward(f, sfForward, false, false);

    // Create vector of independent variables:
    // x > 0
    std::vector<cl::tape_double> X(1);
    X[0] = 2.0;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> cdf(1);

    // Start of tape recording.
    // Mark X as independent variable.
    Independent(X);

    // Analytical formula for CDF:
    //  \frac{1}{2} + \frac{1}{2}\mathrm{erf}\left[\frac{\ln x-\mu}{\sqrt{2}\sigma}\right]
    cdf[0] = boost::math::cdf(lognormalDistribution, X[0]);

    // End of tape recording.
    cl::tape_function<double> f2(X, cdf);

    std::vector<double> sfReverse;
    // Start of differentiation in Reverse mode.
    gradReverse(f2, sfReverse, false, false);

    // Create vector of adjoint derivatives.
    std::vector<double> sfAdjoint;
    sfAdjoint.push_back(sfForward[0]);
    sfAdjoint.push_back(sfForward[1]);
    sfAdjoint.push_back(sfReverse[0]);

    // Create vector of analytical derivatives.
    std::vector<Real> sfAnalytical;

    // Analytical formula for Mean derivative with respect to \sigma:
    // \sigma e^{\mu+\sigma^2/2}
    sfAnalytical.push_back(sigma[0] * std::exp(mu + sigma[0] * sigma[0] / 2));

    // Analytical formula for Variance derivative with respect to \sigma:
    // 2(2e^{\sigma^2}-1)\sigma e^{2\mu+\sigma^2}
    sfAnalytical.push_back(2 * (2 * std::exp(sigma[0] * sigma[0]) - 1)*sigma[0] * std::exp(2 * mu + sigma[0] * sigma[0]));

    // Analytical formula for PDF:
    // \frac{1}{x\sigma\sqrt{2\pi}}\ e^{-\frac{\left(\ln x-\mu\right)^2}{2\sigma^2}}
    sfAnalytical.push_back(boost::math::pdf(lognormalDistribution, X[0]));

    // Check with explicit analytical formulas.
    result = checkWithAnalyticalFormula(sfAdjoint, sfAnalytical);

    return result;
}

// Normal distribution.
bool AdjointDistributionTest::testNormalDistribution()
{
    BOOST_TEST_MESSAGE("Testing adjoint differentiation for normal distribution...");

    bool result = false;

    // Set distribution parameters:
    // Mean: \mu
    // Variance: \sigma^2 > 0
    Real mu = 2.0;
    Real sigma = 3.0;

    boost::math::normal_distribution<Real> normalDistribution(mu, sigma);

    // Create vector of independent variables:
    std::vector<cl::tape_double> X(1);
    X[0] = 3;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> cdf(1);

    // Start of tape recording.
    // Mark X as independent variable.
    Independent(X);

    // Analytical formula for CDF:
    // \frac{1}{2}\left[1 + \operatorname{erf}\left( \frac{x-\mu}{\sigma\sqrt{2}}\right)\right]
    cdf[0] = boost::math::cdf(normalDistribution, X[0]);

    // End of tape recording.
    cl::tape_function<double> f(X, cdf);

    std::vector<double> sfReverse;
    // Start of differentiation in Reverse mode.
    gradReverse(f, sfReverse, false, false);

    // Create vector of adjoint derivatives.
    std::vector<double> sfAdjoint;
    sfAdjoint.push_back(sfReverse[0]);


    // Create vector of analytical derivatives.
    std::vector<Real> sfAnalytical;

    // Analytical formula for PDF:
    // \frac{1}{\sigma\sqrt{2\pi}} e^{-\frac{(x - \mu)^2}{2 \sigma^2}}
    sfAnalytical.push_back(boost::math::pdf(normalDistribution, X[0]));

    // Check with explicit analytical formulas.
    result = checkWithAnalyticalFormula(sfAdjoint, sfAnalytical);
    return result;
}

// Pareto distribution.
bool AdjointDistributionTest::testParetoDistribution()
{
    BOOST_TEST_MESSAGE("Testing adjoint differentiation for pareto distribution...");

    bool result = false;

    // Set distribution parameters:
    // Scale: x_m > 0
    // Shape: \alpha > 0
    Real xm = 2.0;

    std::vector<cl::tape_double> alpha(1);
    alpha[0] = 3.0;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> meanvariance(2);

    // Start of tape recording.
    // Mark \alpha as independent variable.
    Independent(alpha);

    boost::math::pareto_distribution<Real> paretoDistribution(xm, alpha[0]);

    // Analytical formula for Mean:
    // \begin{cases}
    //     \infty & \text{for } \alpha \le 1 \\
        //     \frac { \alpha\, x_m }{\alpha - 1} & \text{for } \alpha > 1
    // \end{cases}
    meanvariance[0] = boost::math::mean(paretoDistribution);

    // Analytical formula for Variance:
    // \begin{cases}
    //     \infty & \text{for }\alpha\in(1,2] \\
        //     \frac{x_\mathrm{m}^2\alpha}{(\alpha-1)^2(\alpha-2)} & \text{for }\alpha>2
    // \end{cases}
    meanvariance[1] = boost::math::variance(paretoDistribution);

    // End of tape recording.
    cl::tape_function<double> f(alpha, meanvariance);

    std::vector<double> sfForward;
    // Start of  differentiation in Forward mode.
    gradForward(f, sfForward, false, false);

    // Create vector of independent variables:
    // x \geq x_m
    std::vector<cl::tape_double> X(1);
    X[0] = 4.0;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> cdf(1);

    // Start of tape recording.
    // Mark X as independent variable.
    Independent(X);

    // Analytical formula for CDF:
    // 1-\left(\frac{x_m}{x}\right)^\alpha
    cdf[0] = boost::math::cdf(paretoDistribution, X[0]);

    // End of tape recording.
    cl::tape_function<double> f2(X, cdf);

    std::vector<double> sfReverse;
    // Start of differentiation in Reverse mode.
    gradReverse(f2, sfReverse, false, false);

    // Create vector of adjoint derivatives.
    std::vector<double> sfAdjoint;
    sfAdjoint.push_back(sfForward[0]);
    sfAdjoint.push_back(sfForward[1]);
    sfAdjoint.push_back(sfReverse[0]);

    // Create vector of analytical derivatives.
    std::vector<Real> sfAnalytical;

    // Analytical formula for Mean derivative with respect to \alpha:
    // -\frac{x_m}{(\alpha-1)^2}
    sfAnalytical.push_back(-xm / ((alpha[0] - 1)*(alpha[0] - 1)));

    // Analytical formula for Variance derivative with respect to \alpha:
    // -\frac{2x_m^2(\alpha^2-\alpha-1)}{(\alpha-2)^2(\alpha-1)^3}
    sfAnalytical.push_back(-(2 * xm*xm*(alpha[0] * alpha[0] - alpha[0] - 1))
                           / ((alpha[0] - 2)*(alpha[0] - 2)*(alpha[0] - 1)*(alpha[0] - 1)*(alpha[0] - 1)));

    // Analytical formula for PDF:
    // \frac{\alpha\,x_m^\alpha}{x^{\alpha+1}}
    sfAnalytical.push_back(boost::math::pdf(paretoDistribution, X[0]));

    // Check with explicit analytical formulas.
    result = checkWithAnalyticalFormula(sfAdjoint, sfAnalytical);
    return result;
}

// Student`s-t distribution.
bool AdjointDistributionTest::testStudentsTDistribution()
{
    BOOST_TEST_MESSAGE("Testing adjoint differentiation for student`s-t  distribution...");

    bool result = false;

    // Set distribution parameter:
    // Degrees od freedom: \nu > 0
    Real nu = 7;

    boost::math::students_t_distribution<Real> studentstDistribution(nu);

    // Create vector of independent variables:
    std::vector<cl::tape_double> X(1);
    X[0] = 3;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> cdf(1);

    // Start of tape recording.
    // Mark X as independent variable.
    Independent(X);

    // Analytical formula for CDF:
    // \frac{1}{2} + x \Gamma \left( \frac{\nu+1}{2} \right)
    // \frac{_2F_1 \left( \frac{1}{2};\frac{\nu+1}{2};\frac{3}{2};-\frac{x^2}{\nu} \right)}
    //     {\sqrt{\pi\nu}\,\Gamma \left(\frac{\nu}{2}\right)}
    // _2F_1() is the hypergeometric function
    // \Gamma() is the gamma function
    cdf[0] = boost::math::cdf(studentstDistribution, X[0]);

    // End of tape recording.
    cl::tape_function<double> f(X, cdf);

    std::vector<double> sfReverse;
    // Start of differentiation in Reverse mode.
    gradReverse(f, sfReverse, false, false);

    // Create vector of adjoint derivatives.
    std::vector<double> sfAdjoint;
    sfAdjoint.push_back(sfReverse[0]);


    // Create vector of analytical derivatives.
    std::vector<Real> sfAnalytical;

    // Analytical formula for PDF:
    // \frac{\Gamma \left(\frac{\nu+1}{2} \right)} {\sqrt{\nu\pi}\,\Gamma \left(\frac{\nu}{2} \right)}/
    // left(1+\frac{x^2}{\nu} \right)^{-\frac{\nu+1}{2}}
    // \Gamma() is the gamma function
    sfAnalytical.push_back(boost::math::pdf(studentstDistribution, X[0]));

    // Check with explicit analytical formulas.
    result = checkWithAnalyticalFormula(sfAdjoint, sfAnalytical);
    return result;
}

// Inverse chi-square distribution.
bool AdjointDistributionTest::testInverseChiSquaredDistribution()
{
    BOOST_TEST_MESSAGE("Testing adjoint differentiation for inverse chi-square distribution...");

    bool result = false;

    // Set distribution parameters:
    // Degrees of freedom: \nu > 0

    std::vector<cl::tape_double> nu(1);
    nu[0] = 5.0;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> meanvariance(2);

    // Start of tape recording.
    // Mark \nu as independent variable.
    Independent(nu);

    boost::math::inverse_chi_squared_distribution<Real> invChiSquareDistribution(nu[0]);

    // Analytical formula for Mean:
    // \frac{1}{\nu-2}
    // \nu > 2
    meanvariance[0] = boost::math::mean(invChiSquareDistribution);

    // Analytical formula for Variance:
    // \frac{2}{(\nu-2)^2 (\nu-4)}
    meanvariance[1] = boost::math::variance(invChiSquareDistribution);

    // End of tape recording.
    cl::tape_function<double> f(nu, meanvariance);

    std::vector<double> sfForward;
    // Start of  differentiation in Forward mode.
    gradForward(f, sfForward, false, false);

    // Create vector of independent variables:
    // x > 0
    std::vector<cl::tape_double> X(1);
    X[0] = 4.0;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> cdf(1);

    // Start of tape recording.
    // Mark X as independent variable.
    Independent(X);

    // Analytical formula for CDF:
    // \Gamma\left(\frac{\nu}{2},\frac{1}{2x}\right)\bigg/ \Gamma\left(\frac{\nu}{2}\right)
    // \Gamma() is the gamma fnction
    cdf[0] = boost::math::cdf(invChiSquareDistribution, X[0]);

    // End of tape recording.
    cl::tape_function<double> f2(X, cdf);

    std::vector<double> sfReverse;
    // Start of differentiation in Reverse mode.
    gradReverse(f2, sfReverse, false, false);

    // Create vector of adjoint derivatives.
    std::vector<double> sfAdjoint;
    sfAdjoint.push_back(sfForward[0]);
    sfAdjoint.push_back(sfForward[1]);
    sfAdjoint.push_back(sfReverse[0]);

    // Create vector of analytical derivatives.
    std::vector<Real> sfAnalytical;

    // Analytical formula for Mean derivative with respect to \nu:
    // -\frac{1}{(\nu-2)^2}
    sfAnalytical.push_back(-1 / ((nu[0] - 2)*(nu[0] - 2)));

    // Analytical formula for Variance derivative with respect to \nu:
    // \frac{20 - 6\nu}{(\nu-4)^2 (\nu-2)^3}
    sfAnalytical.push_back((20 - 6 * nu[0]) / ((nu[0] - 4)*(nu[0] - 4)*(nu[0] - 2)*(nu[0] - 2)*(nu[0] - 2)));

    // Analytical formula for PDF:
    // \frac{2^{-\nu/2}}{\Gamma(\nu/2)}x^{-\nu/2-1}e^{-1/(2 x)}
    // \Gamma() is the gamma fnction
    sfAnalytical.push_back(boost::math::pdf(invChiSquareDistribution, X[0]));

    // Check with explicit analytical formulas.
    result = checkWithAnalyticalFormula(sfAdjoint, sfAnalytical);

    return result;
}

// Inverse gamma distribution.
bool AdjointDistributionTest::testInverseGammaDistribution()
{
    BOOST_TEST_MESSAGE("Testing adjoint differentiation for inverse gamma distribution...");

    bool result = false;

    // Set distribution parameters:
    // Shape: \alpha > 0
    // Scale: \beta > 0

    Real alpha = 3.0;

    std::vector<cl::tape_double> beta(1);
    beta[0] = 2.0;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> meanvariance(2);

    // Start of tape recording.
    // Mark \beta as independent variable.
    Independent(beta);

    boost::math::inverse_gamma_distribution<Real> invGammaDistribution(alpha, beta[0]);

    // Analytical formula for Mean:
    // \frac{\beta}{\alpha-1}
    // \alpha > 1
    meanvariance[0] = boost::math::mean(invGammaDistribution);

    // Analytical formula for Variance:
    // \frac{\beta^2}{(\alpha-1)^2(\alpha-2)}
    // \alpha > 2
    meanvariance[1] = boost::math::variance(invGammaDistribution);

    // End of tape recording.
    cl::tape_function<double> f(beta, meanvariance);

    std::vector<double> sfForward;
    // Start of  differentiation in Forward mode.
    gradForward(f, sfForward, false, false);

    // Create vector of independent variables:
    // x > 0
    std::vector<cl::tape_double> X(1);
    X[0] = 4.0;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> cdf(1);

    // Start of tape recording.
    // Mark X as independent variable.
    Independent(X);

    // Analytical formula for CDF:
    // \frac{\Gamma(\alpha,\beta/x)}{\Gamma(\alpha)}
    // \Gamma() is the gamma fnction
    cdf[0] = boost::math::cdf(invGammaDistribution, X[0]);

    // End of tape recording.
    cl::tape_function<double> f2(X, cdf);

    std::vector<double> sfReverse;
    // Start of differentiation in Reverse mode.
    gradReverse(f2, sfReverse, false, false);

    // Create vector of adjoint derivatives.
    std::vector<double> sfAdjoint;
    sfAdjoint.push_back(sfForward[0]);
    sfAdjoint.push_back(sfForward[1]);
    sfAdjoint.push_back(sfReverse[0]);

    // Create vector of analytical derivatives.
    std::vector<Real> sfAnalytical;

    // Analytical formula for Mean derivative with respect to \beta:
    // \frac{1}{(\alpha - 1}
    sfAnalytical.push_back(1 / (alpha - 1));

    // Analytical formula for Variance derivative with respect to \beta:
    // \frac{2\beta}{(\alpha-1)^2(\alpha-2)}
    sfAnalytical.push_back((2 * beta[0]) / ((alpha - 1)*(alpha - 1)*(alpha - 2)));

    // Analytical formula for PDF:
    // \frac{2^{-\nu/2}}{\Gamma(\nu/2)}x^{-\nu/2-1}e^{-1/(2 x)}
    // \Gamma() is the gamma fnction
    sfAnalytical.push_back(boost::math::pdf(invGammaDistribution, X[0]));

    // Check with explicit analytical formulas.
    result = checkWithAnalyticalFormula(sfAdjoint, sfAnalytical);
    return result;
}

// Inverse gaussian distribution.
bool AdjointDistributionTest::testInverseGaussianDistribution()
{
    BOOST_TEST_MESSAGE("Testing adjoint differentiation for inverse gaussian distribution...");

    bool result = false;

    // Set distribution parameters:
    // Scale: \lambda > 0
    // Mean: \mu > 0

    Real lambda = 3.0;

    std::vector<cl::tape_double> mu(1);
    mu[0] = 2.0;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> meanvariance(2);

    // Start of tape recording.
    // Mark \mu as independent variable.
    Independent(mu);

    boost::math::inverse_gaussian_distribution<Real> invGaussianDistribution(mu[0], lambda);

    // Analytical formula for Mean:
    // \mu
    meanvariance[0] = boost::math::mean(invGaussianDistribution);

    // Analytical formula for Variance:
    // \frac{\mu^3}{\lambda}
    meanvariance[1] = boost::math::variance(invGaussianDistribution);

    // End of tape recording.
    cl::tape_function<double> f(mu, meanvariance);

    std::vector<double> sfForward;
    // Start of  differentiation in Forward mode.
    gradForward(f, sfForward, false, false);

    // Create vector of independent variables:
    // x > 0
    std::vector<cl::tape_double> X(1);
    X[0] = 4.0;

    // Create vector of dependent variables.
    std::vector<cl::tape_double> cdf(1);

    // Start of tape recording.
    // Mark X as independent variable.
    Independent(X);

    // Analytical formula for CDF:
    //  \Phi\left(\sqrt{\frac{\lambda}{x}} \left(\frac{x}{\mu}-1 \right)\right)
    //  +\exp\left(\frac{2 \lambda}{\mu}\right) \Phi\left(-\sqrt{\frac{\lambda}{x}}\left(\frac{x}{\mu}+1 \right)\right)
    //  \Phi() is is the standard normal (standard Gaussian) distribution CDF
    cdf[0] = boost::math::cdf(invGaussianDistribution, X[0]);

    // End of tape recording.
    cl::tape_function<double> f2(X, cdf);

    std::vector<double> sfReverse;
    // Start of differentiation in Reverse mode.
    gradReverse(f2, sfReverse, false, false);

    // Create vector of adjoint derivatives.
    std::vector<double> sfAdjoint;
    sfAdjoint.push_back(sfForward[0]);
    sfAdjoint.push_back(sfForward[1]);
    sfAdjoint.push_back(sfReverse[0]);

    // Create vector of analytical derivatives.
    std::vector<Real> sfAnalytical;

    // Analytical formula for Mean derivative with respect to \mu:
    // 1
    sfAnalytical.push_back(1);

    // Analytical formula for Variance derivative with respect to \mu:
    // \frac{3\mu^2}{\lambda}
    sfAnalytical.push_back(3 * mu[0] * mu[0] / lambda);

    // Analytical formula for PDF:
    //  \left[\frac{\lambda}{2 \pi x^3}\right]^{1/2} \exp{\frac{-\lambda (x-\mu)^2}{2 \mu^2 x}}
    sfAnalytical.push_back(boost::math::pdf(invGaussianDistribution, X[0]));

    // Check with explicit analytical formulas.
    result = checkWithAnalyticalFormula(sfAdjoint, sfAnalytical);
    return result;
}

test_suite* AdjointDistributionTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Adjoint distributions test");
    suite->add(QUANTLIB_TEST_CASE(&AdjointDistributionTest::testBetaDistribution));
    suite->add(QUANTLIB_TEST_CASE(&AdjointDistributionTest::testCauchyDistribution));
    suite->add(QUANTLIB_TEST_CASE(&AdjointDistributionTest::testChiSquaredDistribution));
    suite->add(QUANTLIB_TEST_CASE(&AdjointDistributionTest::testExponentialDistribution));
    suite->add(QUANTLIB_TEST_CASE(&AdjointDistributionTest::testFisherFDistribution));
    suite->add(QUANTLIB_TEST_CASE(&AdjointDistributionTest::testGammaDistribution));
    suite->add(QUANTLIB_TEST_CASE(&AdjointDistributionTest::testLaplaceDistribution));
    suite->add(QUANTLIB_TEST_CASE(&AdjointDistributionTest::testLognormalDistribution));
    suite->add(QUANTLIB_TEST_CASE(&AdjointDistributionTest::testNormalDistribution));
    suite->add(QUANTLIB_TEST_CASE(&AdjointDistributionTest::testParetoDistribution));
    suite->add(QUANTLIB_TEST_CASE(&AdjointDistributionTest::testStudentsTDistribution));
    suite->add(QUANTLIB_TEST_CASE(&AdjointDistributionTest::testInverseChiSquaredDistribution));
    suite->add(QUANTLIB_TEST_CASE(&AdjointDistributionTest::testInverseGammaDistribution));
    suite->add(QUANTLIB_TEST_CASE(&AdjointDistributionTest::testInverseGaussianDistribution));
    return suite;
}


#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_distribution)

BOOST_AUTO_TEST_CASE(testBetaDistribution)
{
    BOOST_CHECK(AdjointDistributionTest::testBetaDistribution());
}

BOOST_AUTO_TEST_CASE(testCauchyDistribution)
{
    BOOST_CHECK(AdjointDistributionTest::testCauchyDistribution());
}

BOOST_AUTO_TEST_CASE(testChiSquaredDistribution)
{
    BOOST_CHECK(AdjointDistributionTest::testChiSquaredDistribution());
}


BOOST_AUTO_TEST_CASE(testExponentialDistribution)
{
    BOOST_CHECK(AdjointDistributionTest::testExponentialDistribution());
}

BOOST_AUTO_TEST_CASE(testFisherFDistribution)
{
    BOOST_CHECK(AdjointDistributionTest::testFisherFDistribution());
}

BOOST_AUTO_TEST_CASE(testGammaDistribution)
{
    BOOST_CHECK(AdjointDistributionTest::testGammaDistribution());
}

BOOST_AUTO_TEST_CASE(testLaplaceDistribution)
{
    BOOST_CHECK(AdjointDistributionTest::testLaplaceDistribution());
}

BOOST_AUTO_TEST_CASE(testLognormalDistribution)
{
    BOOST_CHECK(AdjointDistributionTest::testLognormalDistribution());
}

BOOST_AUTO_TEST_CASE(testNormalDistribution)
{
    BOOST_CHECK(AdjointDistributionTest::testNormalDistribution());
}

BOOST_AUTO_TEST_CASE(testParetoDistribution)
{
    BOOST_CHECK(AdjointDistributionTest::testParetoDistribution());
}

BOOST_AUTO_TEST_CASE(testStudentsTDistribution)
{
    BOOST_CHECK(AdjointDistributionTest::testStudentsTDistribution());
}

BOOST_AUTO_TEST_CASE(testInverseChiSquaredDistribution)
{
    BOOST_CHECK(AdjointDistributionTest::testInverseChiSquaredDistribution());
}

BOOST_AUTO_TEST_CASE(testInverseGammaDistribution)
{
    BOOST_CHECK(AdjointDistributionTest::testInverseGammaDistribution());
}

BOOST_AUTO_TEST_CASE(testInverseGaussianDistribution)
{
    BOOST_CHECK(AdjointDistributionTest::testInverseGaussianDistribution());
}

BOOST_AUTO_TEST_SUITE_END()

#endif
