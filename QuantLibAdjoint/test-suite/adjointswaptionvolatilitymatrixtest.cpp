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

// Based on swaptionvolatilitymatrix.cpp from Quantlib/test - suite.

#include "adjointswaptionvolatilitymatrixtest.hpp"
#include <test-suite/swaptionvolstructuresutilities.hpp>
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include <ql/indexes/swap/euriborswap.hpp>
#include <ql/instruments/makeswaption.hpp>
#include <ql/pricingengines/swaption/blackswaptionengine.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <string>

#define CL_OPEN_GENERAL_OUPUT

using namespace QuantLib;
using namespace boost::unit_test_framework;

namespace
{
    struct ExpVolVariation
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "Expected Vol", "implied Vol"
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, ExpVolVariation& v)
        {
                stm << v.expVol_
                    << ";" << v.implVol_
                    << std::endl;
                return stm;
            }

        Volatility expVol_;
        Volatility implVol_;
    };
    struct CommonVars
    {
        // Setup.
        CommonVars() :
             swapSize_(1)
            , optionSize_(100)
        {
            conventions_.setConventions();
            atm_.setMarketData();
            Settings::instance().evaluationDate() =
                conventions_.calendar.adjust(Date::todaysDate());
            atmVolMatrix_ = RelinkableHandle<SwaptionVolatilityStructure>(
                boost::shared_ptr<SwaptionVolatilityStructure>(new
                SwaptionVolatilityMatrix(conventions_.calendar,
                conventions_.optionBdc,
                atm_.tenors.options,
                atm_.tenors.swaps,
                atm_.volsHandle,
                conventions_.dayCounter)));
            termStructure_.linkTo(
                boost::shared_ptr<YieldTermStructure>(new
                FlatForward(0, conventions_.calendar,
                0.05, Actual365Fixed())));
        }

        // Global data.
        Date referenceDate_;
        SwaptionMarketConventions conventions_;
        AtmVolatility atm_;
        Size swapSize_;
        Size optionSize_;
        RelinkableHandle<YieldTermStructure> termStructure_;
        RelinkableHandle<SwaptionVolatilityStructure> atmVolMatrix_;
        Real tolerance_;
        std::vector<PerformanceTime> performanceTime_;
        std::vector<Real> expVolRandom_;
        // Cleanup.
        SavedSettings backup;

        // Find implied volatility based on expected volatility.
        // Input parameters:
        // - expVol - vector of expected volatilities;
        // - engine - pricing engine for creating swaption.
        // Function returns vector of implied volatilities.
        std::vector<Volatility> findimpliedVol(boost::shared_ptr<BlackSwaptionEngine> engine, std::vector<Volatility>& expVol)
        {
            std::vector<Volatility> implVolVector;
            Real error, tolerance = 0.1;
            Size n = expVol.size() / swapSize_;
            for (Size j = 0; j < swapSize_; j++)
            {
                boost::shared_ptr<SwapIndex> swapIndex(new
                                                       EuriborSwapIsdaFixA(Period(j + 1, Years), termStructure_));

                for (Size i = 0; i < n; ++i)
                {
                    // Create swaption.
                    Swaption swaption =
                        MakeSwaption(swapIndex, Period(i + 1, Months))
                        .withPricingEngine(engine);
                    // Calculate NPV.
                    Real npv = swaption.NPV();
                    // Calculate implied volatility and  put it to the vector.
                    implVolVector.push_back(swaption.impliedVolatility(npv, termStructure_,
                        expVol[j*n + i] * 0.98, 10e-6,
                        100, 10.0e-7, 4.0, 0.0));
                }
            }
            return implVolVector;
        }

        // Use central finite difference for approximation the derivatives
        // of foward value of each Forward Rate Agreement on forward rate.
        // dy/dx = (y(x+h) - y(x-h))/(2*h)
        // Input parameters:
        //  - engine - pricing engine for swaption
        //  - expVol - vector of expected volatilities;
        //  - sf_Analytical - vector for calculated derivatives.
        //  - h - step size for finite difference method;
        // Function returns time for calculation derivatives.

        double calculateCentralFinDiff(boost::shared_ptr<BlackSwaptionEngine> engine, std::vector<Volatility>& expVol, std::vector<Real>& sf_Analytical, double h)
        {
            Size n = expVol.size();
            sf_Analytical.resize(n*n);
            std::vector<Volatility> implVolVectorRight, implVolVectorLeft;
            implVolVectorRight.resize(n);
            implVolVectorLeft.resize(n);
            // Start timing for calculating derivatives by central finite difference.
            boost::timer timer;
            for (Size i = 0; i < n; i++)
            {
                // Find implied volatility with shifting expected volatility with a step size of +h.
                expVol[i] += h;
                implVolVectorRight = findimpliedVol(engine, expVol);
                // Find implied volatility with shifting expected volatility with a step size of -h.
                expVol[i] -= 2 * h;
                implVolVectorLeft = findimpliedVol(engine, expVol);
                //Evaluate derivatives using central finite difference
                sf_Analytical[i*n + i] = (implVolVectorRight[i] - implVolVectorLeft[i]) / (2 * h);
                expVol[i] += h;
            }
            // Return calculated time for calculating derivatives by central finite difference.
            return timer.elapsed();
        }

        // Initialize vector of expected volatilities.

        void initializeExpVol(std::vector<Real>& expVol)
        {
            for (Size i = 0; i < swapSize_; ++i)
                expVol.push_back(0.005 + i*(0.200 - 0.005) / 100);
        }

        // Function returns true if adjoint derivatives equals to finite difference derivatives
        // and expected and implied volatilities are cohesive for all datasets; otherwise, it throws exception and returns false.

        bool makeCoherenceTest(
            const std::string& description,
            const boost::shared_ptr<SwaptionVolatilityDiscrete>& vol)
        {
            // Create pricing engine
            boost::shared_ptr<BlackSwaptionEngine> engine(new
                                                          BlackSwaptionEngine(termStructure_,
                                                          Handle<SwaptionVolatilityStructure>(vol)));
            boost::timer timer;
            double tol = 1e-4;
            double h = 1e-3;
            double timeTapeRecording;
            double timeAdjoint;
            double timeAnalytical;
            performanceTime_.clear();
            std::vector<AdjointTime> performanceAdjointTime;
            std::vector<TapeSize> tapeMemory_;
            std::vector<cl::TapeDouble> implVolVector;
            std::vector<cl::TapeDouble> expVol;
            Size maxSize = 100;
            bool result = false;

            // Plots streams.
            cl::AdjointTestOutput outPerform(
                "AdjointSwaptionVolatilityMatrix//" + description
                , { { "filename", "AdjointPerformance" }
                , { "not_clear", "Not" }
                , { "title", "Implied volatility differentiation performance with respect to number of volatilities for " + description }
                , { "ylabel", "Time (s)" }
                , { "xlabel", "Number of expected volatilities" }
                , { "line_box_width", "-5" } });

            cl::AdjointTestOutput outAdjoint(
                "AdjointSwaptionVolatilityMatrix//" + description
                , { { "filename", "Adjoint" }
                , { "not_clear", "Not" }
                , { "title", "Implied volatility adjoint performance  with respect to number of volatilities for " + description }
                , { "ylabel", "Time (s)" }
                , { "xlabel", "Number of expected volatilities" } });

            cl::AdjointTestOutput outSize(
                "AdjointSwaptionVolatilityMatrix//" + description
                , { { "filename", "TapeSize" }
                , { "not_clear", "Not" }
                , { "title", "Tape size dependence on number of number of expected volatilities" }
                , { "ylabel", "Memory (MB)" }
                , { "xlabel", "Number of expected volatilities" } });

#ifdef CL_GRAPH_GEN
            Size startSize = 0;
#else
            Size startSize = maxSize - 1;
            for (Size n = 0; n < startSize; n++)
                initializeExpVol(expVol);
#endif
            //Change number of expected volatilities.
            for (Size n = startSize; n < maxSize; n++)
            {
                // Add expected volatilities.
                outPerform.log() << "\nNumber of expected volatilities: n = " << n << std::endl;
                initializeExpVol(expVol);
                // Start timing of tape recording.
                outPerform.log() << "Start taping : " << currentTime() << std::endl;
                timer.restart();

                // Start taping. Declare expected volatilities as independent variables.
                Independent(expVol);

                // Find implied volatilities.
                outPerform.log() << "Find implied volatilities." << std::endl;
                implVolVector = findimpliedVol(engine, expVol);

                // End of tape recording. Declare vector of implied volatilities as dependent variables.
                // Differentiaion will be held with respect to the independent variables vector.
                cl::TapeFunction<double> f(expVol, implVolVector);
                timeTapeRecording = timer.elapsed();
                outPerform.log() << "End of tape recording. Time for tape recording : " << timeTapeRecording << std::endl;
                // Store size of tape.
                 tapeMemory_.push_back(TapeSize { n + 1, f.Memory() });

                std::vector<double> sf_Forward, sf_Reverse;
                // Start differentiation in Forward mode.
                outPerform.log() << "Start differentiation in Forward mode : " << currentTime() << std::endl;
                double timeForward = gradForward(f, sf_Forward, false, false);
                outPerform.log() << "Time for differentiation in Forward mode : " << timeForward << " s" << std::endl;

                // Start differentiation in Reverse mode.
                outPerform.log() <<"Start differentiation in Reverse mode : " << currentTime() << std::endl;
                double timeReverse = gradReverse(f, sf_Reverse, false, false);
                outPerform.log() << "Time for differentiation in Reverse mode : " << timeReverse << " s" << std::endl;
                timeAdjoint = timeForward < timeReverse ? timeForward : timeReverse;

                // Start differentiation using central finite differences.
                outPerform.log() << "Start differentiation using central finite differences : " << currentTime() << std::endl;
                std::vector<Real> sf_Analytical;
                timeAnalytical = calculateCentralFinDiff(engine, expVol, sf_Analytical, h);
                outPerform.log() << "Time for differentiation using central finite differences : " << timeAnalytical << " s" << std::endl;

                // Check derivatives calculated by forward, reverse and finite differences methods.
                outPerform.log() << "Check derivatives calculated by forward, reverse and finite differences methods." << std::endl;
                result = checkWithFiniteDiff(sf_Forward, sf_Reverse, sf_Analytical, tol, tol);

                // Adding new data to the performace result vector.
                performanceTime_.push_back(PerformanceTime { timeTapeRecording, timeAdjoint, timeAnalytical, n + 1 });
            }

            //Output differentiation performance results to csv file and create plot.          
            outPerform << performanceTime_;

            //Output adjoint differentiation performance results to csv file and create plot.
            for (Size n = 0; n < maxSize - startSize; n += 5)
                performanceAdjointTime.push_back(AdjointTime { performanceTime_[n].timeAdjoint_, performanceTime_[n].indepVarNumber_ });

           
            outAdjoint << performanceAdjointTime;

            //Output tape size dependence to csv file and create plot.
            outSize << tapeMemory_;
            return result;
        }

        // This function shows dependence of implied volatilities on variation of expected volatilities.
        // Input parameters:
        // - description - description about used type of market data;
        // - vol - object with swaption volatility matrix.
        void expVolVariance(const std::string& description,
                            const boost::shared_ptr<SwaptionVolatilityDiscrete>& vol)
        {
            std::vector<Volatility> expVol = { 0.005 };
            std::vector<Volatility> implVolVector;
            // Create pricing engine.
            boost::shared_ptr<BlackSwaptionEngine> engine(new
                                                          BlackSwaptionEngine(termStructure_,
                                                          Handle<SwaptionVolatilityStructure>(vol)));
            std::vector<ExpVolVariation> result;
            // Variate value of expected volatility.
            for (Size i = 0; i < 100; i++)
            {
                expVol[0] += (0.200 - 0.05) / 100;
                implVolVector = findimpliedVol(engine, expVol);
                result.push_back(ExpVolVariation { expVol[0], implVolVector[0] });
            }
            // Output results.
            cl::AdjointTestOutput output("AdjointSwaptionVolatilityMatrix//"+ description +"//Output", { { "filename", description }, { "not_clear", "Not" }, { "title", "implied volatility dependence on expected volatility " + description }, { "ylabel", "implied Volatility" }, { "xlabel", "Expected Volatility" } });
            output << result;
        }
    };
}

// This method tests adjoint differentiation of implied volatilities on expected volatility
// of set of swaptions, that represented in matrix form. 

// Swaption volatility matrix provides the at - the - money volatility for a given
// swaption by interpolating a volatility matrix whose elements
// are the market volatilities of a set of swaptions given
// option reference date and swap lengths.

// The volatility matrix M must be defined so that :
//  -the number of rows equals the number of option dates;
//  -the number of columns equals the number of swap tenors;
//  -M[i][j] contains the volatility corresponding
// to the i-th option and j -th tenor.

// TestSwaptionVolMatrixCoherence uses four different datasets that form swaption matrix : 
// floating reference date and floating market data, fixed reference date and floating market data,
// floating reference date and fixed market data, fixed reference date and fixed market data.
// Method returns true if adjoint derivatives equals to finite difference derivatives
// and expected and implied volatilities are cohesive for all datasets; otherwise, it throws exception and returns false.

bool AdjointSwaptionVolatilityMatrixTest::testSwaptionVolMatrixCoherence()
{

    BOOST_TEST_MESSAGE("Testing swaption volatility matrix...");
    bool result = false;

#ifdef CL_TAPE_CPPAD

    CommonVars vars;

    boost::shared_ptr<SwaptionVolatilityMatrix> vol;
    std::string description;


    // Create swaptions with floating reference date, floating market data.
    description = "swaption volatility matrix with floating reference date and floating market data";
    vol = boost::shared_ptr<SwaptionVolatilityMatrix>(new
                                                      SwaptionVolatilityMatrix(vars.conventions_.calendar,
                                                      vars.conventions_.optionBdc,
                                                      vars.atm_.tenors.options,
                                                      vars.atm_.tenors.swaps,
                                                      vars.atm_.volsHandle,
                                                      vars.conventions_.dayCounter));
    result = vars.makeCoherenceTest(description, vol);

    // Create swaptions with fixed reference date, floating market data.
    description = "swaption volatility matrix with fixed reference date and floating market data";
    vol = boost::shared_ptr<SwaptionVolatilityMatrix>(new
                                                      SwaptionVolatilityMatrix(Settings::instance().evaluationDate(),
                                                      vars.conventions_.calendar,
                                                      vars.conventions_.optionBdc,
                                                      vars.atm_.tenors.options,
                                                      vars.atm_.tenors.swaps,
                                                      vars.atm_.volsHandle,
                                                      vars.conventions_.dayCounter));
    result &= vars.makeCoherenceTest(description, vol);

    // Create swaptions with floating reference date, fixed market data.
    description = "swaption volatility matrix with floating reference date and fixed market data";
    vol = boost::shared_ptr<SwaptionVolatilityMatrix>(new
                                                      SwaptionVolatilityMatrix(vars.conventions_.calendar,
                                                      vars.conventions_.optionBdc,
                                                      vars.atm_.tenors.options,
                                                      vars.atm_.tenors.swaps,
                                                      vars.atm_.volsHandle,
                                                      vars.conventions_.dayCounter));
    result &= vars.makeCoherenceTest(description, vol);

    // Create swaptions with fixed reference date, fixed market data.
    description = "swaption volatility matrix with fixed reference date and fixed market data";
    vol = boost::shared_ptr<SwaptionVolatilityMatrix>(new
                                                      SwaptionVolatilityMatrix(Settings::instance().evaluationDate(),
                                                      vars.conventions_.calendar,
                                                      vars.conventions_.optionBdc,
                                                      vars.atm_.tenors.options,
                                                      vars.atm_.tenors.swaps,
                                                      vars.atm_.volsHandle,
                                                      vars.conventions_.dayCounter));
    result &= vars.makeCoherenceTest(description, vol);
#endif
    return result;
}

// This function shows dependence of implied volatilities on variation of expected volatility.
// This function uses four different datasets : floating reference date and floating market data,
// fixed reference date and floating market data, floating reference date and fixed market data, 
// fixed reference date and fixed market data.

void outputSwaptionVolMatrixCoherence()
{

    BOOST_TEST_MESSAGE("Testing swaption volatility matrix...");
    CommonVars vars;

    boost::shared_ptr<SwaptionVolatilityMatrix> vol;
    std::string description;


    // Create swaptions with floating reference date, floating market data.
    description = "swaption volatility matrix with floating reference date and floating market data";
    vol = boost::shared_ptr<SwaptionVolatilityMatrix>(new
                                                      SwaptionVolatilityMatrix(vars.conventions_.calendar,
                                                      vars.conventions_.optionBdc,
                                                      vars.atm_.tenors.options,
                                                      vars.atm_.tenors.swaps,
                                                      vars.atm_.volsHandle,
                                                      vars.conventions_.dayCounter));
    // Test dependence of implied volatilities on variation of expected volatility.
    vars.expVolVariance(description, vol);

    // Create swaptions with fixed reference date, floating market data.
    description = "swaption volatility matrix with fixed reference date and floating market data";
    vol = boost::shared_ptr<SwaptionVolatilityMatrix>(new
                                                      SwaptionVolatilityMatrix(Settings::instance().evaluationDate(),
                                                      vars.conventions_.calendar,
                                                      vars.conventions_.optionBdc,
                                                      vars.atm_.tenors.options,
                                                      vars.atm_.tenors.swaps,
                                                      vars.atm_.volsHandle,
                                                      vars.conventions_.dayCounter));
    // Test dependence of implied volatilities on variation of expected volatility.
    vars.expVolVariance(description, vol);

    // Create swaptions with floating reference date, fixed market data.
    description = "swaption volatility matrix with floating reference date and fixed market data";
    vol = boost::shared_ptr<SwaptionVolatilityMatrix>(new
                                                      SwaptionVolatilityMatrix(vars.conventions_.calendar,
                                                      vars.conventions_.optionBdc,
                                                      vars.atm_.tenors.options,
                                                      vars.atm_.tenors.swaps,
                                                      vars.atm_.volsHandle,
                                                      vars.conventions_.dayCounter));
    // Test dependence of implied volatilities on variation of expected volatility.
    vars.expVolVariance(description, vol);

    // Create swaptions with fixed reference date, fixed market data.
    description = "swaption volatility matrix with fixed reference date and fixed market data";
    vol = boost::shared_ptr<SwaptionVolatilityMatrix>(new
                                                      SwaptionVolatilityMatrix(Settings::instance().evaluationDate(),
                                                      vars.conventions_.calendar,
                                                      vars.conventions_.optionBdc,
                                                      vars.atm_.tenors.options,
                                                      vars.atm_.tenors.swaps,
                                                      vars.atm_.volsHandle,
                                                      vars.conventions_.dayCounter));
    // Test dependence of implied volatilities on variation of expected volatility.
    vars.expVolVariance(description, vol);
}

test_suite* AdjointSwaptionVolatilityMatrixTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Swaption Volatility Matrix tests");

    suite->add(QUANTLIB_TEST_CASE(
        &AdjointSwaptionVolatilityMatrixTest::testSwaptionVolMatrixCoherence));
    suite->add(QUANTLIB_TEST_CASE(outputSwaptionVolMatrixCoherence));


    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_default_curves)

BOOST_AUTO_TEST_CASE(testSwaptionVolMatrixCoherence)
{
    BOOST_CHECK(AdjointSwaptionVolatilityMatrixTest::testSwaptionVolMatrixCoherence());
}

BOOST_AUTO_TEST_SUITE_END()

#endif