/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2006 Allen Kuo
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

#include "adjointfraportfoliotest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include <ql/quantlib.hpp>

using namespace QuantLib;
using namespace std;
using namespace boost::unit_test_framework;

namespace
{
    struct CommonVars
    {
        CommonVars() :
        tolerance_(1e-4)
        , h_(1.0e-10)
        {
        }

        double tolerance_;
        double h_;
        double timeTapeRecording_;
        double timeAdjoint_;
        double timeAnalytical_;
        vector<boost::shared_ptr<ForwardRateAgreement>> fraPortfolio_;
        vector<PerformanceTime> performanceTime_;
    };

    struct RateVariation
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "Forward Rate", "Forward Value"
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, RateVariation& v)
        {
                stm << v.rate_
                    << ";" << v.forwardValue_
                    << std::endl;
                return stm;
            }

        Real rate_;
        Real forwardValue_;
    };
}

// This function creates portfolio of Forward Rate Agreements (FRA)
// based on values of forward rates. Each forward rate refers to
// 3 month term FRA quotes (index refers to months to start) are considered.
// Input parameters:
//  - FraRate - vector of forward rate for every FRA in portfolio;
// Function returns vector of FRAs.

vector<boost::shared_ptr<ForwardRateAgreement>> CreateFraPortfolio(vector<Real>& FraRate)
{
    Size n = FraRate.size();
    vector<boost::shared_ptr<ForwardRateAgreement>> fraPortfolio;

    // Set Market Data
    RelinkableHandle<YieldTermStructure> euriborTermStructure;
    boost::shared_ptr<IborIndex> euribor(
        new Euribor3M(euriborTermStructure));

    Date todaysDate = Date::todaysDate();
    Settings::instance().evaluationDate() = todaysDate;

    Calendar calendar = euribor->fixingCalendar();
    Integer fixingDays = euribor->fixingDays();
    Date settlementDate = calendar.advance(todaysDate, fixingDays, Days);

    // Quotes
    // Forward Rate Agreements are quoted in the format AxB, with (A) representing
    // the number of months until the loan is set to begin, and
    // (B) representing the number of months until the loan ends.
    // 3 month term FRA quotes (index refers to monthsToStart) are considered.


    // Forward rate depends of Forward Rate Agreement quote.
    vector<boost::shared_ptr<SimpleQuote> > fraRate(n);
    for (Size i = 0; i < n; ++i)
        fraRate[i] = boost::make_shared<SimpleQuote>(FraRate[i]);

    //Create relinkable handle to rate for Forward Rate Agreement.
    vector<RelinkableHandle<Quote> > quoteHandles(n);
    for (Size i = 0; i < n; ++i)
        quoteHandles[i] = RelinkableHandle<Quote>(fraRate[i]);

    // Rate helpers

    // RateHelpers are built from the above quotes together with
    // other instrument dependant infos.  Quotes are passed in
    // relinkable handles which could be relinked to some other
    // data source later.

    DayCounter fraDayCounter = euribor->dayCounter();
    BusinessDayConvention convention = euribor->businessDayConvention();
    bool endOfMonth = euribor->endOfMonth();
    vector<boost::shared_ptr<RateHelper>> fra(n);
    for (Size i = 0; i < n; ++i)
        fra[i] = boost::make_shared<FraRateHelper>(
        quoteHandles[i], (i + 1), (i + 4), fixingDays, calendar, convention,
        endOfMonth, fraDayCounter);


    // Curve building

    // Any DayCounter would be fine.
    // ActualActual::ISDA ensures that 30 years is 30.0.
    DayCounter termStructureDayCounter = ActualActual(ActualActual::ISDA);

    double tolerance = 1.0e-15;

    // Create Forward Rate Agreement curve.
    boost::shared_ptr<YieldTermStructure> fraTermStructure(
        new PiecewiseYieldCurve<Discount, LogLinear>(
        settlementDate, fra,
        termStructureDayCounter, tolerance));

    // Term structures used for pricing/discounting.
    RelinkableHandle<YieldTermStructure> discountingTermStructure;
    discountingTermStructure.linkTo(fraTermStructure);

    // Construct Forward Rate Agreements.

    // Set position of Forward Rate Agreement as long (Buyer - the one who pays fixed interest).
    Position::Type fraFwdType = Position::Long;
    // Set the notional amount for which FRA is contracted.
    Real fraNotional = 100.0;
    // Set length of the loan.
    const Integer FraTermMonths = 3;

    euriborTermStructure.linkTo(fraTermStructure);

    for (Size i = 0; i < n; i++)
    {
        // Set value (settelment) date (The date on which the notional loan becomes effective;
        // FRA gets cash settled on this date).
        Date fraValueDate = calendar.advance(
            settlementDate, i + 1, Months,
            convention);

        // Set maturity date (The date on which the notional loan matures).
        Date fraMaturityDate = calendar.advance(
            fraValueDate, FraTermMonths, Months,
            convention);

        boost::shared_ptr<ForwardRateAgreement> FRA = boost::make_shared<ForwardRateAgreement>(
            fraValueDate, fraMaturityDate,
            fraFwdType, FraRate[i],
            fraNotional, euribor,
            discountingTermStructure);

        // Add new Forward Rate Agreement to Forward Rate Agreement Portfolio.
        fraPortfolio.push_back(FRA);
    }

    return fraPortfolio;
}

// Use central finite difference for approximation the derivatives
// of foward value of each Forward Rate Agreement on forward rate.
// dy/dx = (y(x+h) - y(x-h))/(2*h)
// Input parameters:
//  - FraRate - vector of forward rate for every FRA in portfolio;
//  - h - step size for finite difference method;
//  - sf_Finite - calculated derivatives.
// Function returns time for calculation derivatives.

double calculateCentralFinDiff(vector<Real>& FraRate, double h, vector<Real>& sf_Finite)
{
    // Start timing for calculating derivatives by central finite difference.
    boost::timer timer;
    Size n = FraRate.size();
    sf_Finite.resize(n);
    vector<boost::shared_ptr<ForwardRateAgreement>> fraPortfolioRight_;
    vector<boost::shared_ptr<ForwardRateAgreement>> fraPortfolioLeft_;
    for (Size i = 0; i < n; i++)
    {
        // Create Forward Rate Agreement Portfolio with shifting rates with a step size of +h.
        FraRate[i] += h;
        fraPortfolioRight_ = CreateFraPortfolio(FraRate);
        // Create Forward Rate Agreement Portfolio with shifting rates with a step size of -h.
        FraRate[i] -= 2 * h;
        fraPortfolioLeft_ = CreateFraPortfolio(FraRate);
        //Evaluate derivatives using central finite difference
        sf_Finite[i] = (fraPortfolioRight_[i]->forwardValue() - fraPortfolioLeft_[i]->forwardValue()) / (2 * h);

    }
    // Return calculated time for calculating derivatives by central finite difference.
    return timer.elapsed();
}

// This function shows dependence of Forward Value of
// Forward Rate Agreement on variation of forward rate.

void testRateVariation()
{
    std::vector<RateVariation> result;
    vector<Real> FraRate(1);
    vector<boost::shared_ptr<ForwardRateAgreement>> fraPortfolio;

    for (Size i = 0; i < 400; i += 4)
    {
        // Change rate on the value 0.00001 on each iteration.
        FraRate[0] = 0.030 + i*0.00001;
        // Create Forward Rate Agreement with new rate.
        fraPortfolio = CreateFraPortfolio(FraRate);
        // Save results in structure.
        result.push_back(RateVariation { FraRate[0], fraPortfolio[0]->forwardValue() });
    }
    // Output results.
     cl::AdjointTestOutput out("AdjointFRAPortfolio//Output"
                               , { { "filename", "Output" }
                               , { "not_clear", "Not" }
                               , { "title", "Forward value dependence on FRA rate" }
                               , { "ylabel", "Forward Value" }
                               , { "xlabel", "FRA Rate" } });
                               out << result;

}

// This example tests adjoint differentiation of forward value
// of Forward Rate Agreement portfolio on initial forward rates.
// Method returns true if adjoint derivatives equals to finite difference derivative; otherwise, it throws exception and returns false.

bool AdjointFRAPortfolioTest::testFRAPortfolio()
{
    BOOST_TEST_MESSAGE("Testing the CppAD differentiation of forward value of FRA portfolio to a change in initial forward rates...");
    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVars var;
    boost::timer timer;
    double timeForward, timeReverse;
    vector<double> sf_Forward, sf_Reverse;
    vector<Real> sf_Finite;
    std::vector<cl::TapeDouble> forwardValuePortfolio(1, 0.0);
    std::vector<cl::TapeDouble> FraRate;
    vector<AdjointTime> performanceAdjointTime_;
    vector<TapeSize> tapeMemory_;
    Size maxSize = 100;

    // Plots streams.
    cl::AdjointTestOutput outPerform("AdjointFRAPortfolio"
                                    , { { "filename", "AdjointPerformance" }
                                    , { "not_clear", "Not" }
                                    , { "title", "Foward value differentiation performance with respect to number of forward rates" }
                                    , { "ylabel", "Time (s)" }
                                    , { "xlabel", "Number of forward rates" }
                                    , { "line_box_width", "-5" } });

    cl::AdjointTestOutput outAdjoint("AdjointFRAPortfolio"
                                    , { { "filename", "Adjoint" }
                                    , { "not_clear", "Not" }
                                    , { "title", "Foward value adjoint differentiation performance with respect to number of forward rates" }
                                    , { "ylabel", "Time (s)" }
                                    , { "xlabel", "Number of forward rates" } });

    cl::AdjointTestOutput outSize("AdjointFRAPortfolio"
                                    , { { "filename", "TapeSize" }
                                    , { "not_clear", "Not" }
                                    , { "title", "Tape size dependence on number of forward rates" }
                                    , { "ylabel", "Memory (MB)" }
                                    , { "xlabel", "Number of forward rates" } });

#ifdef CL_GRAPH_GEN
    Size startSize = 0;
#else
    Size startSize = maxSize - 1;
    // Initialize vector.
    for (Size n = 0; n < startSize; n++)
        FraRate.push_back(0.030 + n*0.0000001);
#endif
    // Change number of Forward Rate Agreements in portfolio.
    for (Size n = startSize; n < maxSize; n++)
    {
        outPerform.log() << "\nNumber of Forward Rate Agreements in portfolio: n = " << n << std::endl;
        // Add new forward rate that refers to 3 month term Forward Rate Agreement quote that begin in n-th month.
        FraRate.push_back(0.030 + n*0.0000001);

        // Start timing of tape recording.
        outPerform.log() << "Start taping : " << currentTime() << std::endl;
        timer.restart();

        // Start taping. Declare forward rates as independent variables.
        Independent(FraRate);

        // Create Portfolio with n Forward Rate Agreement.
        outPerform.log() << "Create Portfolio with n Forward Rate Agreement." << std::endl;
        var.fraPortfolio_ = CreateFraPortfolio(FraRate);

        // Calculate total forward value as sum of forward value of every Forward Rate Agreement in portfolio.
        // forward value = notional amount * (1 + rate per compounding period)^(number of compounding periods)
        outPerform.log() << "Calculate total forward value as sum of forward value of every Forward Rate Agreement in portfolio." << std::endl;
        for (Size i = 0; i < n + 1; i++)
            forwardValuePortfolio[0] += var.fraPortfolio_[i]->forwardValue();

        // End of tape recording. Declare total forward value as dependent variable.
        // Differentiaion will be held with respect to the independent variables vector.
        cl::TapeFunction<double> f(FraRate, forwardValuePortfolio);
        var.timeTapeRecording_ = timer.elapsed();
        outPerform.log() << "End of tape recording. Time for tape recording : " << var.timeTapeRecording_ << std::endl;
        //Store size of tape.
        tapeMemory_.push_back(TapeSize { n + 1, f.Memory() });

        //Start differentiation in Forward mode.
        outPerform.log() << "Start differentiation in Forward mode : " << currentTime() << std::endl;
        timeForward = gradForward(f, sf_Forward, false, false);
        outPerform.log() << "Time for differentiation in Forward mode : " << timeForward  <<" s" << std::endl;

        //Start differentiation in Reverse mode.
        outPerform.log() << "Start differentiation in Reverse mode : " << currentTime() << std::endl;
        timeReverse = gradReverse(f, sf_Reverse, false, false);
        outPerform.log() << "Time for differentiation in Reverse mode : " << timeReverse << " s" << std::endl;
        var.timeAdjoint_ = timeForward < timeReverse ? timeForward : timeReverse;

        //Start differentiation using central finite differences.
        outPerform.log() << "Start differentiation using central finite differences : " << currentTime() << std::endl;
        var.timeAnalytical_ = calculateCentralFinDiff(FraRate, var.h_, sf_Finite);
        outPerform.log() << "Time for differentiation using central finite differences : " << var.timeAnalytical_ << " s" << std::endl;

        // Check derivatives calculated by forward, reverse and finite differences methods.
        outPerform.log() << "Check derivatives calculated by forward, reverse and finite differences methods." << std::endl;
        result = checkWithFiniteDiff(sf_Forward, sf_Reverse, sf_Finite, var.tolerance_, 0.5);

        //Adding new data to the performace result vector.
        var.performanceTime_.push_back(PerformanceTime { var.timeTapeRecording_, var.timeAdjoint_, var.timeAnalytical_, n + 1 });
    }

    //Output differentiation performance results to csv file and create plot.
    outPerform << var.performanceTime_;

    //Output adjoint differentiation performance results to csv file and create plot.
    for (Size n = 0; n < maxSize - startSize; n += 5)
        performanceAdjointTime_.push_back(AdjointTime { var.performanceTime_[n].timeAdjoint_, var.performanceTime_[n].indepVarNumber_ });
    outAdjoint << performanceAdjointTime_;

    //Output tape size dependence to csv file and create plot.
    outSize << tapeMemory_;
#endif

    return result;
}

test_suite* AdjointFRAPortfolioTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("AD FRA Portfolio test");
    suite->add(QUANTLIB_TEST_CASE(&AdjointFRAPortfolioTest::testFRAPortfolio));
    suite->add(QUANTLIB_TEST_CASE(testRateVariation));
    return suite;
}


#if defined CL_ENABLE_BOOST_TEST_ADAPTER && defined CL_OPEN_EXCLUDED

BOOST_AUTO_TEST_SUITE(fra_portfolio)

BOOST_AUTO_TEST_CASE(testFRAPortfolio)
{
    BOOST_CHECK(AdjointFRAPortfolioTest::testFRAPortfolio());
}

BOOST_AUTO_TEST_SUITE_END()

#endif
