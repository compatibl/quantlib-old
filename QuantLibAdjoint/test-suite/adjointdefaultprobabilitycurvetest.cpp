/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2008, 2009 StatPro Italia srl
Copyright (C) 2009 Ferdinando Ametrano
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

//based on defaultirobabilitycurve.cpp from test-suite

#include "adjointdefaultprobabilitycurvetest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include <ql/termstructures/credit/piecewisedefaultcurve.hpp>
#include <ql/termstructures/credit/defaultprobabilityhelpers.hpp>
#include <ql/termstructures/credit/flathazardrate.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/instruments/creditdefaultswap.hpp>
#include <ql/pricingengines/credit/midpointcdsengine.hpp>
#include <ql/math/interpolations/linearinterpolation.hpp>
#include <ql/math/interpolations/backwardflatinterpolation.hpp>
#include <ql/math/interpolations/loginterpolation.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/utilities/dataformatters.hpp>
#include <iomanip>
#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/daycounters/thirty360.hpp>

using namespace QuantLib;
using namespace std;
using namespace boost::unit_test_framework;

struct CommonVarPC
{
    CommonVarPC() :
          dayCounter_(Actual360())
        , calendar_(TARGET())
        , today_(Settings::instance().evaluationDate())
        , startDate_(today_)
        , endDate_(startDate_)
        , frequency_(Quarterly)
        , convention_(Following)
        , rule_(DateGeneration::TwentiethIMM)
        , settlementDays_(1)
        , recoveryRate_(0.4)
        , fixedRate_(0.05)
        , tolerance_(1.0e-10)
        , h_(1e-4)
    {
        discountCurve.linkTo(boost::shared_ptr<YieldTermStructure>(
            new FlatForward(today_, 0.06, Actual360())));
    }

    DayCounter dayCounter_;
    Calendar calendar_;

    Frequency frequency_;
    BusinessDayConvention convention_;
    DateGeneration::Rule rule_;
    Integer settlementDays_;
    Real recoveryRate_;
    Real fixedRate_;
    RelinkableHandle<YieldTermStructure> discountCurve;

    Date today_;
    Date startDate_;
    Date endDate_;
    double tolerance_;
    double timeTapeRecording_;
    double timeAdjoint_;
    double timeAnalytical_;
    vector<PerformanceTime> performanceTime_;
    vector<TapeSize> tapeMemory_;

    boost::timer timer;
    double h_;
};

bool AdjointDefaultProbabilityCurveTest::testFlatHazardRate()
{
    BOOST_TEST_MESSAGE("Testing flat hazard rate...");
    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVarPC vars;
    cl::AdjointTestOutput outPerform("AdjointDefaultProbabilityCurve//FlatHazardRate");

    Size n = 20;
    //Running Adjoint Differentiation with respect to hazard rate
    std::vector<cl::TapeDouble> hazardRate(1);
    hazardRate[0] = 0.0100;

    std::vector<cl::TapeDouble> computedProbability(n);
   
    //Running Adjoint Differentiation with respect to forward price
    outPerform.log() << "\nNumber of forward prices : k = " << 1 << std::endl;

    // Start timing of tape recording.
    outPerform.log() << "Start taping : " << currentTime() << std::endl;
    vars.timer.restart();

    // Start taping. Declare forward rates as independent variables.
    Independent(hazardRate);

    // Compute probability.
    outPerform.log() << "Compute probability." << std::endl;
    Handle<Quote> hazardRateQuote = Handle<Quote>(
        boost::shared_ptr<Quote>(new SimpleQuote(hazardRate[0])));

    FlatHazardRate flatHazardRate(vars.today_, hazardRateQuote, vars.dayCounter_);

    for (Size i = 0; i < n; i++)
    {
        vars.endDate_ = vars.calendar_.advance(vars.endDate_, 1, Years);
        Time t = vars.dayCounter_.yearFraction(vars.startDate_, vars.endDate_);
        Probability probability = 1.0 - std::exp(-hazardRate[0] * t);
        computedProbability[i] = flatHazardRate.defaultProbability(t);

        if (std::fabs(probability - computedProbability[i]) > vars.tolerance_)
            BOOST_ERROR(
            "Failed to reproduce probability for flat hazard rate\n"
            << std::setprecision(10)
            << "    calculated probability: " << computedProbability[i] << "\n"
            << "    expected probability:   " << probability);
    }

    cl::TapeFunction<double> f(hazardRate, computedProbability);
    //End of tape recording. Differentiaion wiil be held with respect to  the independent variables vector
    vars.timeTapeRecording_ = vars.timer.elapsed();
    //End of tape recording. Differentiaion wiil be held with respect to  the independent variables vector
    vars.timeTapeRecording_ = vars.timer.elapsed();
    outPerform.log() << "End of tape recording. Time for tape recording : " << vars.timeTapeRecording_ << std::endl;
    //Store size of tape.
    vars.tapeMemory_.push_back(TapeSize { 1, f.Memory() });

    //Start differentiation in Forward mode
    vector<double> sf_Forward;
    outPerform.log() << "Start differentiation in Forward mode : " << currentTime() << std::endl;
    vars.timeAdjoint_ = gradForward(f, sf_Forward, false, false);
    outPerform.log() << "Time for differentiation in Forward mode : " << vars.timeAdjoint_ << " s" << std::endl;

    //Start differentiation using central finite differences.
    outPerform.log() << "Start differentiation using central finite differences : " << currentTime() << std::endl;
    vector<Real> sf_Analytical(n);
    vars.endDate_ = vars.startDate_;
    for (Size i = 0; i < n; i++)
    {
        vars.endDate_ = vars.calendar_.advance(vars.endDate_, 1, Years);
        Time t = vars.dayCounter_.yearFraction(vars.startDate_, vars.endDate_);
        sf_Analytical[i] = t*std::exp(-hazardRate[0] * t);
    }
    outPerform.log() << "Time for differentiation using central finite differences : " << vars.timeAnalytical_ << " s" << std::endl;

    // Check derivatives calculated by forward, reverse and finite differences methods.
    outPerform.log() << "Check derivatives calculated by forward and finite differences methods." << std::endl;
    result = checkWithFiniteDiff(sf_Forward, sf_Analytical, vars.tolerance_, vars.tolerance_);

    //Adding new data to the performace result vector
    vars.performanceTime_.push_back(PerformanceTime { vars.timeTapeRecording_, vars.timeAdjoint_, vars.timeAnalytical_, 1 });
#endif
    return result;
}

bool AdjointDefaultProbabilityCurveTest::testFlatHazardRateTime()
{
    BOOST_TEST_MESSAGE("Testing flat hazard rate...");
    bool result = false;

#ifdef CL_TAPE_CPPAD
    CommonVarPC vars;
    cl::AdjointTestOutput outPerform("AdjointDefaultProbabilityCurve//FlatHazardRate");

    Real hazardRate = 0.0100;
    Handle<Quote> hazardRateQuote = Handle<Quote>(
        boost::shared_ptr<Quote>(new SimpleQuote(hazardRate)));

    //Running Adjoint Differentiation with respect to various number of time
    std::vector<cl::TapeDouble> t;
    for (Size n = 0; n < 200; n++)
    {
        outPerform.log() << "\nNumber of time periods : k = " << n + 1 << std::endl;
        // Add new time.
        vars.endDate_ = vars.calendar_.advance(vars.endDate_, 1, Months);
        t.push_back(vars.dayCounter_.yearFraction(vars.startDate_, vars.endDate_));

        std::vector<cl::TapeDouble> computedProbability(n + 1);

        // Start timing of tape recording.
        outPerform.log() << "Start taping : " << currentTime() << std::endl;
        vars.timer.restart();

        // Start taping. Declare forward rates as independent variables.
        Independent(t);

        // Compute probability.
        outPerform.log() << "Compute probability." << std::endl;
        FlatHazardRate flatHazardRate(vars.today_, hazardRateQuote, vars.dayCounter_);

        for (Size i = 0; i < n + 1; i++)
        {
            computedProbability[i] = flatHazardRate.defaultProbability(t[i]);
        }

        cl::TapeFunction<double> f(t, computedProbability);
        //End of tape recording. Differentiaion wiil be held with respect to  the independent variables vector
        vars.timeTapeRecording_ = vars.timer.elapsed();
        //End of tape recording. Differentiaion wiil be held with respect to  the independent variables vector
        vars.timeTapeRecording_ = vars.timer.elapsed();
        outPerform.log() << "End of tape recording. Time for tape recording : " << vars.timeTapeRecording_ << std::endl;
        //Store size of tape.
        vars.tapeMemory_.push_back(TapeSize { n + 1, f.Memory() });

        vector<double> sf_Forward, sf_Reverse;

        //Start differentiation in Forward mode.
        outPerform.log() << "Start differentiation in Forward mode : " << currentTime() << std::endl;
        double timeForward = gradForward(f, sf_Forward, false, false);
        outPerform.log() << "Time for differentiation in Forward mode : " << timeForward << " s" << std::endl;

        //Start differentiation in Reverse mode.
        outPerform.log() << "Start differentiation in Forward mode : " << currentTime() << std::endl;
        double timeReverse = gradReverse(f, sf_Reverse, false, false);
        outPerform.log() << "Time for differentiation in Forward mode : " << timeReverse << " s" << std::endl;;

        vars.timeAdjoint_ = timeForward < timeReverse ? timeForward : timeReverse;

        //Start differentiation using analytical formula.
        outPerform.log() << "Start differentiation using analytical formula : " << currentTime() << std::endl;
        vars.timer.restart();
        vector<Real> sf_Analytical((n + 1)*(n + 1), 0);

        for (Size i = 0; i < n + 1; i++)
        {
            for (Size j = 0; j < n + 1; j++)
            {
                if (i == j)
                    sf_Analytical[i*(n + 1) + j] = (flatHazardRate.defaultProbability(t[j] + vars.h_) - flatHazardRate.defaultProbability(t[j] - vars.h_)) / (2 * vars.h_);
            }
        }

        vars.timeAnalytical_ = vars.timer.elapsed();
        outPerform.log() << "Time for differentiation using central finite differences : " << vars.timeAnalytical_ << " s" << std::endl;

        // Check derivatives calculated by forward, reverse and finite differences methods.
        outPerform.log() << "Check derivatives calculated by forward, reverse and finite differences methods." << std::endl;
        double tol = 1e-4;

        //Check results
        result = checkWithFiniteDiff(sf_Forward, sf_Reverse, sf_Analytical, tol, tol);

        //Adding new data to the performace result vector
        vars.performanceTime_.push_back(PerformanceTime { vars.timeTapeRecording_, vars.timeAdjoint_, vars.timeAnalytical_, n + 1 });
  }

#endif
    return result;
}

namespace
{
    template <class T, class I>
    Real computeFairRateFromSpread(CommonVarPC& vars, Real quote, Integer n)
    {
        std::vector<boost::shared_ptr<DefaultProbabilityHelper> > helpers;
        vars.settlementDays_ = 1;
        helpers.push_back(
            boost::shared_ptr<DefaultProbabilityHelper>(
            new SpreadCdsHelper(quote, Period(n, Months),
            vars.settlementDays_, vars.calendar_,
            vars.frequency_, vars.convention_, vars.rule_,
            vars.dayCounter_, vars.recoveryRate_,
            vars.discountCurve)));

        RelinkableHandle<DefaultProbabilityTermStructure> piecewiseCurve;
        piecewiseCurve.linkTo(
            boost::shared_ptr<DefaultProbabilityTermStructure>(
            new PiecewiseDefaultCurve<T, I>(vars.today_, helpers,
            Thirty360())));

        Real notional = 1.0;

        // ensure apple-to-apple comparison
        SavedSettings backup;
        Settings::instance().includeTodaysCashFlows() = true;


        Date protectionStart = vars.today_ + vars.settlementDays_;
        vars.startDate_ = vars.calendar_.adjust(protectionStart, vars.convention_);
        vars.endDate_ = vars.today_ + n * Months;

        Schedule schedule(vars.startDate_, vars.endDate_, Period(vars.frequency_), vars.calendar_,
                          vars.convention_, Unadjusted, vars.rule_, false);

        CreditDefaultSwap cds(Protection::Buyer, notional, quote,
                              schedule, vars.convention_, vars.dayCounter_,
                              true, true, protectionStart);
        cds.setPricingEngine(boost::shared_ptr<PricingEngine>(
            new MidPointCdsEngine(piecewiseCurve, vars.recoveryRate_,
            vars.discountCurve)));

        return cds.fairSpread();
    }

    template <class T, class I>
    Real computeFairRateFromUpfront(CommonVarPC& vars, Real quote, Integer n)
    {
        std::vector<boost::shared_ptr<DefaultProbabilityHelper> > helpers;
        vars.settlementDays_ = 0;

        helpers.push_back(
            boost::shared_ptr<DefaultProbabilityHelper>(
            new UpfrontCdsHelper(quote, vars.fixedRate_, Period(n, Months),
            vars.settlementDays_, vars.calendar_,
            vars.frequency_, vars.convention_, vars.rule_,
            vars.dayCounter_, vars.recoveryRate_,
            vars.discountCurve)));

        RelinkableHandle<DefaultProbabilityTermStructure> piecewiseCurve;
        piecewiseCurve.linkTo(
            boost::shared_ptr<DefaultProbabilityTermStructure>(
            new PiecewiseDefaultCurve<T, I>(vars.today_, helpers,
            Thirty360())));

        Real notional = 1.0;

        // ensure apple-to-apple comparison
        SavedSettings backup;
        Settings::instance().includeTodaysCashFlows() = true;


        Date protectionStart = vars.today_ + vars.settlementDays_;
        vars.startDate_ = vars.calendar_.adjust(protectionStart, vars.convention_);
        vars.endDate_ = vars.today_ + n * Months;

        Schedule schedule(vars.startDate_, vars.endDate_, Period(vars.frequency_), vars.calendar_,
                          vars.convention_, Unadjusted, vars.rule_, false);

        CreditDefaultSwap cds(Protection::Buyer, notional, quote, vars.fixedRate_,
                              schedule, vars.convention_, vars.dayCounter_,
                              true, true, protectionStart);
        cds.setPricingEngine(boost::shared_ptr<PricingEngine>(
            new MidPointCdsEngine(piecewiseCurve, vars.recoveryRate_,
            vars.discountCurve, true)));

        return cds.fairUpfront();
    }

    template <class T, class I>
    bool testBootstrapFromSpread()
    {
        CommonVarPC vars;

        std::vector<Integer> n;
        std::vector<cl::TapeDouble> quote;
        std::vector<cl::TapeDouble> computedRate;
        bool result = false;
        for (Size k = 0; k < 20; k++)
        {
            computedRate.resize(k + 1);
            n.push_back(k + 1);
            quote.push_back((k + 1)*0.001);
            vars.timer.restart();

            Independent(quote);

            for (Size i = 0; i < (k + 1); i++)
            {
                Rate inputRate = quote[i];
                computedRate[i] = computeFairRateFromSpread<T, I>(vars, quote[i], n[i]);

                if (std::fabs(inputRate - computedRate[i]) > vars.tolerance_)
                    BOOST_ERROR(
                    "\nFailed to reproduce fair spread for " << n[i] <<
                    "Y credit-default swaps\n"
                    << std::setprecision(10)
                    << "    computed rate: " << io::rate(computedRate[i]) << "\n"
                    << "    input rate:    " << io::rate(inputRate));
            }

            cl::TapeFunction<double> f(quote, computedRate);
            vars.timeTapeRecording_ = vars.timer.elapsed();
            vector<double> sf_Forward, sf_Reverse;

            //Start differentiation in Forward mode
            double timeForward = gradForward(f, sf_Forward, false, false);

            //Start differentiation in Reverse mode
            double timeReverse = gradReverse(f, sf_Reverse, false, false);
            vars.timeAdjoint_ = timeForward < timeReverse ? timeForward : timeReverse;

            //Central finite difference scheme
            vector<Real> sf_Analytical((k + 1)*(k + 1), 0.0);
            Real rightValue, leftValue;
            vars.timer.restart();

            for (Size i = 0; i < (k + 1); i++)
            {
                rightValue = computeFairRateFromSpread<T, I>(vars, quote[i] + vars.h_, n[i]);
                leftValue = computeFairRateFromSpread<T, I>(vars, quote[i] - vars.h_, n[i]);
                sf_Analytical[i*(k + 1) + i] = (rightValue - leftValue) / (2 * vars.h_);
            }
            vars.timeAnalytical_ = vars.timer.elapsed();

            //Check results
            double tol = 1e-4;
            result &= checkWithFiniteDiff(sf_Forward, sf_Reverse, sf_Analytical, tol, tol);
      }
        return result;
    }


    template <class T, class I>
    bool testBootstrapFromUpfront()
    {
        CommonVarPC vars;

        std::vector<Integer> n;
        std::vector<cl::TapeDouble> quote;
        std::vector<cl::TapeDouble> computedRate;
        bool result = false;
        for (Size k = 0; k < 20; k++)
        {
            computedRate.resize(k + 1);
            n.push_back(k + 2);
            quote.push_back((k + 1)*0.001);
            vars.timer.restart();

            Independent(quote);

            for (Size i = 0; i < (k + 1); i++)
            {
                Rate inputRate = quote[i];
                computedRate[i] = computeFairRateFromUpfront<T, I>(vars, quote[i], n[i]);

                if (std::fabs(inputRate - computedRate[i]) > vars.tolerance_)
                    BOOST_ERROR(
                    "\nFailed to reproduce fair spread for " << n[i] <<
                    "Y credit-default swaps\n"
                    << std::setprecision(10)
                    << "    computed rate: " << io::rate(computedRate[i]) << "\n"
                    << "    input rate:    " << io::rate(inputRate));
            }

            cl::TapeFunction<double> f(quote, computedRate);
            vars.timeTapeRecording_ = vars.timer.elapsed();
            vector<double> sf_Forward, sf_Reverse;

            //Start differentiation in Forward mode
            double timeForward = gradForward(f, sf_Forward, false, false);

            //Start differentiation in Reverse mode
            double timeReverse = gradReverse(f, sf_Reverse, false, false);
            vars.timeAdjoint_ = timeForward < timeReverse ? timeForward : timeReverse;

            //Central finite difference scheme
            vector<Real> sf_Analytical((k + 1)*(k + 1), 0.0);
            Real rightValue, leftValue;
            vars.timer.restart();

            for (Size i = 0; i < (k + 1); i++)
            {
                rightValue = computeFairRateFromUpfront<T, I>(vars, quote[i] + vars.h_, n[i]);
                leftValue = computeFairRateFromUpfront<T, I>(vars, quote[i] - vars.h_, n[i]);
                sf_Analytical[i*(k + 1) + i] = (rightValue - leftValue) / (2 * vars.h_);
            }
            vars.timeAnalytical_ = vars.timer.elapsed();

            //Check results
            double tol = 1e-4;
            result &= checkWithFiniteDiff(sf_Forward, sf_Reverse, sf_Analytical, tol, tol);
        }
        return result;
    }

}

bool AdjointDefaultProbabilityCurveTest::testFlatHazardConsistency()
{
    BOOST_TEST_MESSAGE("Testing piecewise-flat hazard-rate consistency...");
    bool result = false;

#ifdef CL_TAPE_CPPAD

    result = testBootstrapFromSpread<HazardRate, BackwardFlat>();
    result = testBootstrapFromUpfront<HazardRate, BackwardFlat>();
#endif
    return result;
}

bool AdjointDefaultProbabilityCurveTest::testFlatDensityConsistency()
{
    BOOST_TEST_MESSAGE("Testing piecewise-flat default-density consistency...");
    bool result = false;

#ifdef CL_TAPE_CPPAD

    result = testBootstrapFromSpread<DefaultDensity, BackwardFlat>();
    result &= testBootstrapFromUpfront<DefaultDensity, BackwardFlat>();
#endif
    return result;
}

bool AdjointDefaultProbabilityCurveTest::testLinearDensityConsistency()
{
    BOOST_TEST_MESSAGE("Testing piecewise-linear default-density consistency...");
    bool result = false;

#ifdef CL_TAPE_CPPAD

    result = testBootstrapFromSpread<DefaultDensity, Linear>();
    result &= testBootstrapFromUpfront<DefaultDensity, Linear>();
#endif
    return result;
}

bool AdjointDefaultProbabilityCurveTest::testLogLinearSurvivalConsistency()
{
    BOOST_TEST_MESSAGE("Testing log-linear survival-probability consistency...");
    bool result = false;

#ifdef CL_TAPE_CPPAD

    result = testBootstrapFromSpread<SurvivalProbability, LogLinear>();
    result &= testBootstrapFromUpfront<SurvivalProbability, LogLinear>();
#endif
    return result;
}


test_suite* AdjointDefaultProbabilityCurveTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Default-probability curve tests");
    suite->add(QUANTLIB_TEST_CASE(
    &AdjointDefaultProbabilityCurveTest::testFlatHazardRate));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointDefaultProbabilityCurveTest::testFlatHazardRateTime));
   suite->add(QUANTLIB_TEST_CASE(
         &AdjointDefaultProbabilityCurveTest::testFlatHazardConsistency));
         suite->add(QUANTLIB_TEST_CASE(
         &AdjointDefaultProbabilityCurveTest::testFlatDensityConsistency));
         suite->add(QUANTLIB_TEST_CASE(
         &AdjointDefaultProbabilityCurveTest::testLinearDensityConsistency));
         suite->add(QUANTLIB_TEST_CASE(
         &AdjointDefaultProbabilityCurveTest::testLogLinearSurvivalConsistency));
    return suite;
}


#if defined CL_ENABLE_BOOST_TEST_ADAPTER && defined CL_OPEN_EXCLUDED

BOOST_AUTO_TEST_SUITE(ad_default_curves)

BOOST_AUTO_TEST_CASE(testFlatHazardRate)
{
    BOOST_CHECK(AdjointDefaultProbabilityCurveTest::testFlatHazardRate());
}

BOOST_AUTO_TEST_CASE(testFlatHazardRateTime)
{
    BOOST_CHECK(AdjointDefaultProbabilityCurveTest::testFlatHazardRateTime());
}
BOOST_AUTO_TEST_CASE(testFlatHazardConsistency)
{
    BOOST_CHECK(AdjointDefaultProbabilityCurveTest::testFlatHazardConsistency());
}
BOOST_AUTO_TEST_CASE(testFlatDensityConsistency)
{
    BOOST_CHECK(AdjointDefaultProbabilityCurveTest::testFlatDensityConsistency());
}
BOOST_AUTO_TEST_CASE(testLinearDensityConsistency)
{
    BOOST_CHECK(AdjointDefaultProbabilityCurveTest::testLinearDensityConsistency());
}
BOOST_AUTO_TEST_CASE(testLogLinearSurvivalConsistency)
{
    BOOST_CHECK(AdjointDefaultProbabilityCurveTest::testLogLinearSurvivalConsistency());
}
BOOST_AUTO_TEST_SUITE_END()

#endif