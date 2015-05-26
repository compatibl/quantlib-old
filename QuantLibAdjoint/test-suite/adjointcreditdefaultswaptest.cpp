/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2008, 2009 StatPro Italia srl
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

// Based on creditdefaultswap.cpp file from test-suite.

#include "adjointcreditdefaultswaptest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include <ql/instruments/creditdefaultswap.hpp>
#include <ql/pricingengines/credit/midpointcdsengine.hpp>
#include <ql/pricingengines/credit/integralcdsengine.hpp>
#include <ql/termstructures/credit/flathazardrate.hpp>
#include <ql/termstructures/credit/interpolatedhazardratecurve.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/yield/discountcurve.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/math/interpolations/backwardflatinterpolation.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/calendars/unitedstates.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/daycounters/thirty360.hpp>
#include <iomanip>
#include <iostream>

using namespace QuantLib;
using namespace boost::unit_test_framework;


// Function std::vector<Real> jacobian(cl::TapeFunction<double>& f, std::vector<Real> const& x)
// Adapter for Jacobian.
std::vector<Real> jacobian(cl::TapeFunction<double>& f, std::vector<Real> const& x)
{
    std::vector<double> xd, yd;
    std::for_each(x.cbegin(), x.cend(),
        [&](Real const& a){ xd.push_back(double(a)); }
    );
    yd = f.Jacobian(xd);
    std::vector<Real> y;
    std::for_each(yd.cbegin(), yd.cend(),
        [&](double const& a){ y.push_back(Real(a)); }
    );
    return y;
}


// Function std::vector<Real> jacobian(cl::TapeFunction<double>& f, std::vector<Real> const& x, size_t iterNum, double& time)
// Adapter with a timer for Jacobian.
std::vector<Real> jacobian(cl::TapeFunction<double>& f, std::vector<Real> const& x, size_t iterNum, double& time)
{
    std::vector<double> xd, yd;
    std::for_each(x.cbegin(), x.cend(),
        [&](Real const& a){ xd.push_back(double(a)); }
    );

    boost::timer timer;
    for (size_t i = 0; i < iterNum; ++i)
        yd = f.Jacobian(xd);
    time = timer.elapsed() / iterNum;

    std::vector<Real> y;
    std::for_each(yd.cbegin(), yd.cend(),
        [&](double const& a){ y.push_back(Real(a)); }
    );
    return y;
}


struct RealPair
{
    Real x, y;

    static std::deque<std::string> get_columns()
    {
        static std::deque<std::string> columns = { " ", "" };
        return columns;
    }

    template <typename stream_type>
    friend inline stream_type& operator << (stream_type& stm, RealPair& p)
    {
        stm << p.x << ";" << p.y << std::endl;
        return stm;
    }
};


// Function bool AdjointCreditDefaultSwapTest::testDefaultProbabilities()
// Tests adjoint derivatives of NPV and fairSpread of CreditDefaultSwap by default probabilities for hazardRates.
bool AdjointCreditDefaultSwapTest::testDefaultProbabilities()
{
   BOOST_TEST_MESSAGE("Testing adjoint derivatives of NPV and fairSpread by default probabilities"
        " for hazardRates...");

   const Real relativeError = 1.0e-4,
       absoluteError = 1.0e-5,
       delta = 1.0e-6;

   // Minimal and maximal sizes of independent vectors and size step
   // pointsNum Points number in graphs  
   const size_t minSize = 12,
       step = 12,
#if defined CL_GRAPH_GEN
       maxSize = 600,
       pointsNum = 100;
#else
       maxSize = minSize,
       pointsNum = 1;
#endif
   
   bool result = true;

   std::vector<PerformanceTime> performTimes;
   std::vector<AdjointTime> adjointTimes;
   std::vector<TapeSize> tapeSizes;
   std::vector<RealPair> graphNPV, graphFairSpread;

   size_t vSize = (maxSize - minSize) / step + 1;
   performTimes.reserve(vSize);
   adjointTimes.reserve(vSize);
   tapeSizes.reserve(vSize);
   graphNPV.reserve(pointsNum);
   graphFairSpread.reserve(pointsNum);

   SavedSettings backup;
   Settings::instance().evaluationDate() = Date(9, June, 2006);
   Date evalDate = Settings::instance().evaluationDate();
   Calendar calendar = UnitedStates();
   std::vector<Date> discountDates = { 
       evalDate,
       calendar.advance(evalDate, 1, Weeks, ModifiedFollowing),
       calendar.advance(evalDate, 1, Months, ModifiedFollowing),
       calendar.advance(evalDate, 2, Months, ModifiedFollowing),
       calendar.advance(evalDate, 3, Months, ModifiedFollowing),
       calendar.advance(evalDate, 6, Months, ModifiedFollowing),
       calendar.advance(evalDate, 1, Years, ModifiedFollowing),
       calendar.advance(evalDate, 2, Years, ModifiedFollowing),
       calendar.advance(evalDate, 3, Years, ModifiedFollowing),
       calendar.advance(evalDate, 4, Years, ModifiedFollowing),
       calendar.advance(evalDate, 5, Years, ModifiedFollowing),
       calendar.advance(evalDate, 6, Years, ModifiedFollowing),
       calendar.advance(evalDate, 7, Years, ModifiedFollowing),
       calendar.advance(evalDate, 8, Years, ModifiedFollowing),
       calendar.advance(evalDate, 9, Years, ModifiedFollowing),
       calendar.advance(evalDate, 10, Years, ModifiedFollowing),
       calendar.advance(evalDate, 15, Years, ModifiedFollowing)
   };
   std::vector<DiscountFactor> dfs = {
       1.0,
       0.99901513757687310,
       0.99570502636871183,
       0.99118260474528685,
       0.98661167950906203,
       0.97325929533593880,
       0.94724424481038083,
       0.89844996737120875,
       0.85216647839921411,
       0.80775477692556874,
       0.76517289234200347,
       0.72401019553182933,
       0.68503909569219212,
       0.64797499814013748,
       0.61263171936255534,
       0.57919423507487910,
       0.43518868769953606
   };

   const DayCounter& curveDayCounter = Actual360();

   RelinkableHandle<YieldTermStructure> discountCurve;
   discountCurve.linkTo(
       boost::shared_ptr<YieldTermStructure>(new DiscountCurve(discountDates, dfs, curveDayCounter)));

   std::vector<Date> datesAll(maxSize);
   for (size_t i = 0; i < maxSize; ++i)
       datesAll[i] = evalDate + i * Weeks;

   std::vector<Probability> defaultProbabilitiesAll(maxSize);
   for (size_t i = 0; i < maxSize; ++i)
       defaultProbabilitiesAll[i] = 0.4 * i / maxSize;

   // Computations for different sizes of vectors of dates and defaultProbabilities
   for (size_t size = minSize; size <= maxSize; size += step)
   {
       std::vector<Date> dates(datesAll.begin(), datesAll.begin() + size);
       // Independent vector
       std::vector<Probability> defaultProbabilities(
           defaultProbabilitiesAll.begin(), defaultProbabilitiesAll.begin() + size);
       boost::timer timer;
       // Tape recording begin
       Independent(defaultProbabilities);

       std::vector<Real> hazardRates;
       DayCounter dayCounter = Thirty360();
       hazardRates.push_back(0.0);
       for (Size i = 1; i < dates.size(); ++i)
       {
           Time t1 = dayCounter.yearFraction(dates[0], dates[i - 1]);
           Time t2 = dayCounter.yearFraction(dates[0], dates[i]);
           Probability S1 = 1.0 - defaultProbabilities[i - 1];
           Probability S2 = 1.0 - defaultProbabilities[i];
           hazardRates.push_back(std::log(S1 / S2) / (t2 - t1));
       }

       RelinkableHandle<DefaultProbabilityTermStructure> piecewiseFlatHazardRate;
       piecewiseFlatHazardRate.linkTo(boost::shared_ptr<DefaultProbabilityTermStructure>(
           new InterpolatedHazardRateCurve<BackwardFlat>(dates, hazardRates, Thirty360())));

       // Testing credit default swap

       // Build the schedule
       Date issueDate(20, March, 2006);
       Date maturity = issueDate + size * Weeks;
       Frequency cdsFrequency = Weekly;
       BusinessDayConvention cdsConvention = ModifiedFollowing;

       Schedule schedule(issueDate, maturity, Period(cdsFrequency), calendar,
           cdsConvention, cdsConvention, DateGeneration::Forward, false);

       // Build the CDS
       Real recoveryRate = 0.25;
       Rate fixedRate = 0.0224;
       DayCounter dayCount = Actual360();
       Real cdsNotional = 100.0;

       CreditDefaultSwap cds(Protection::Seller, cdsNotional, fixedRate,
           schedule, cdsConvention, dayCount, true, true);
       cds.setPricingEngine(boost::shared_ptr<PricingEngine>(
           new MidPointCdsEngine(piecewiseFlatHazardRate, recoveryRate, discountCurve)));

       Real npv = cds.NPV();
       Real fairRate = cds.fairSpread();

       std::vector<Real> calculated(2);
       calculated[0] = npv;
       calculated[1] = fairRate;

       // Build Adjoint function
       // Tape recording end
       cl::TapeFunction<double> f(defaultProbabilities, calculated);
       double timeTapeRecording = timer.elapsed();
       size_t tapeSize = f.Memory();

       // Compute derivatives using Jacobian
       double timeAdjoint;
       std::vector<Real> jacob = jacobian(f, defaultProbabilities, 1000, timeAdjoint);

       // Compute finite difference derivatives using step-forward formula
       // Function derivative = (step-forward function value - current function value) / step
       size_t n = defaultProbabilities.size();
       std::vector<Real> jacobFinDiff(2 * n);
       timer.restart();
       const size_t iterNum = 10;
       for (size_t k = 0; k < iterNum; ++k)
       {
           for (size_t j = 0; j < n; ++j)
           {
               defaultProbabilities[j] += delta;
               std::vector<Real> hazardRates_;
               hazardRates_.push_back(0.0);
               for (Size i = 1; i < dates.size(); ++i)
               {
                   Time t1 = dayCounter.yearFraction(dates[0], dates[i - 1]);
                   Time t2 = dayCounter.yearFraction(dates[0], dates[i]);
                   Probability S1 = 1.0 - defaultProbabilities[i - 1];
                   Probability S2 = 1.0 - defaultProbabilities[i];
                   hazardRates_.push_back(std::log(S1 / S2) / (t2 - t1));
               }

               RelinkableHandle<DefaultProbabilityTermStructure> piecewiseFlatHazardRate_;
               piecewiseFlatHazardRate_.linkTo(
                   boost::shared_ptr<DefaultProbabilityTermStructure>(
                   new InterpolatedHazardRateCurve<BackwardFlat>(dates, hazardRates_, Thirty360())));

               cds.setPricingEngine(boost::shared_ptr<PricingEngine>(
                   new MidPointCdsEngine(piecewiseFlatHazardRate_, recoveryRate, discountCurve)));
               jacobFinDiff[j] = (cds.NPV() - npv) / delta;
               jacobFinDiff[j + n] = (cds.fairSpread() - fairRate) / delta;
               defaultProbabilities[j] -= delta;
           }
       }
       double timeAnalytical = timer.elapsed() / iterNum;

       // Adding new data to performace result vectors
       performTimes.push_back(PerformanceTime{ timeTapeRecording, timeAdjoint, timeAnalytical, size });
       adjointTimes.push_back(AdjointTime{ timeAdjoint, size });
       tapeSizes.push_back(TapeSize{ size, tapeSize });

       if (size == minSize)
       {
           Real dpDelta = (defaultProbabilities[2] - defaultProbabilities[0]) / (pointsNum + 1);
           defaultProbabilities[1] = defaultProbabilities[0];
           for (size_t j = 0; j < pointsNum; ++j)
           {
               defaultProbabilities[1] += dpDelta;
               for (Size i = 1; i <= 2; ++i)
               {
                   Time t1 = dayCounter.yearFraction(dates[0], dates[i - 1]);
                   Time t2 = dayCounter.yearFraction(dates[0], dates[i]);
                   Probability S1 = 1.0 - defaultProbabilities[i - 1];
                   Probability S2 = 1.0 - defaultProbabilities[i];
                   hazardRates[i] = std::log(S1 / S2) / (t2 - t1);
               }
               RelinkableHandle<DefaultProbabilityTermStructure> piecewiseFlatHazardRate_;
               piecewiseFlatHazardRate_.linkTo(
                   boost::shared_ptr<DefaultProbabilityTermStructure>(
                   new InterpolatedHazardRateCurve<BackwardFlat>(dates, hazardRates, Thirty360())));
               cds.setPricingEngine(boost::shared_ptr<PricingEngine>(
                   new MidPointCdsEngine(piecewiseFlatHazardRate_, recoveryRate, discountCurve)));
               graphNPV.push_back(RealPair{ defaultProbabilities[1], cds.NPV() });
               graphFairSpread.push_back(RealPair{ defaultProbabilities[1], cds.fairSpread() });
           }
       }

       // Check results
       result &= checkWithFiniteDiff(jacob, jacobFinDiff, relativeError, absoluteError);
   }

   cl::AdjointTestOutput outPerform("AdjointCreditDefaultSwap//testDefaultProbabilities", {
       { "filename", "AdjointPerformance" },
       { "not_clear", "Not" },
       { "title", "CreditDefaultSwap NPV and fair spread differentiation performance with respect to default probabilities" },
       { "xlabel", "Number of default probabilities" },
       { "ylabel", "Time (s)" } });
   outPerform << performTimes;

   cl::AdjointTestOutput outAdjoint("AdjointCreditDefaultSwap//testDefaultProbabilities", {
       { "filename", "Adjoint" },
       { "not_clear", "Not" },
       { "title", "CreditDefaultSwap NPV and fair spread adjoint differentiation performance with respect to default probabilities" },
       { "xlabel", "Number of default probabilities" },
       { "ylabel", "Time (s)" } });
   outAdjoint << adjointTimes;

   cl::AdjointTestOutput outSize("AdjointCreditDefaultSwap//testDefaultProbabilities", {
       { "filename", "TapeSize" },
       { "not_clear", "Not" },
       { "title", "Tape size dependence on number of elements in default probabilities vector" },
       { "xlabel", "Number of default probabilities" },
       { "ylabel", "Memory (MB)" } });
   outSize << tapeSizes;

   cl::AdjointTestOutput outNPV("AdjointCreditDefaultSwap//testDefaultProbabilities", {
       { "filename", "NPV on default probability" },
       { "not_clear", "Not" },
       { "title", "CreditDefaultSwap NPV dependence on default probability" },
       { "xlabel", "Default probability" },
       { "ylabel", "NPV" } });
   outNPV << graphNPV;

   cl::AdjointTestOutput outFairSpread("AdjointCreditDefaultSwap//testDefaultProbabilities", {
       { "filename", "Fair spread on default probability" },
       { "not_clear", "Not" },
       { "title", "CreditDefaultSwap fair spread dependence on default probability" },
       { "xlabel", "Default probability" },
       { "ylabel", "Fair spread" } });
   outFairSpread << graphFairSpread;

   return result;
}


// Function bool AdjointCreditDefaultSwapTest::testDiscountFactor()
// Tests adjoint derivatives of NPV and fairSpread of CreditDefaultSwap by discounts of the discount curve.
bool AdjointCreditDefaultSwapTest::testDiscountFactor()
{
    BOOST_TEST_MESSAGE("Testing adjoint derivatives of NPV and fairSpread by "
        "discounts of the discount curve...");

    const Real relativeError = 1.0e-4,
        absoluteError = 1.0e-5,
        delta = 1.0e-6;

    // Minimal and maximal sizes of independent vectors and size step
    // pointsNum Points number in graphs  
    const size_t minSize = 12,
        step = 12,
#if defined CL_GRAPH_GEN
        maxSize = 510,
        pointsNum = 100;
#else
        maxSize = minSize,
        pointsNum = 1;
#endif

    bool result = true;

    std::vector<PerformanceTime> performTimes;
    std::vector<AdjointTime> adjointTimes;
    std::vector<TapeSize> tapeSizes;
    std::vector<RealPair> graphNPV, graphFairSpread;

    size_t vSize = (maxSize - minSize) / step + 1;
    performTimes.reserve(vSize);
    adjointTimes.reserve(vSize);
    tapeSizes.reserve(vSize);
    graphNPV.reserve(pointsNum);
    graphFairSpread.reserve(pointsNum);

    SavedSettings backup;
    Settings::instance().evaluationDate() = Date(9, June, 2006);
    Date evalDate = Settings::instance().evaluationDate();
    Calendar calendar = UnitedStates();
    std::vector<Date> dates = {
        evalDate,
        calendar.advance(evalDate, 6, Months, ModifiedFollowing),
        calendar.advance(evalDate, 1, Years, ModifiedFollowing),
        calendar.advance(evalDate, 2, Years, ModifiedFollowing),
        calendar.advance(evalDate, 3, Years, ModifiedFollowing),
        calendar.advance(evalDate, 4, Years, ModifiedFollowing),
        calendar.advance(evalDate, 5, Years, ModifiedFollowing),
        calendar.advance(evalDate, 7, Years, ModifiedFollowing),
        calendar.advance(evalDate, 10, Years, ModifiedFollowing),
        calendar.advance(evalDate, 12, Years, ModifiedFollowing)
    };
    std::vector<Probability> defaultProbabilities =
    { 0.0000, 0.0047, 0.0093, 0.0286, 0.0619, 0.0953, 0.1508, 0.2288, 0.3666, 0.5 };

    std::vector<Date> discountDatesAll(maxSize + 1);
    for (size_t i = 0; i < discountDatesAll.size(); ++i)
        discountDatesAll[i] = evalDate + i * Weeks;

    std::vector<DiscountFactor> dfsAll(maxSize);
    Real coef = -1.0 / maxSize;
    for (size_t i = 0; i < maxSize; ++i)
        dfsAll[i] = std::exp(coef * (i + 1));

    // Computations for different sizes of vectors of discount dates and factors
    for (size_t size = minSize; size <= maxSize; size += step)
    {
        std::vector<Date> discountDates(discountDatesAll.begin(), discountDatesAll.begin() + size + 1);
        // Independent vector
        std::vector<Real> x(dfsAll.begin(), dfsAll.begin() + size);
        boost::timer timer;
        // Tape recording begin
        Independent(x);
        std::vector<DiscountFactor> dfs;
        dfs.push_back(1.0);
        dfs.insert(dfs.end(), x.begin(), x.end());

        const DayCounter& curveDayCounter = Actual360();

        RelinkableHandle<YieldTermStructure> discountCurve;
        discountCurve.linkTo(boost::shared_ptr<YieldTermStructure>(
            new DiscountCurve(discountDates, dfs, curveDayCounter)));

        DayCounter dayCounter = Thirty360();

        std::vector<Real> hazardRates;
        hazardRates.push_back(0.0);
        for (Size i = 1; i < dates.size(); ++i)
        {
            Time t1 = dayCounter.yearFraction(dates[0], dates[i - 1]);
            Time t2 = dayCounter.yearFraction(dates[0], dates[i]);
            Probability S1 = 1.0 - defaultProbabilities[i - 1];
            Probability S2 = 1.0 - defaultProbabilities[i];
            hazardRates.push_back(std::log(S1 / S2) / (t2 - t1));
        }

        RelinkableHandle<DefaultProbabilityTermStructure> piecewiseFlatHazardRate;
        piecewiseFlatHazardRate.linkTo(
            boost::shared_ptr<DefaultProbabilityTermStructure>(
            new InterpolatedHazardRateCurve<BackwardFlat>(dates, hazardRates, Thirty360())));

        // Testing credit default swap

        // Build the schedule
        Date issueDate(20, March, 2006);
        Date maturity = issueDate + size * Weeks;
        //Frequency cdsFrequency = Semiannual;
        Frequency cdsFrequency = Weekly;
        BusinessDayConvention cdsConvention = ModifiedFollowing;

        Schedule schedule(issueDate, maturity, Period(cdsFrequency), calendar,
            cdsConvention, cdsConvention, DateGeneration::Forward, false);

        // Build the CDS
        Real recoveryRate = 0.25;
        Rate fixedRate = 0.0224;
        DayCounter dayCount = Actual360();
        Real cdsNotional = 100.0;

        CreditDefaultSwap cds(Protection::Seller, cdsNotional, fixedRate,
            schedule, cdsConvention, dayCount, true, true);
        cds.setPricingEngine(boost::shared_ptr<PricingEngine>(
            new MidPointCdsEngine(piecewiseFlatHazardRate, recoveryRate, discountCurve)));

        Real npv = cds.NPV();
        Real fairRate = cds.fairSpread();
        std::vector<Real> y(2);
        y[0] = npv;
        y[1] = fairRate;

        // Build Adjoint function
        // Tape recording end
        cl::TapeFunction<double> f(x, y);
        double timeTapeRecording = timer.elapsed();
        size_t tapeSize = f.Memory();

        // Compute derivatives using Jacobian
        double timeAdjoint;
        std::vector<Real> jacob = jacobian(f, x, 1000, timeAdjoint);

        // Compute finite difference derivatives using step-forward formula
        // Function derivative = (step-forward function value - current function value) / step
        size_t n = x.size();
        std::vector<Real> jacobFinDiff(2 * n);
        timer.restart();
        const size_t iterNum = 10;
        for (size_t k = 0; k < iterNum; ++k)
        {
            for (size_t i = 0; i < n; ++i)
            {
                dfs[i + 1] += delta;
                RelinkableHandle<YieldTermStructure> discountCurve_;
                discountCurve_.linkTo(
                    boost::shared_ptr<YieldTermStructure>(
                    new DiscountCurve(discountDates, dfs, curveDayCounter)));
                cds.setPricingEngine(boost::shared_ptr<PricingEngine>(
                    new MidPointCdsEngine(piecewiseFlatHazardRate,
                    recoveryRate, discountCurve_)));
                jacobFinDiff[i] = (cds.NPV() - npv) / delta;
                jacobFinDiff[i + n] = (cds.fairSpread() - fairRate) / delta;
                dfs[i + 1] -= delta;
            }
        }
        double timeAnalytical = timer.elapsed() / iterNum;

        // Adding new data to performace result vectors
        performTimes.push_back(PerformanceTime{ timeTapeRecording, timeAdjoint, timeAnalytical, size });
        adjointTimes.push_back(AdjointTime{ timeAdjoint, size });
        tapeSizes.push_back(TapeSize{ size, tapeSize });

        if (size == minSize)
        {
            Real dfsDelta = (dfs[0] - dfs[2]) / (pointsNum + 1);
            dfs[1] = dfs[2];
            for (size_t i = 0; i < pointsNum; ++i)
            {
                dfs[1] += dfsDelta;
                RelinkableHandle<YieldTermStructure> discountCurve_;
                discountCurve_.linkTo(
                    boost::shared_ptr<YieldTermStructure>(
                    new DiscountCurve(discountDates, dfs, curveDayCounter)));
                cds.setPricingEngine(boost::shared_ptr<PricingEngine>(
                    new MidPointCdsEngine(piecewiseFlatHazardRate,
                    recoveryRate, discountCurve_)));
                graphNPV.push_back(RealPair{ dfs[1], cds.NPV() });
                graphFairSpread.push_back(RealPair{ dfs[1], cds.fairSpread() });
            }
        }

        // Check results
        result &= checkWithFiniteDiff(jacob, jacobFinDiff, relativeError, absoluteError);
    }

    cl::AdjointTestOutput outPerform("AdjointCreditDefaultSwap//testDiscountFactor", {
        { "filename", "AdjointPerformance" },
        { "not_clear", "Not" },
        { "title", "CreditDefaultSwap NPV and fair spred differentiation performance with respect to discount factors" },
        { "xlabel", "Number of discount factors" },
        { "ylabel", "Time (s)" } });
    outPerform << performTimes;

    cl::AdjointTestOutput outAdjoint("AdjointCreditDefaultSwap//testDiscountFactor", {
        { "filename", "Adjoint" },
        { "not_clear", "Not" },
        { "title", "CreditDefaultSwap NPV and fair spred adjoint differentiation performance with respect to discount factors" },
        { "xlabel", "Number of discount factors" },
        { "ylabel", "Time (s)" } });
    outAdjoint << adjointTimes;

    cl::AdjointTestOutput outSize("AdjointCreditDefaultSwap//testDiscountFactor", {
        { "filename", "TapeSize" },
        { "not_clear", "Not" },
        { "title", "Tape size dependence on number of elements in discount factors vector" },
        { "xlabel", "Number of discount factors" },
        { "ylabel", "Memory (MB)" } });
    outSize << tapeSizes;

    cl::AdjointTestOutput outNPV("AdjointCreditDefaultSwap//testDiscountFactor", {
        { "filename", "NPV on discount factor" },
        { "not_clear", "Not" },
        { "title", "CreditDefaultSwap NPV dependence on discount factor" },
        { "xlabel", "Discount factor" },
        { "ylabel", "NPV" } });
    outNPV << graphNPV;

    cl::AdjointTestOutput outFairSpread("AdjointCreditDefaultSwap//testDiscountFactor", {
        { "filename", "Fair spred on discount factor" },
        { "not_clear", "Not" },
        { "title", "CreditDefaultSwap fair spred dependence on discount factor" },
        { "xlabel", "Discount factor" },
        { "ylabel", "Fair spred" } });
    outFairSpread << graphFairSpread;

    return result;
}


// Function bool AdjointCreditDefaultSwapTest::testHazardRate()
// Testing adjoint derivatives of NPV and impliedHazardRate of CreditDefaultSwap by hazardRates.
bool AdjointCreditDefaultSwapTest::testHazardRate()
{
    BOOST_TEST_MESSAGE("Testing adjoint derivatives of NPV and impliedHazardRate by hazardRates...");

    const Real relativeError = 1.0e-4,
        absoluteError = 1.0e-5,
        delta = 1.0e-6;

    // Minimal and maximal sizes of independent vectors and size step
    // pointsNum Points number in graphs  
    const size_t minSize = 12,
        step = 12,
#if defined CL_GRAPH_GEN
        maxSize = 610,
        pointsNum = 100;
#else
        maxSize = minSize,
        pointsNum = 2;
#endif

    bool result = true;

    std::vector<PerformanceTime> performTimes;
    std::vector<AdjointTime> adjointTimes;
    std::vector<TapeSize> tapeSizes;
    std::vector<RealPair> graphNPV, graphImpliedHazardRate;

    size_t vSize = (maxSize - minSize) / step + 1;
    performTimes.reserve(vSize);
    adjointTimes.reserve(vSize);
    tapeSizes.reserve(vSize);
    graphNPV.reserve(pointsNum);
    graphImpliedHazardRate.reserve(pointsNum);

    SavedSettings backup;
    Calendar calendar = TARGET();
    Date today = calendar.adjust(Date::todaysDate());
    Settings::instance().evaluationDate() = today;

    std::vector<Date> datesAll(maxSize);
    for (size_t i = 0; i < maxSize; ++i)
        datesAll[i] = today + i * Weeks;

    std::vector<Real> hazardRatesAll(maxSize);
    for (size_t i = 0; i < maxSize; ++i)
        hazardRatesAll[i] = 0.3 + 0.4 * i / maxSize;

    // Computations for different sizes of vectors of dates and hazardRates
    for (size_t size = minSize; size <= maxSize; size += step)
    {
        std::vector<Date> dates(datesAll.begin(), datesAll.begin() + size);
        // Independent vector
        std::vector<Real> hazardRates(hazardRatesAll.begin(), hazardRatesAll.begin() + size);
        boost::timer timer;
        // Tape recording begin
        Independent(hazardRates);

        // Initialize curves
        DayCounter dayCounter = Actual365Fixed();

        RelinkableHandle<DefaultProbabilityTermStructure> probabilityCurve;
        probabilityCurve.linkTo(boost::shared_ptr<DefaultProbabilityTermStructure>(
            new InterpolatedHazardRateCurve<BackwardFlat>(dates, hazardRates, dayCounter)));

        RelinkableHandle<YieldTermStructure> discountCurve;
        discountCurve.linkTo(boost::shared_ptr<YieldTermStructure>(
            new FlatForward(today, 0.03, Actual360())));

        Frequency frequency = Weekly;
        BusinessDayConvention convention = ModifiedFollowing;

        // Build the CDS
        Date issueDate = calendar.advance(today, -6, Months);
        Rate fixedRate = 0.0120;
        DayCounter cdsDayCount = Actual360();
        Real notional = 10000.0;
        Real recoveryRate = 0.4;

        Date maturity = calendar.advance(issueDate, size / 2, Weeks);
        Schedule schedule(issueDate, maturity, Period(frequency), calendar,
            convention, convention, DateGeneration::Forward, false);

        CreditDefaultSwap cds(Protection::Seller, notional, fixedRate,
            schedule, convention, cdsDayCount, true, true);
        cds.setPricingEngine(boost::shared_ptr<PricingEngine>(
            new MidPointCdsEngine(probabilityCurve, recoveryRate, discountCurve)));

        std::vector<Real> computed(2);
        computed[0] = cds.NPV();
        computed[1] = cds.impliedHazardRate(computed[0], discountCurve, dayCounter, recoveryRate);

        // Build Adjoint function
        // Tape recording end
        cl::TapeFunction<double> f(hazardRates, computed);
        double timeTapeRecording = timer.elapsed();
        size_t tapeSize = f.Memory();

        // Compute derivatives using Jacobian
        double timeAdjoint;
        std::vector<Real> jacob = jacobian(f, hazardRates, 1000, timeAdjoint);

        // Compute finite difference derivatives using step-forward formula
        // Function derivative = (step-forward function value - current function value) / step
        std::vector<Real> jacobFinDiff(2 * size);
        timer.restart();
        const size_t iterNum = 10;
        for (size_t k = 0; k < iterNum; ++k)
        {
            for (size_t i = 0; i < size; ++i)
            {
                hazardRates[i] += delta;
                RelinkableHandle<DefaultProbabilityTermStructure> probabilityCurve_;
                probabilityCurve_.linkTo(boost::shared_ptr<DefaultProbabilityTermStructure>(
                    new InterpolatedHazardRateCurve<BackwardFlat>(dates, hazardRates, dayCounter)));
                cds.setPricingEngine(boost::shared_ptr<PricingEngine>(
                    new MidPointCdsEngine(probabilityCurve_, recoveryRate, discountCurve)));
                Real npv_ = cds.NPV();
                jacobFinDiff[i] = (npv_ - computed[0]) / delta;
                jacobFinDiff[i + size] =
                    (cds.impliedHazardRate(npv_, discountCurve, dayCounter, recoveryRate) - computed[1]) / delta;
                hazardRates[i] -= delta;
            }
        }
        double timeAnalytical = timer.elapsed() / iterNum;

        // Adding new data to performace result vectors
        performTimes.push_back(PerformanceTime{ timeTapeRecording, timeAdjoint, timeAnalytical, size });
        adjointTimes.push_back(AdjointTime{ timeAdjoint, size });
        tapeSizes.push_back(TapeSize{ size, tapeSize });

        if (size == minSize)
        {
            Real hrDelta = (hazardRates[1] - hazardRates[0]) / (pointsNum - 1);
            for (size_t i = 0; i < pointsNum; ++i)
            {
                RelinkableHandle<DefaultProbabilityTermStructure> probabilityCurve_;
                probabilityCurve_.linkTo(boost::shared_ptr<DefaultProbabilityTermStructure>(
                    new InterpolatedHazardRateCurve<BackwardFlat>(dates, hazardRates, dayCounter)));
                cds.setPricingEngine(boost::shared_ptr<PricingEngine>(
                    new MidPointCdsEngine(probabilityCurve_, recoveryRate, discountCurve)));
                Real npv = cds.NPV();
                graphNPV.push_back(RealPair{ hazardRates[0], npv });
                graphImpliedHazardRate.push_back(RealPair{ hazardRates[0], cds.impliedHazardRate(npv, discountCurve, dayCounter, recoveryRate) });
                hazardRates[0] += hrDelta;
            }
        }

        // Check results
        result &= checkWithFiniteDiff(jacob, jacobFinDiff, relativeError, absoluteError);
    }

    cl::AdjointTestOutput outPerform("AdjointCreditDefaultSwap//testHazardRate", {
        { "filename", "AdjointPerformance" },
        { "not_clear", "Not" },
        { "title", "CreditDefaultSwap NPV and implied hazard rate differentiation performance with respect to hazard rates" },
        { "xlabel", "Number of hazard rates" },
        { "ylabel", "Time (s)" } });
    outPerform << performTimes;

    cl::AdjointTestOutput outAdjoint("AdjointCreditDefaultSwap//testHazardRate", {
        { "filename", "Adjoint" },
        { "not_clear", "Not" },
        { "title", "CreditDefaultSwap NPV and implied hazard rate adjoint differentiation performance with respect to hazard rates" },
        { "xlabel", "Number of hazard rates" },
        { "ylabel", "Time (s)" } });
    outAdjoint << adjointTimes;

    cl::AdjointTestOutput outSize("AdjointCreditDefaultSwap//testHazardRate", {
        { "filename", "TapeSize" },
        { "not_clear", "Not" },
        { "title", "Tape size dependence on number of elements in hazard rates vector" },
        { "xlabel", "Number of hazard rates" },
        { "ylabel", "Memory (MB)" } });
    outSize << tapeSizes;

    cl::AdjointTestOutput outNPV("AdjointCreditDefaultSwap//testHazardRate", {
        { "filename", "NPV on hazard rate" },
        { "not_clear", "Not" },
        { "title", "CreditDefaultSwap NPV dependence on hazard rate" },
        { "xlabel", "Hazard rate" },
        { "ylabel", "NPV" } });
    outNPV << graphNPV;

    cl::AdjointTestOutput outImpliedHazardRate("AdjointCreditDefaultSwap//testHazardRate", {
        { "filename", "Implied hazard rate on hazard rate" },
        { "not_clear", "Not" },
        { "title", "CreditDefaultSwap implied hazard rate dependence on hazard rate" },
        { "xlabel", "Hazard rate" },
        { "ylabel", "Implied hazard rate" } });
    outImpliedHazardRate << graphImpliedHazardRate;

    return result;
}


// Function bool AdjointCreditDefaultSwapTest::testNotionalSpredRate()
// Tests adjoint derivatives of NPV and fairSpread of CreditDefaultSwap by notional, spred, and recoveryRate of the pricing engine.
bool AdjointCreditDefaultSwapTest::testNotionalSpredRate()
{
    BOOST_TEST_MESSAGE("Testing adjoint derivatives of NPV and fairSpread by notional, spred, and "
        "recoveryRate of the pricing engine...");

    const Real relativeError = 1.0e-4,
        absoluteError = 1.0e-5,
        delta = 1.0e-4;

    // Initialize curves
    SavedSettings backup;
    Settings::instance().evaluationDate() = Date(9, June, 2006);
    Date today = Settings::instance().evaluationDate();
    Calendar calendar = TARGET();

    Handle<Quote> hazardRate = Handle<Quote>(
        boost::shared_ptr<Quote>(new SimpleQuote(0.01234)));
    RelinkableHandle<DefaultProbabilityTermStructure> probabilityCurve;
    probabilityCurve.linkTo(
        boost::shared_ptr<DefaultProbabilityTermStructure>(
        new FlatHazardRate(0, calendar, hazardRate, Actual360())));

    RelinkableHandle<YieldTermStructure> discountCurve;

    discountCurve.linkTo(boost::shared_ptr<YieldTermStructure>(
        new FlatForward(today, 0.06, Actual360())));

    // Build the schedule
    Date issueDate = calendar.advance(today, -1, Years);
    Date maturity = calendar.advance(issueDate, 10, Years);
    Frequency frequency = Semiannual;
    BusinessDayConvention convention = ModifiedFollowing;

    Schedule schedule(issueDate, maturity, Period(frequency), calendar,
        convention, convention, DateGeneration::Forward, false);

    // Independent vector
    std::vector<Real> x(3), y(4);
    x[0] = 0.0120;
    x[1] = 10000.0;
    x[2] = 0.4;
    // Tape recording begin
    Independent(x);

    // Build the CDS
    Real notional = x[0];
    Rate fixedRate = x[1];
    Real recoveryRate = x[2];
    DayCounter dayCount = Actual360();
    CreditDefaultSwap cds(Protection::Seller, notional, fixedRate, schedule, convention, dayCount, true, true);

    auto midPointCdsEngine = boost::shared_ptr<PricingEngine>(
        new MidPointCdsEngine(probabilityCurve, recoveryRate, discountCurve));
    cds.setPricingEngine(midPointCdsEngine);
    y[0] = cds.NPV();
    y[1] = cds.fairSpread();

    auto integralCdsEngine = boost::shared_ptr<PricingEngine>(
        new IntegralCdsEngine(1 * Days, probabilityCurve, recoveryRate, discountCurve));
    cds.setPricingEngine(integralCdsEngine);
    y[2] = cds.NPV();
    y[3] = cds.fairSpread();

    // Build Adjoint function
    // Tape recording end
    cl::TapeFunction<double> f(x, y);
    // Compute derivatives using Jacobian
    std::vector<Real> jacob = jacobian(f, x);

    // Compute finite difference derivatives using step-forward formula
    // Function derivative = (step-forward function value - current function value) / step
    std::vector<Real> jacobFinDiff;
    std::vector<CreditDefaultSwap> cdsv;
    cdsv.emplace_back(Protection::Seller, notional + delta, fixedRate, schedule, convention, dayCount, true, true);
    cdsv.emplace_back(Protection::Seller, notional, fixedRate + delta, schedule, convention, dayCount, true, true);
    cdsv.emplace_back(Protection::Seller, notional, fixedRate, schedule, convention, dayCount, true, true);

    cdsv[0].setPricingEngine(midPointCdsEngine);
    cdsv[1].setPricingEngine(midPointCdsEngine);
    cdsv[2].setPricingEngine(boost::shared_ptr<PricingEngine>(
        new MidPointCdsEngine(probabilityCurve, recoveryRate + delta, discountCurve)));
    for (int i = 0; i < 3; ++i)
        jacobFinDiff.push_back((cdsv[i].NPV() - y[0]) / delta);
    for (int i = 0; i < 3; ++i)
        jacobFinDiff.push_back((cdsv[i].fairSpread() - y[1]) / delta);

    cdsv[0].setPricingEngine(integralCdsEngine);
    cdsv[1].setPricingEngine(integralCdsEngine);
    cdsv[2].setPricingEngine(boost::shared_ptr<PricingEngine>(
        new IntegralCdsEngine(1 * Days, probabilityCurve, recoveryRate + delta, discountCurve)));
    for (int i = 0; i < 3; ++i)
        jacobFinDiff.push_back((cdsv[i].NPV() - y[2]) / delta);
    for (int i = 0; i < 3; ++i)
        jacobFinDiff.push_back((cdsv[i].fairSpread() - y[3]) / delta);

    // Check results
    return checkWithFiniteDiff(jacob, jacobFinDiff, relativeError, absoluteError);
}


test_suite* AdjointCreditDefaultSwapTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Adjoint credit-default swap tests");
    suite->add(QUANTLIB_TEST_CASE(&testDefaultProbabilities));
    suite->add(QUANTLIB_TEST_CASE(&testDiscountFactor));
    suite->add(QUANTLIB_TEST_CASE(&testHazardRate));
    suite->add(QUANTLIB_TEST_CASE(&testNotionalSpredRate));
    return suite;
}


#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_credit_default_swap)

BOOST_AUTO_TEST_CASE(testCreditDefaultSwapDefaultProbabilities)
{
    BOOST_CHECK(AdjointCreditDefaultSwapTest::testDefaultProbabilities());
}
BOOST_AUTO_TEST_CASE(testCreditDefaultSwapDiscountFactor)
{
    BOOST_CHECK(AdjointCreditDefaultSwapTest::testDiscountFactor());
}
BOOST_AUTO_TEST_CASE(testCreditDefaultSwapHazardRate)
{
    BOOST_CHECK(AdjointCreditDefaultSwapTest::testHazardRate());
}
BOOST_AUTO_TEST_CASE(testCreditDefaultSwapNotionalSpredRate)
{
    BOOST_CHECK(AdjointCreditDefaultSwapTest::testNotionalSpredRate());
}

BOOST_AUTO_TEST_SUITE_END()

#endif


