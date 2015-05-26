/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2014 StatPro Italia srl
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

// Based on piecewisezerospreadedtermstructure.cpp from Quantlib/test-suite.

#include "adjointpiecewisezerospreadedtermstructuretest.hpp"
#include "utilities.hpp"
#include <ql/termstructures/yield/piecewisezerospreadedtermstructure.hpp>
#include <ql/termstructures/yield/zerocurve.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/termstructures/yield/ratehelpers.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/thirty360.hpp>
#include <ql/math/interpolations/all.hpp>
#include <boost/make_shared.hpp>
#include "adjointtestutilities.hpp"

using namespace QuantLib;
using namespace boost::unit_test_framework;

void outputFlatInterpolationLeft(cl::AdjointTestOutput & logOut);
void outputFlatInterpolationRight(cl::AdjointTestOutput & logOut);
void outputLinearInterpolationMultipleSpreads(cl::AdjointTestOutput & logOut);
void outputLinearInterpolation(cl::AdjointTestOutput & logOut);
void outputForwardFlatInterpolation(cl::AdjointTestOutput & logOut);
void outputBackwardFlatInterpolation(cl::AdjointTestOutput & logOut);
void outputDefaultInterpolation(cl::AdjointTestOutput & logOut);
void outputSetInterpolationFactory(cl::AdjointTestOutput & logOut);
void outputQuoteChanging(cl::AdjointTestOutput & logOut);

namespace {

#if defined CL_GRAPH_GEN
    Size MIN_RATE_SIZE = 8;
    Size MAX_RATE_SIZE = 110;
    Real MIN_RATE_RANGE = 0.0;
    int repeatNumber = 10000; //Jacobian matrix calculated to fast and timer can't notice the time
#else
    Size MIN_RATE_SIZE = 30;
    Size MAX_RATE_SIZE = 30;
    Real MIN_RATE_RANGE = 1.0;
    int repeatNumber = 1;
#endif

    

    Real MAX_RATE_RANGE = 1.0;
    Real RATE_STEP = 0.01;

    std::vector<double> getDoubleVector(Size size)
    {
        std::vector<double> vec =
        {
              0.04, 0.09, 0.12, 0.06, 0.05, 0.19, 0.52, 0.66, 0.45, 0.82, 0.73, 0.33, 0.82
            , 0.66, 0.39, 0.74, 0.58, 0.55, 0.36, 0.91, 0.11, 0.24, 0.5,  0.22, 0.43, 0.5
            , 0.99, 0.64, 0.58, 0.35, 0.46, 0.86, 0.71, 0.23, 0.3,  0.07, 0.03, 0.82, 0.12
            , 0.55, 0.82, 0.95, 0.11, 0.91, 0.16, 0.73, 0.24, 0.96, 0.31, 0.32, 0.46, 0.88
            , 0.71, 0.03, 0.55, 0.58, 0.1,  0.63, 0.23, 0.02, 0.62, 0.7,  0.37, 0.21, 0.19
            , 0.13, 0.2,  0.46, 0.41, 0.9,  0.99, 0.01, 0.73, 0.19, 0.2,  0.58, 0.59, 0.83
            , 0.1,  0.85, 0.76, 0.9,  0.49, 0.01, 0.65, 0.99, 0.5,  0.63, 0.52, 0.43, 0.79
            , 0.49, 0.69, 0.01, 0.44, 0.22, 0.64, 0.88, 0.83, 0.35, 0.54, 0.16, 0.91, 0.18
            , 0.35, 0.45, 0.26, 0.48, 0.67, 0.94, 0.03, 0.22, 0.64, 0.88, 0.83, 0.35, 0.64
        };
        return std::vector<double>(vec.begin(), vec.begin() + size);
    }

    std::vector<Integer> getIntegerVector(Size size)
    {
        std::vector<Integer> vec = {
            18, 51, 98, 120, 170, 201, 296, 297, 315, 383, 394, 460, 504, 588, 620, 691, 753
            , 776, 781, 806, 904, 993, 1062, 1107, 1126, 1196, 1199, 1236, 1301, 1347, 1382
            , 1475, 1486, 1539, 1576, 1664, 1740, 1793, 1799, 1837, 1906, 1942, 1982, 2070
            , 2089, 2115, 2138, 2196, 2210, 2244, 2293, 2368, 2425, 2445, 2494, 2572, 2588, 2674
            , 2770, 2790, 2792, 2802, 2810, 2903, 2973, 2990, 3032, 3108, 3193, 3196, 3289
            , 3346, 3349, 3356, 3415, 3468, 3492, 3538, 3578, 3635, 3664, 3760, 3781, 3875
            , 3930, 4025, 4047, 4072, 4150, 4216, 4295, 4388, 4441, 4521, 4572, 4594, 4615, 4692
            , 4766, 4778, 4834, 4856, 4895, 4936, 4965, 4995, 5012, 5096, 5145, 5185, 5215, 5305
        };
        return std::vector<Integer>(vec.begin(), vec.begin() + size);
    }

    struct RateDependence
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "Rate", ""
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, RateDependence& v)
        {
                stm << v.inputRate_
                    << ";" << v.interpolatedZeroRate_ << std::endl;

                return stm;
            }

        Real inputRate_;
        Real interpolatedZeroRate_;
    };

    struct Datum {
        Integer n;
        TimeUnit units;
        Rate rate;
    };

    struct CommonVars {
        // common data
        Calendar calendar;
        Natural settlementDays;
        DayCounter dayCount;
        Compounding compounding;
        boost::shared_ptr<YieldTermStructure> termStructure;
        Date today;
        Date settlementDate;

        // cleanup
        SavedSettings backup;

        // setup
        CommonVars(std::vector<Rate> rates, std::vector<Integer> ts) {
            calendar = TARGET();
            settlementDays = 2;
            today = Date(9, June, 2009);
            compounding = Continuous;
            dayCount = Actual360();
            settlementDate = calendar.advance(today, settlementDays, Days);

            Settings::instance().evaluationDate() = today;

            std::vector<Date> dates(1, settlementDate);
            for (Size i = 0; i < rates.size() - 1; ++i) {
                dates.push_back(calendar.advance(today, ts[i], Days));
            }
            termStructure = boost::make_shared<ZeroCurve>(dates, rates, dayCount);
        }
    };

}

// Method testFlatInterpolationLeft()
// Testing flat interpolation before the first spreaded date.
// zeroRate is interpolated in a date point by the above mentioned method using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays) + spread quote.
// The output quantities (interpolatedZeroRate and expectedRate) are differentiated
// with respect to input yield rate vector values.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same plot.
// Flat Interpolation Left returns the output value corresponding to the node value that is immediately less than the input value.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it throws an exception and returns false.
bool AdjointPiecewiseZeroSpreadedTermStructureTest::testFlatInterpolationLeft() {

    BOOST_MESSAGE("Testing flat interpolation before the first spreaded date...");
    std::vector<PerformanceTime> pTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<AdjointTime> adjTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<TapeSize> tapeSize(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    cl::AdjointTestOutput outPerform("AdjointPiecwZeroSprTermStruct//FlatInterpLeft", { { "filename", "AdjointPerformance" }
                                                                                   , { "not_clear", "Not" }
                                                                                   , { "title", "Interpolated zero rate (by Flat Interp. Left) differentiation performance with respect to yield rate" }
                                                                                   , { "ylabel", "Time (s)" }
                                                                                   , { "xlabel", "Number of yield rates" }
                                                                                   , { "line_box_width", "-5" } });
    cl::AdjointTestOutput outSize("AdjointPiecwZeroSprTermStruct//FlatInterpLeft", { { "filename", "TapeSize" }
                                                                                   , { "not_clear", "Not" }
                                                                                   , { "title", "Tape size dependence on number of yield rates" }
                                                                                   , { "xlabel", "Number of yield rates" }
                                                                                   , { "ylabel", "Memory(MB)" } });
    cl::AdjointTestOutput outAdjoint("AdjointPiecwZeroSprTermStruct//FlatInterpLeft", { { "filename", "Adjoint" }
                                                                                  , { "not_clear", "Not" }
                                                                                  , { "title", "Interpolated zero rate (by Flat Interp. Left) adjoint differentiation with respect to yield rate" }
                                                                                  , { "xlabel", "Number of yield rates" }
                                                                                  , { "ylabel", "Adjoint calculations time (s)" } });
    bool ok = true;

    for (Size rate_size = MIN_RATE_SIZE; rate_size <= MAX_RATE_SIZE; rate_size++)
    {
        pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        tapeSize[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        std::vector<double> rates_double = getDoubleVector(rate_size);
        std::vector<Integer> ts = getIntegerVector(rates_double.size() - 1);
        std::vector<Rate> rates(rates_double.size());

        for (Size i = 0; i < rates.size(); i++)
        {
            rates[i] = Real(rates_double[i]);
        }
        outPerform.log() << std::endl << "Start of tape recording " << currentTime() << std::endl;
        boost::timer timer;
        Independent(rates);
        CommonVars vars(rates, ts);

        std::vector<Handle<Quote> > spreads;
        boost::shared_ptr<SimpleQuote> spread1 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread2 = boost::make_shared<SimpleQuote>(0.03);
        spreads.push_back(Handle<Quote>(spread1));
        spreads.push_back(Handle<Quote>(spread2));

        std::vector<Date> spreadDates;
        spreadDates.push_back(vars.calendar.advance(vars.today, 8, Months));
        spreadDates.push_back(vars.calendar.advance(vars.today, 15, Months));

        Date interpolationDate = vars.calendar.advance(vars.today, 6, Months);

        boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
            boost::make_shared<PiecewiseZeroSpreadedTermStructure>(
            Handle<YieldTermStructure>(vars.termStructure),
            spreads, spreadDates);

        Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
        Rate interpolatedZeroRate = spreadedTermStructure->zeroRate(t, vars.compounding);
        Real expectedRate = vars.termStructure->zeroRate(t, vars.compounding) +
            spread1->value();
        pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ = timer.elapsed();

        outPerform.log() << "End of tape recording " << currentTime() << std::endl;

        std::vector<Rate> resultRates = { interpolatedZeroRate, expectedRate };
        cl::TapeFunction<double> f(rates, resultRates);
        tapeSize[rate_size - MIN_RATE_SIZE].memory_ = f.Memory();
        outPerform.log() << "Tape Memory in bytes: " << f.Memory() << std::endl;
        outPerform.log() << "Time for tape recording: " << pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ << std::endl;
        //Compute derivatives in forward and reverse mode
        std::vector<double> dy_forward(resultRates.size() * rates.size());
        std::vector<double> dy_reverce(resultRates.size() * rates.size());
        pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = std::min(gradReverse(f, dy_reverce, outPerform, false, false, repeatNumber)
                                                               , gradForward(f, dy_forward, outPerform, false, false, repeatNumber));
        double h = 1.0e-6;
        std::vector<Real> dy_diff(resultRates.size() * rates.size());
        timer.restart();
        //Finite differences
        for (Size i = 0; i < rates.size(); i++)
        {
            rates[i] += h;
            CommonVars vars(rates, ts);

            std::vector<Handle<Quote> > spreads;
            boost::shared_ptr<SimpleQuote> spread1 = boost::make_shared<SimpleQuote>(0.02);
            boost::shared_ptr<SimpleQuote> spread2 = boost::make_shared<SimpleQuote>(0.03);
            spreads.push_back(Handle<Quote>(spread1));
            spreads.push_back(Handle<Quote>(spread2));

            std::vector<Date> spreadDates;
            spreadDates.push_back(vars.calendar.advance(vars.today, 8, Months));
            spreadDates.push_back(vars.calendar.advance(vars.today, 15, Months));

            Date interpolationDate = vars.calendar.advance(vars.today, 6, Months);

            boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
                boost::make_shared<PiecewiseZeroSpreadedTermStructure>(
                Handle<YieldTermStructure>(vars.termStructure),
                spreads, spreadDates);

            Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
            Rate interpolatedZeroRate_h = spreadedTermStructure->zeroRate(t, vars.compounding);


            Real expectedRate_h = vars.termStructure->zeroRate(t, vars.compounding) +
                spread1->value();

            dy_diff[i] = (interpolatedZeroRate_h - interpolatedZeroRate) / h;
            dy_diff[rates.size() + i] = (expectedRate_h - expectedRate) / h;
            rates[i] -= h;
        }
        pTime[rate_size - MIN_RATE_SIZE].timeAnalytical_ = timer.elapsed();
        adjTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_;
        adjTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_;
        ok &= checkWithFiniteDiff(dy_forward, dy_reverce, dy_diff, outPerform, 1e-7, 1e-7);

    }
    outPerform << pTime;
    outPerform.log() << "Flat Interpolation Left Performance Time plot has been generated " << currentTime() << std::endl;
    outAdjoint << adjTime;
    outPerform.log() << "Flat Interpolation Left Adjoint Performance plot has been generated " << currentTime() << std::endl;
    outSize << tapeSize;
    outPerform.log() << "Tape size dependence on number of yield rates ploy has been generated" << std::endl;
    outputFlatInterpolationLeft(outPerform);
    return ok;
}

// Method testFlatInterpolationRight()
// Testing flat interpolation after the last spreaded date.
// zeroRate is interpolated in a date point by the above mentioned method using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays) + spread quote.
// The output quantities (interpolatedZeroRate and expectedRate) are differentiated
// with respect to input yield rate vector values.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same plot.
// Flat Interpolation Right returns the output value corresponding to the node value that is immediately greater than the input value.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it throws an exception and returns false.
bool AdjointPiecewiseZeroSpreadedTermStructureTest::testFlatInterpolationRight() {

    BOOST_MESSAGE("Testing flat interpolation after the last spreaded date...");
    std::vector<PerformanceTime> pTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<AdjointTime> adjTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<TapeSize> tapeSize(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    cl::AdjointTestOutput outPerform("AdjointPiecwZeroSprTermStruct//FlatInterpRight", { { "filename", "AdjointPerformance" }
                                                                                    , { "not_clear", "Not" }
                                                                                    , { "title", "Interpolated zero rate(by Flat Interp. Right) differentiation performance with respect to yield rate" }
                                                                                    , { "ylabel", "Time (s)" }
                                                                                    , { "xlabel", "Number of yield rates" }
                                                                                    , { "line_box_width", "-5" } });
    cl::AdjointTestOutput outAdjoint("AdjointPiecwZeroSprTermStruct//FlatInterpRight", { { "filename", "Adjoint" }
                                                                                   , { "not_clear", "Not" }
                                                                                   , { "title", "Interpolated zero rate (by Flat Interp. Right) adjoint differentiation with respect to yield rate" }
                                                                                   , { "xlabel", "Number of yield rates" }
                                                                                   , { "ylabel", "Adjoint calculations time (s)" } });
    cl::AdjointTestOutput outSize("AdjointPiecwZeroSprTermStruct//FlatInterpRight", { { "filename", "TapeSize" }
                                                                                    , { "not_clear", "Not" }
                                                                                    , { "title", "Tape size dependence on number of yield rates" }
                                                                                    , { "xlabel", "Number of yield rates" }
                                                                                    , { "ylabel", "Memory(MB)" } });

    bool ok = true;

    for (Size rate_size = MIN_RATE_SIZE; rate_size <= MAX_RATE_SIZE; rate_size++)
    {
        pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        tapeSize[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        std::vector<double> rates_double = getDoubleVector(rate_size);
        std::vector<Integer> ts = getIntegerVector(rates_double.size() - 1);
        std::vector<Rate> rates(rates_double.size());
        for (Size i = 0; i < rates.size(); i++)
        {
            rates[i] = Real(rates_double[i]);
        }
        outPerform.log() << std::endl << "Start of tape recording " << currentTime() << std::endl;
        boost::timer timer;
        Independent(rates);
        CommonVars vars(rates, ts);

        std::vector<Handle<Quote> > spreads;
        boost::shared_ptr<SimpleQuote> spread1 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread2 = boost::make_shared<SimpleQuote>(0.03);
        spreads.push_back(Handle<Quote>(spread1));
        spreads.push_back(Handle<Quote>(spread2));

        std::vector<Date> spreadDates;
        spreadDates.push_back(vars.calendar.advance(vars.today, 8, Months));
        spreadDates.push_back(vars.calendar.advance(vars.today, 15, Months));

        Date interpolationDate = vars.calendar.advance(vars.today, 6, Months);

        boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
            boost::make_shared<PiecewiseZeroSpreadedTermStructure>(
            Handle<YieldTermStructure>(vars.termStructure),
            spreads, spreadDates);
        spreadedTermStructure->enableExtrapolation();

        Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
        Rate interpolatedZeroRate = spreadedTermStructure->zeroRate(t, vars.compounding);


        Real expectedRate = vars.termStructure->zeroRate(t, vars.compounding) +
            spread2->value();
        pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ = timer.elapsed();
        outPerform.log() << "End of tape recording " << currentTime() << std::endl;
        std::vector<Rate> resultRates = { interpolatedZeroRate, expectedRate };
        cl::TapeFunction<double> f(rates, resultRates);
        tapeSize[rate_size - MIN_RATE_SIZE].memory_ = f.Memory();
        outPerform.log() << "Tape Memory in bytes: " << f.Memory() << std::endl;
        outPerform.log() << "Time for tape recording: " << pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ << std::endl;
        //Compute derivatives in forward and reverse mode
        std::vector<double> jacobian;
        timer.restart();
        std::vector<double> dy_forward(resultRates.size() * rates.size());
        std::vector<double> dy_reverce(resultRates.size() * rates.size());
        pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = std::min(gradReverse(f, dy_reverce, outPerform, false, false, repeatNumber)
                                                               , gradForward(f, dy_forward, outPerform, false, false, repeatNumber));
        double h = 1.0e-6;
        std::vector<Real> dy_diff(resultRates.size() * rates.size());
        timer.restart();
        //Finite differences
        for (Size i = 0; i < rates.size(); i++)
        {
            rates[i] += h;
            CommonVars vars(rates, ts);

            std::vector<Date> spreadDates;
            spreadDates.push_back(vars.calendar.advance(vars.today, 8, Months));
            spreadDates.push_back(vars.calendar.advance(vars.today, 15, Months));

            Date interpolationDate = vars.calendar.advance(vars.today, 6, Months);

            boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
                boost::make_shared<PiecewiseZeroSpreadedTermStructure>(
                Handle<YieldTermStructure>(vars.termStructure),
                spreads, spreadDates);
            spreadedTermStructure->enableExtrapolation();

            Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
            Rate interpolatedZeroRate_h = spreadedTermStructure->zeroRate(t, vars.compounding);


            Real expectedRate_h = vars.termStructure->zeroRate(t, vars.compounding) +
                spread2->value();

            dy_diff[i] = (interpolatedZeroRate_h - interpolatedZeroRate) / h;
            dy_diff[rates.size() + i] = (expectedRate_h - expectedRate) / h;
            rates[i] -= h;
        }
        pTime[rate_size - MIN_RATE_SIZE].timeAnalytical_ = timer.elapsed();
        adjTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_;
        adjTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_;
        ok &= checkWithFiniteDiff(dy_forward, dy_reverce, dy_diff, outPerform, 1e-7, 1e-7);
    }
    outPerform << pTime;
    outPerform.log() << "Flat Interpolation Right Performance Time plot has been generated " << currentTime() << std::endl;
    outAdjoint << adjTime;
    outPerform.log() << "Flat Interpolation Right Adjoint Performance plot has been generated " << currentTime() << std::endl;
    outSize << tapeSize;
    outPerform.log() << "Tape size dependence on number of yield rates plot has been generated " << currentTime() << std::endl;
    outputFlatInterpolationRight(outPerform);
    return ok;
}

// Method testLinearInterpolationMultipleSpreads()
// Testing linear interpolation between two dates.
// zeroRate is interpolated in a date point by the above mentioned method using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays) + spread quote.
// The output quantities (interpolatedZeroRate and expectedRate) are differentiated
// with respect to input yield rate vector values.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same plot.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it throws an exception and returns false.
bool AdjointPiecewiseZeroSpreadedTermStructureTest::testLinearInterpolationMultipleSpreads() {

    BOOST_MESSAGE("Testing linear interpolation with more than two spreaded dates...");
    std::vector<PerformanceTime> pTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<AdjointTime> adjTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<TapeSize> tapeSize(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    cl::AdjointTestOutput outPerform("AdjointPiecwZeroSprTermStruct//LinInterpMultSpreads", { { "filename", "AdjointPerformance" }
                                                                                         , { "not_clear", "Not" }
                                                                                         , { "title", "Interpolated zero rate(by Linear Interp.) differentiation performance with respect to yield rate" }
                                                                                         , { "ylabel", "Time (s)" }
                                                                                         , { "xlabel", "Number of yield rates" }
                                                                                         , { "line_box_width", "-5" } });
    cl::AdjointTestOutput outAdjoint("AdjointPiecwZeroSprTermStruct//LinInterpMultSpreads", { { "filename", "Adjoint" }
                                                                                        , { "not_clear", "Not" }
                                                                                        , { "title", "Interpolated zero rate (by Linear Interp.) adjoint differentiation with respect to yield rate" }
                                                                                        , { "xlabel", "Number of yield rates" }
                                                                                        , { "ylabel", "Adjoint calculations time (s)" } });
    cl::AdjointTestOutput outSize("AdjointPiecwZeroSprTermStruct//LinInterpMultSpreads", { { "filename", "TapeSize" }
                                                                                         , { "not_clear", "Not" }
                                                                                         , { "title", "Tape size dependence on number of yield rates" }
                                                                                         , { "xlabel", "Number of yield rates" }
                                                                                         , { "ylabel", "Memory(MB)" } });
    bool ok = true;

    for (Size rate_size = MIN_RATE_SIZE; rate_size <= MAX_RATE_SIZE; rate_size++)
    {
        pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        tapeSize[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        std::vector<double> rates_double = getDoubleVector(rate_size);
        std::vector<Integer> ts = getIntegerVector(rates_double.size() - 1);
        std::vector<Rate> rates(rates_double.size());

        for (Size i = 0; i < rates.size(); i++)
        {
            rates[i] = Real(rates_double[i]);
        }
        outPerform.log() << std::endl << "Start of tape recording " << currentTime() << std::endl;
        boost::timer timer;
        Independent(rates);
        CommonVars vars(rates, ts);

        std::vector<Handle<Quote> > spreads;
        boost::shared_ptr<SimpleQuote> spread1 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread2 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread3 = boost::make_shared<SimpleQuote>(0.035);
        boost::shared_ptr<SimpleQuote> spread4 = boost::make_shared<SimpleQuote>(0.04);
        spreads.push_back(Handle<Quote>(spread1));
        spreads.push_back(Handle<Quote>(spread2));
        spreads.push_back(Handle<Quote>(spread3));
        spreads.push_back(Handle<Quote>(spread4));

        std::vector<Date> spreadDates;
        spreadDates.push_back(vars.calendar.advance(vars.today, 90, Days));
        spreadDates.push_back(vars.calendar.advance(vars.today, 150, Days));
        spreadDates.push_back(vars.calendar.advance(vars.today, 30, Months));
        spreadDates.push_back(vars.calendar.advance(vars.today, 40, Months));

        Date interpolationDate = vars.calendar.advance(vars.today, 120, Days);

        boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
            boost::make_shared<PiecewiseZeroSpreadedTermStructure>(
            Handle<YieldTermStructure>(vars.termStructure),
            spreads, spreadDates);

        Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
        Rate interpolatedZeroRate = spreadedTermStructure->zeroRate(t, vars.compounding);
        Real expectedRate = vars.termStructure->zeroRate(t, vars.compounding) +
            spread1->value();
        outPerform.log() << "End of tape recording " << currentTime() << std::endl;
        pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ = timer.elapsed();
        std::vector<Rate> resultRates = { interpolatedZeroRate, expectedRate };
        cl::TapeFunction<double> f(rates, resultRates);
        tapeSize[rate_size - MIN_RATE_SIZE].memory_ = f.Memory();
        outPerform.log() << "Tape Memory in bytes: " << f.Memory() << std::endl;
        outPerform.log() << "Time for tape recording: " << pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ << std::endl;
        //Compute derivatives in forward and reverse mode
        std::vector<double> dy_forward(resultRates.size() * rates.size());
        std::vector<double> dy_reverce(resultRates.size() * rates.size());
        pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = std::min(gradReverse(f, dy_reverce, outPerform, false, false, repeatNumber)
            , gradForward(f, dy_forward, outPerform, false, false, repeatNumber));
        double h = 1.0e-6;
        std::vector<Real> dy_diff(resultRates.size() * rates.size());
        timer.restart();
        //Finite differences
        for (Size i = 0; i < rates.size(); i++)
        {
            rates[i] += h;
            CommonVars vars(rates, ts);

            std::vector<Date> spreadDates;
            spreadDates.push_back(vars.calendar.advance(vars.today, 90, Days));
            spreadDates.push_back(vars.calendar.advance(vars.today, 150, Days));
            spreadDates.push_back(vars.calendar.advance(vars.today, 30, Months));
            spreadDates.push_back(vars.calendar.advance(vars.today, 40, Months));

            Date interpolationDate = vars.calendar.advance(vars.today, 120, Days);

            boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
                boost::make_shared<PiecewiseZeroSpreadedTermStructure>(
                Handle<YieldTermStructure>(vars.termStructure),
                spreads, spreadDates);

            Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
            Rate interpolatedZeroRate_h = spreadedTermStructure->zeroRate(t, vars.compounding);

            Real expectedRate_h = vars.termStructure->zeroRate(t, vars.compounding) +
                spread1->value();

            dy_diff[i] = (interpolatedZeroRate_h - interpolatedZeroRate) / h;
            dy_diff[rates.size() + i] = (expectedRate_h - expectedRate) / h;
            rates[i] -= h;
        }
        pTime[rate_size - MIN_RATE_SIZE].timeAnalytical_ = timer.elapsed();
        adjTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_;
        adjTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_;
        ok &= checkWithFiniteDiff(dy_forward, dy_reverce, dy_diff, outPerform, 1e-7, 1e-7);
    }
    outPerform << pTime;
    outPerform.log() << "Linear Interpolation Multiple Spreads Performance Time plot has been generated " << currentTime() << std::endl;
    outAdjoint << adjTime;
    outPerform.log() << "Linear Interpolation Multiple Spreads Adjoint Performance plot has been generated " << currentTime() << std::endl;
    outSize << tapeSize;
    outPerform.log() << "Tape size dependence on number of yield rates plot has been generated " << currentTime() << std::endl;
    outputLinearInterpolationMultipleSpreads(outPerform);
    return ok;
}

// Method testLinearInterpolation()
// Testing linear interpolation between two dates.
// zeroRate is interpolated in a date point by the above mentioned method using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays) + spread quote.
// The output quantities (interpolatedZeroRate and expectedRate) are differentiated
// with respect to input yield rate vector values.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same plot.
// Linear Interpolation fits a line between the adjacent nodes, and returns the point on that line corresponding to the input x-value.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it throws an exception and returns false.
bool AdjointPiecewiseZeroSpreadedTermStructureTest::testLinearInterpolation() {

    BOOST_MESSAGE("Testing linear interpolation between two dates...");
    std::vector<PerformanceTime> pTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<AdjointTime> adjTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<TapeSize> tapeSize(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    cl::AdjointTestOutput outPerform("AdjointPiecwZeroSprTermStruct//LinInterp", { { "filename", "AdjointPerformance" }
                                                                              , { "not_clear", "Not" }
                                                                              , { "title", "Interpolated zero rate(by Linear Interp.) differentiation performance with respect to yield rate" }
                                                                              , { "ylabel", "Time (s)" }
                                                                              , { "xlabel", "Number of yield rates" }
                                                                              , { "line_box_width", "-5" } });
    cl::AdjointTestOutput outAdjoint("AdjointPiecwZeroSprTermStruct//LinInterp", { { "filename", "Adjoint" }
                                                                           , { "not_clear", "Not" }
                                                                           , { "title", "Interpolated zero rate (by Linear Interp.) adjoint differentiation with respect to yield rate" }
                                                                           , { "xlabel", "Number of yield rates" }
                                                                           , { "ylabel", "Adjoint calculations time (s)" } });
    cl::AdjointTestOutput outSize("AdjointPiecwZeroSprTermStruct//LinInterp", { { "filename", "TapeSize" }
                                                                              , { "not_clear", "Not" }
                                                                              , { "title", "Tape size dependence on number of yield rates" }
                                                                              , { "xlabel", "Number of yield rates" }
                                                                              , { "ylabel", "Memory(MB)" } });
    bool ok = true;

    for (Size rate_size = MIN_RATE_SIZE; rate_size <= MAX_RATE_SIZE; rate_size++)
    {
        pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        tapeSize[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        std::vector<double> rates_double = getDoubleVector(rate_size);
        std::vector<Integer> ts = getIntegerVector(rates_double.size() - 1);
        std::vector<Rate> rates(rates_double.size());

        for (Size i = 0; i < rates.size(); i++)
        {
            rates[i] = Real(rates_double[i]);
        }
        outPerform.log() << std::endl << "Start of tape recording " << currentTime() << std::endl;
        boost::timer timer;
        Independent(rates);
        CommonVars vars(rates, ts);

        std::vector<Handle<Quote> > spreads;
        boost::shared_ptr<SimpleQuote> spread1 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread2 = boost::make_shared<SimpleQuote>(0.03);
        spreads.push_back(Handle<Quote>(spread1));
        spreads.push_back(Handle<Quote>(spread2));

        std::vector<Date> spreadDates;
        spreadDates.push_back(vars.calendar.advance(vars.today, 100, Days));
        spreadDates.push_back(vars.calendar.advance(vars.today, 150, Days));

        Date interpolationDate = vars.calendar.advance(vars.today, 120, Days);

        boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
            boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<Linear> >(
            Handle<YieldTermStructure>(vars.termStructure),
            spreads, spreadDates);

        Date d0 = vars.calendar.advance(vars.today, 100, Days);
        Date d1 = vars.calendar.advance(vars.today, 150, Days);
        Date d2 = vars.calendar.advance(vars.today, 120, Days);

        Real m = (0.03 - 0.02) / vars.dayCount.yearFraction(d0, d1);
        Real expectedRate = m * vars.dayCount.yearFraction(d0, d2) + 0.054;

        Time t = vars.dayCount.yearFraction(vars.settlementDate, interpolationDate);
        Rate interpolatedZeroRate = spreadedTermStructure->zeroRate(t, vars.compounding);
        pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ = timer.elapsed();
        outPerform.log() << "End of tape recording " << currentTime() << std::endl;
        std::vector<Rate> resultRates = { interpolatedZeroRate, expectedRate };
        cl::TapeFunction<double> f(rates, resultRates);
        tapeSize[rate_size - MIN_RATE_SIZE].memory_ = f.Memory();
        outPerform.log() << "Tape Memory in bytes: " << f.Memory() << std::endl;
        outPerform.log() << "Time for tape recording: " << pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ << std::endl;
        //Compute derivatives in forward and reverse mode
        std::vector<double> dy_forward(resultRates.size() * rates.size());
        std::vector<double> dy_reverce(resultRates.size() * rates.size());
        pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = std::min(gradReverse(f, dy_reverce, outPerform, false, false, repeatNumber)
            , gradForward(f, dy_forward, outPerform, false, false, repeatNumber));
        double h = 1.0e-6;
        std::vector<Real> dy_diff(resultRates.size() * rates.size());
        timer.restart();
        //Finite differences
        for (Size i = 0; i < rates.size(); i++)
        {
            rates[i] += h;
            CommonVars vars(rates, ts);

            std::vector<Date> spreadDates;
            spreadDates.push_back(vars.calendar.advance(vars.today, 100, Days));
            spreadDates.push_back(vars.calendar.advance(vars.today, 150, Days));

            Date interpolationDate = vars.calendar.advance(vars.today, 120, Days);

            boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
                boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<Linear> >(
                Handle<YieldTermStructure>(vars.termStructure),
                spreads, spreadDates);

            Date d0 = vars.calendar.advance(vars.today, 100, Days);
            Date d1 = vars.calendar.advance(vars.today, 150, Days);
            Date d2 = vars.calendar.advance(vars.today, 120, Days);

            Real m = (0.03 - 0.02) / vars.dayCount.yearFraction(d0, d1);
            Real expectedRate_h = m * vars.dayCount.yearFraction(d0, d2) + 0.054;

            Time t = vars.dayCount.yearFraction(vars.settlementDate, interpolationDate);
            Rate interpolatedZeroRate_h = spreadedTermStructure->zeroRate(t, vars.compounding);

            dy_diff[i] = (interpolatedZeroRate_h - interpolatedZeroRate) / h;
            dy_diff[rates.size() + i] = (expectedRate_h - expectedRate) / h;
            rates[i] -= h;
        }
        pTime[rate_size - MIN_RATE_SIZE].timeAnalytical_ = timer.elapsed();
        adjTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_;
        adjTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_;
        ok &= checkWithFiniteDiff(dy_forward, dy_reverce, dy_diff, outPerform, 1e-7, 1e-7);
    }
    outPerform << pTime;
    outPerform.log() << "Linear Interpolation Performance Time plot has been generated " << currentTime() << std::endl;
    outAdjoint << adjTime;
    outPerform.log() << "Linear Interpolation Adjoint Performance plot has been generated " << currentTime() << std::endl;
    outSize << tapeSize;
    outPerform.log() << "Tape size dependence on number of yield rates plot has been generated " << currentTime() << std::endl;
    outputLinearInterpolation(outPerform);
    return ok;
}

// Method testForwardFlatInterpolation()
// Testing forward flat interpolation between two dates.
// zeroRate is interpolated in a date point by the above mentioned method using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays) + spread quote.
// The output quantities (interpolatedZeroRate and expectedRate) are differentiated
// with respect to input yield rate vector values.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same plot.
// Forward Flat Interpolation returns the value of the left node.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it throws an exception and returns false.
bool AdjointPiecewiseZeroSpreadedTermStructureTest::testForwardFlatInterpolation() {

    BOOST_MESSAGE("Testing forward flat interpolation between two dates...");
    std::vector<PerformanceTime> pTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<AdjointTime> adjTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<TapeSize> tapeSize(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    cl::AdjointTestOutput outPerform("AdjointPiecwZeroSprTermStruct//ForwardFlatInterp", { { "filename", "AdjointPerformance" }
                                                                                      , { "not_clear", "Not" }
                                                                                      , { "title", "Interpolated zero rate(by Forward Flat Interp.) differentiation performance with respect to yield rate" }
                                                                                      , { "ylabel", "Time (s)" }
                                                                                      , { "xlabel", "Number of yield rates" }
                                                                                      , { "line_box_width", "-5" } });
    cl::AdjointTestOutput outAdjoint("AdjointPiecwZeroSprTermStruct//ForwardFlatInterp", { { "filename", "Adjoint" }
                                                                                     , { "not_clear", "Not" }
                                                                                     , { "title", "Interpolated zero rate (by Forward Flat Interp.) adjoint differentiation with respect to yield rate" }
                                                                                     , { "xlabel", "Number of yield rates" }
                                                                                     , { "ylabel", "Adjoint calculations time (s)" } });
    cl::AdjointTestOutput outSize("AdjointPiecwZeroSprTermStruct//ForwardFlatInterp", { { "filename", "TapeSize" }
                                                                                      , { "not_clear", "Not" }
                                                                                      , { "title", "Tape size dependence on number of yield rates" }
                                                                                      , { "xlabel", "Number of yield rates" }
                                                                                      , { "ylabel", "Memory(MB)" } });
    bool ok = true;

    for (Size rate_size = MIN_RATE_SIZE; rate_size <= MAX_RATE_SIZE; rate_size++)
    {
        pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        tapeSize[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        std::vector<double> rates_double = getDoubleVector(rate_size);
        std::vector<Integer> ts = getIntegerVector(rates_double.size() - 1);
        std::vector<Rate> rates(rates_double.size());

        for (Size i = 0; i < rates.size(); i++)
        {
            rates[i] = Real(rates_double[i]);
        }
        outPerform.log() << std::endl << "Start of tape recording " << currentTime() << std::endl;
        boost::timer timer;
        Independent(rates);
        CommonVars vars(rates, ts);

        std::vector<Handle<Quote> > spreads;
        boost::shared_ptr<SimpleQuote> spread1 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread2 = boost::make_shared<SimpleQuote>(0.03);
        spreads.push_back(Handle<Quote>(spread1));
        spreads.push_back(Handle<Quote>(spread2));

        std::vector<Date> spreadDates;
        spreadDates.push_back(vars.calendar.advance(vars.today, 75, Days));
        spreadDates.push_back(vars.calendar.advance(vars.today, 260, Days));

        Date interpolationDate = vars.calendar.advance(vars.today, 100, Days);

        boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
            boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<ForwardFlat> >(
            Handle<YieldTermStructure>(vars.termStructure),
            spreads, spreadDates);

        Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
        Rate interpolatedZeroRate = spreadedTermStructure->zeroRate(t, vars.compounding);
        Real expectedRate = vars.termStructure->zeroRate(t, vars.compounding) +
            spread1->value();
        outPerform.log() << "End of tape recording " << currentTime() << std::endl;
        pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ = timer.elapsed();
        std::vector<Rate> resultRates = { interpolatedZeroRate, expectedRate };
        cl::TapeFunction<double> f(rates, resultRates);
        tapeSize[rate_size - MIN_RATE_SIZE].memory_ = f.Memory();
        outPerform.log() << "Tape Memory in bytes: " << f.Memory() << std::endl;
        outPerform.log() << "Time for tape recording: " << pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ << std::endl;
        //Compute derivatives in forward and reverse mode
        std::vector<double> dy_forward(resultRates.size() * rates.size());
        std::vector<double> dy_reverce(resultRates.size() * rates.size());
        pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = std::min(gradReverse(f, dy_reverce, outPerform, false, false, repeatNumber)
            , gradForward(f, dy_forward, outPerform, false, false, repeatNumber));
        double h = 1.0e-6;
        std::vector<Real> dy_diff(resultRates.size() * rates.size());
        timer.restart();
        //Finite differences
        for (Size i = 0; i < rates.size(); i++)
        {
            rates[i] += h;
            CommonVars vars(rates, ts);

            std::vector<Date> spreadDates;
            spreadDates.push_back(vars.calendar.advance(vars.today, 75, Days));
            spreadDates.push_back(vars.calendar.advance(vars.today, 260, Days));

            Date interpolationDate = vars.calendar.advance(vars.today, 100, Days);

            boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
                boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<ForwardFlat> >(
                Handle<YieldTermStructure>(vars.termStructure),
                spreads, spreadDates);

            Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
            Rate interpolatedZeroRate_h = spreadedTermStructure->zeroRate(t, vars.compounding);


            Real expectedRate_h = vars.termStructure->zeroRate(t, vars.compounding) +
                spread1->value();

            dy_diff[i] = (interpolatedZeroRate_h - interpolatedZeroRate) / h;
            dy_diff[rates.size() + i] = (expectedRate_h - expectedRate) / h;
            rates[i] -= h;
        }
        pTime[rate_size - MIN_RATE_SIZE].timeAnalytical_ = timer.elapsed();
        adjTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_;
        adjTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_;
        ok &= checkWithFiniteDiff(dy_forward, dy_reverce, dy_diff, outPerform, 1e-7, 1e-7);
    }
    outPerform << pTime;
    outPerform.log() << "Forward Flat Interpolation Performance Time plot has been generated " << currentTime() << std::endl;
    outAdjoint << adjTime;
    outPerform.log() << "Forward Flat Interpolation Adjoint Performance plot has been generated " << currentTime() << std::endl;
    outSize << tapeSize;
    outPerform.log() << "Tape size dependence on number of yield rates plot has been generated " << currentTime() << std::endl;
    outputForwardFlatInterpolation(outPerform);
    return ok;
}

// Method testBackwardFlatInterpolation()
// Testing backward flat interpolation between two dates.
// zeroRate is interpolated in a date point by the above mentioned method using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays) + spread quote.
// The output quantities (interpolatedZeroRate and expectedRate) are differentiated
// with respect to input yield rate vector values.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same plot.
// Backward Flat Interpolation returns the value of the right node.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it throws an exception and returns false.
bool AdjointPiecewiseZeroSpreadedTermStructureTest::testBackwardFlatInterpolation() {

    BOOST_MESSAGE("Testing backward flat interpolation between two dates...");
    std::vector<PerformanceTime> pTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<AdjointTime> adjTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<TapeSize> tapeSize(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    cl::AdjointTestOutput outPerform("AdjointPiecwZeroSprTermStruct//BackwardFlatInterp", { { "filename", "AdjointPerformance" }
                                                                                       , { "not_clear", "Not" }
                                                                                       , { "title", "Interpolated zero rate(by Backward Flat Interp.) differentiation performance with respect to yield rate" }
                                                                                       , { "ylabel", "Time (s)" }
                                                                                       , { "xlabel", "Number of yield rates" }
                                                                                       , { "line_box_width", "-5" } });
    cl::AdjointTestOutput outAdjoint("AdjointPiecwZeroSprTermStruct//BackwardFlatInterp", { { "filename", "Adjoint" }
                                                                                      , { "not_clear", "Not" }
                                                                                      , { "title", "Interpolated zero rate (by Backward Flat Interp.) adjoint differentiation with respect to yield rate" }
                                                                                      , { "xlabel", "Number of yield rates" }
                                                                                      , { "ylabel", "Adjoint calculations time (s)" } });
    cl::AdjointTestOutput outSize("AdjointPiecwZeroSprTermStruct//BackwardFlatInterp", { { "filename", "TapeSize" }
                                                                                       , { "not_clear", "Not" }
                                                                                       , { "title", "Tape size dependence on number of yield rates" }
                                                                                       , { "xlabel", "Number of yield rates" }
                                                                                       , { "ylabel", "Memory(MB)" } });
    bool ok = true;

    for (Size rate_size = MIN_RATE_SIZE; rate_size <= MAX_RATE_SIZE; rate_size++)
    {
        pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        tapeSize[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        std::vector<double> rates_double = getDoubleVector(rate_size);
        std::vector<Integer> ts = getIntegerVector(rates_double.size() - 1);
        std::vector<Rate> rates(rates_double.size());

        for (Size i = 0; i < rates.size(); i++)
        {
            rates[i] = Real(rates_double[i]);
        }
        outPerform.log() << std::endl << "Start of tape recording " << currentTime() << std::endl;
        boost::timer timer;
        Independent(rates);
        CommonVars vars(rates, ts);

        std::vector<Handle<Quote> > spreads;
        boost::shared_ptr<SimpleQuote> spread1 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread2 = boost::make_shared<SimpleQuote>(0.03);
        boost::shared_ptr<SimpleQuote> spread3 = boost::make_shared<SimpleQuote>(0.04);
        spreads.push_back(Handle<Quote>(spread1));
        spreads.push_back(Handle<Quote>(spread2));
        spreads.push_back(Handle<Quote>(spread3));

        std::vector<Date> spreadDates;
        spreadDates.push_back(vars.calendar.advance(vars.today, 100, Days));
        spreadDates.push_back(vars.calendar.advance(vars.today, 200, Days));
        spreadDates.push_back(vars.calendar.advance(vars.today, 300, Days));

        Date interpolationDate = vars.calendar.advance(vars.today, 110, Days);

        boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
            boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<BackwardFlat> >(
            Handle<YieldTermStructure>(vars.termStructure),
            spreads, spreadDates);

        Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
        Rate interpolatedZeroRate = spreadedTermStructure->zeroRate(t, vars.compounding);
        Real expectedRate = vars.termStructure->zeroRate(t, vars.compounding) +
            spread2->value();
        outPerform.log() << "End of tape recording " << currentTime() << std::endl;
        pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ = timer.elapsed();
        std::vector<Rate> resultRates = { interpolatedZeroRate, expectedRate };
        cl::TapeFunction<double> f(rates, resultRates);
        tapeSize[rate_size - MIN_RATE_SIZE].memory_ = f.Memory();
        outPerform.log() << "Tape Memory in bytes: " << f.Memory() << std::endl;
        outPerform.log() << "Time for tape recording: " << pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ << std::endl;
        //Compute derivatives in forward and reverse mode
        std::vector<double> dy_forward(resultRates.size() * rates.size());
        std::vector<double> dy_reverce(resultRates.size() * rates.size());
        pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = std::min(gradReverse(f, dy_reverce, outPerform, false, false, repeatNumber)
            , gradForward(f, dy_forward, outPerform, false, false, repeatNumber));
        double h = 1.0e-6;
        std::vector<Real> dy_diff(resultRates.size() * rates.size());
        timer.restart();
        //Finite differences
        for (Size i = 0; i < rates.size(); i++)
        {
            rates[i] += h;
            CommonVars vars(rates, ts);

            std::vector<Date> spreadDates;
            spreadDates.push_back(vars.calendar.advance(vars.today, 100, Days));
            spreadDates.push_back(vars.calendar.advance(vars.today, 200, Days));
            spreadDates.push_back(vars.calendar.advance(vars.today, 300, Days));

            Date interpolationDate = vars.calendar.advance(vars.today, 110, Days);

            boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
                boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<BackwardFlat> >(
                Handle<YieldTermStructure>(vars.termStructure),
                spreads, spreadDates);

            Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
            Rate interpolatedZeroRate_h = spreadedTermStructure->zeroRate(t, vars.compounding);


            Real expectedRate_h = vars.termStructure->zeroRate(t, vars.compounding) +
                spread2->value();

            dy_diff[i] = (interpolatedZeroRate_h - interpolatedZeroRate) / h;
            dy_diff[rates.size() + i] = (expectedRate_h - expectedRate) / h;
            rates[i] -= h;
        }
        pTime[rate_size - MIN_RATE_SIZE].timeAnalytical_ = timer.elapsed();
        adjTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_;
        adjTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_;
        ok &= checkWithFiniteDiff(dy_forward, dy_reverce, dy_diff, outPerform, 1e-7, 1e-7);
    }
    outPerform << pTime;
    outPerform.log() << "Backward Flat Interpolation Performance Time plot has been generated " << currentTime() << std::endl;
    outAdjoint << adjTime;
    outPerform.log() << "Backward Flat Interpolation Adjoint Performance plot has been generated " << currentTime() << std::endl;
    outSize << tapeSize;
    outPerform.log() << "Tape size dependence on number of yield rates plot has been generated " << currentTime() << std::endl;
    outputBackwardFlatInterpolation(outPerform);
    return ok;
}

// Method testDefaultInterpolation()
// Testing default interpolation between two dates.
// Linear Interpolation is used as a Default Interpolation.
// zeroRate is interpolated in a date point by the above mentioned method using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays) + spread quote.
// The output quantities (interpolatedZeroRate and expectedRate) are differentiated
// with respect to input yield rate vector values.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same plot.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it throws an exception and returns false.
bool AdjointPiecewiseZeroSpreadedTermStructureTest::testDefaultInterpolation() {

    BOOST_MESSAGE("Testing default interpolation between two dates...");
    std::vector<PerformanceTime> pTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<AdjointTime> adjTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<TapeSize> tapeSize(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    cl::AdjointTestOutput outPerform("AdjointPiecwZeroSprTermStruct//DefaultInterp", { { "filename", "AdjointPerformance" }
                                                                                  , { "not_clear", "Not" }
                                                                                  , { "title", "Interpolated zero rate(by Default (Linear) Interp.) differentiation performance with respect to yield rate" }
                                                                                  , { "ylabel", "Time (s)" }
                                                                                  , { "xlabel", "Number of yield rates" }
                                                                                  , { "line_box_width", "-5" } });
    cl::AdjointTestOutput outAdjoint("AdjointPiecwZeroSprTermStruct//DefaultInterp", { { "filename", "Adjoint" }
                                                                                 , { "not_clear", "Not" }
                                                                                 , { "title", "Interpolated zero rate (by Default (Linear) Interp.) adjoint differentiation with respect to yield rate" }
                                                                                 , { "xlabel", "Number of yield rates" }
                                                                                 , { "ylabel", "Adjoint calculations time (s)" } });
    cl::AdjointTestOutput outSize("AdjointPiecwZeroSprTermStruct//DefaultInterp", { { "filename", "TapeSize" }
                                                                                  , { "not_clear", "Not" }
                                                                                  , { "title", "Tape size dependence on number of yield rates" }
                                                                                  , { "xlabel", "Number of yield rates" }
                                                                                  , { "ylabel", "Memory(MB)" } });
    bool ok = true;

    for (Size rate_size = MIN_RATE_SIZE; rate_size <= MAX_RATE_SIZE; rate_size++)
    {
        pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        tapeSize[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        std::vector<double> rates_double = getDoubleVector(rate_size);
        std::vector<Integer> ts = getIntegerVector(rates_double.size() - 1);
        std::vector<Rate> rates(rates_double.size());

        for (Size i = 0; i < rates.size(); i++)
        {
            rates[i] = Real(rates_double[i]);
        }
        outPerform.log() << std::endl << "Start of tape recording " << currentTime() << std::endl;
        boost::timer timer;
        Independent(rates);
        CommonVars vars(rates, ts);

        std::vector<Handle<Quote> > spreads;
        boost::shared_ptr<SimpleQuote> spread1 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread2 = boost::make_shared<SimpleQuote>(0.02);
        spreads.push_back(Handle<Quote>(spread1));
        spreads.push_back(Handle<Quote>(spread2));

        std::vector<Date> spreadDates;
        spreadDates.push_back(vars.calendar.advance(vars.today, 75, Days));
        spreadDates.push_back(vars.calendar.advance(vars.today, 160, Days));

        Date interpolationDate = vars.calendar.advance(vars.today, 100, Days);

        boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
            boost::make_shared<PiecewiseZeroSpreadedTermStructure>(
            Handle<YieldTermStructure>(vars.termStructure),
            spreads, spreadDates);

        Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
        Rate interpolatedZeroRate = spreadedTermStructure->zeroRate(t, vars.compounding);
        Real expectedRate = vars.termStructure->zeroRate(t, vars.compounding) +
            spread1->value();
        outPerform.log() << "End of tape recording " << currentTime() << std::endl;
        pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ = timer.elapsed();
        std::vector<Rate> resultRates = { interpolatedZeroRate, expectedRate };
        cl::TapeFunction<double> f(rates, resultRates);
        tapeSize[rate_size - MIN_RATE_SIZE].memory_ = f.Memory();
        outPerform.log() << "Tape Memory in bytes: " << f.Memory() << std::endl;
        outPerform.log() << "Time for tape recording: " << pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ << std::endl;
        //Compute derivatives in forward and reverse mode
        std::vector<double> dy_forward(resultRates.size() * rates.size());
        std::vector<double> dy_reverce(resultRates.size() * rates.size());
        pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = std::min(gradReverse(f, dy_reverce, outPerform, false, false, repeatNumber)
            , gradForward(f, dy_forward, outPerform, false, false, repeatNumber));
        double h = 1.0e-6;
        std::vector<Real> dy_diff(resultRates.size() * rates.size());
        timer.restart();
        //Finite differences
        for (Size i = 0; i < rates.size(); i++)
        {
            rates[i] += h;
            CommonVars vars(rates, ts);

            std::vector<Date> spreadDates;
            spreadDates.push_back(vars.calendar.advance(vars.today, 75, Days));
            spreadDates.push_back(vars.calendar.advance(vars.today, 160, Days));

            Date interpolationDate = vars.calendar.advance(vars.today, 100, Days);

            boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
                boost::make_shared<PiecewiseZeroSpreadedTermStructure>(
                Handle<YieldTermStructure>(vars.termStructure),
                spreads, spreadDates);

            Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
            Rate interpolatedZeroRate_h = spreadedTermStructure->zeroRate(t, vars.compounding);


            Real expectedRate_h = vars.termStructure->zeroRate(t, vars.compounding) +
                spread1->value();

            dy_diff[i] = (interpolatedZeroRate_h - interpolatedZeroRate) / h;
            dy_diff[rates.size() + i] = (expectedRate_h - expectedRate) / h;
            rates[i] -= h;
        }
        pTime[rate_size - MIN_RATE_SIZE].timeAnalytical_ = timer.elapsed();
        adjTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_;
        adjTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_;
        ok &= checkWithFiniteDiff(dy_forward, dy_reverce, dy_diff, outPerform, 1e-7, 1e-7);
    }
    outPerform << pTime;
    outPerform.log() << "Default Interpolation Performance Time plot has been generated " << currentTime() << std::endl;
    outAdjoint << adjTime;
    outPerform.log() << "Default Interpolation Adjoint Performance plot has been generated " << currentTime() << std::endl;
    outSize << tapeSize;
    outPerform.log() << "Tape size dependence on number of yield rates plot has been generated " << currentTime() << std::endl;
    outputDefaultInterpolation(outPerform);
    return ok;
}

// Method testSetInterpolationFactory()
// Testing factory constructor with additional parameters.
// As a parameter for interpolation factory used Cubic Spline Interpolation -
// interpolation with special piecewise cubic polynomial.
// zeroRate is interpolated in a date point by the above mentioned method using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays) + spread quote.
// The output quantities (interpolatedZeroRate and expectedRate) are differentiated
// with respect to input yield rate vector values.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same plot.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it throws an exception and returns false.
bool AdjointPiecewiseZeroSpreadedTermStructureTest::testSetInterpolationFactory() {

    BOOST_MESSAGE("Testing factory constructor with additional parameters...");
    std::vector<PerformanceTime> pTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<AdjointTime> adjTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<TapeSize> tapeSize(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    cl::AdjointTestOutput outPerform("AdjointPiecwZeroSprTermStruct//SetInterpFactory", { { "filename", "AdjointPerformance" }
                                                                                     , { "not_clear", "Not" }
                                                                                     , { "title", "Interpolated zero rate(by Cubic Spline Interp.) differentiation performance with respect to yield rate" }
                                                                                     , { "ylabel", "Time (s)" }
                                                                                     , { "xlabel", "Number of yield rates" }
                                                                                     , { "line_box_width", "-5" } });
    cl::AdjointTestOutput outAdjoint("AdjointPiecwZeroSprTermStruct//SetInterpFactory", { { "filename", "Adjoint" }
                                                                                    , { "not_clear", "Not" }
                                                                                    , { "title", "Interpolated zero rate (by Cubic Spline Interp.) adjoint differentiation with respect to yield rate" }
                                                                                    , { "xlabel", "Number of yield rates" }
                                                                                    , { "ylabel", "Adjoint calculations time (s)" } });
    cl::AdjointTestOutput outSize("AdjointPiecwZeroSprTermStruct//SetInterpFactory", { { "filename", "TapeSize" }
                                                                                     , { "not_clear", "Not" }
                                                                                     , { "title", "Tape size dependence on number of yield rates" }
                                                                                     , { "xlabel", "Number of yield rates" }
                                                                                     , { "ylabel", "Memory(MB)" } });
    bool ok = true;

    for (Size rate_size = MIN_RATE_SIZE; rate_size <= MAX_RATE_SIZE; rate_size++)
    {
        pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        tapeSize[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        std::vector<double> rates_double = getDoubleVector(rate_size);
        std::vector<Integer> ts = getIntegerVector(rates_double.size() - 1);
        std::vector<Rate> rates(rates_double.size());

        for (Size i = 0; i < rates.size(); i++)
        {
            rates[i] = Real(rates_double[i]);
        }
        outPerform.log() << std::endl << "Start of tape recording " << currentTime() << std::endl;
        boost::timer timer;
        Independent(rates);
        CommonVars vars(rates, ts);

        std::vector<Handle<Quote> > spreads;
        boost::shared_ptr<SimpleQuote> spread1 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread2 = boost::make_shared<SimpleQuote>(0.03);
        boost::shared_ptr<SimpleQuote> spread3 = boost::make_shared<SimpleQuote>(0.01);
        spreads.push_back(Handle<Quote>(spread1));
        spreads.push_back(Handle<Quote>(spread2));
        spreads.push_back(Handle<Quote>(spread3));

        std::vector<Date> spreadDates;
        spreadDates.push_back(vars.calendar.advance(vars.today, 8, Months));
        spreadDates.push_back(vars.calendar.advance(vars.today, 15, Months));
        spreadDates.push_back(vars.calendar.advance(vars.today, 25, Months));

        Date interpolationDate = vars.calendar.advance(vars.today, 11, Months);

        boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure;

        Frequency freq = NoFrequency;

        Cubic factory;
        factory = Cubic(CubicInterpolation::Spline, false);

        spreadedTermStructure =
            boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<Cubic> >(
            Handle<YieldTermStructure>(vars.termStructure),
            spreads, spreadDates, vars.compounding,
            freq, vars.dayCount, factory);

        Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
        Rate interpolatedZeroRate = spreadedTermStructure->zeroRate(t, vars.compounding);
        Real expectedRate = vars.termStructure->zeroRate(t, vars.compounding) +
            Real(0.026065770863);
        outPerform.log() << "End of tape recording " << currentTime() << std::endl;
        pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ = timer.elapsed();
        std::vector<Rate> resultRates = { interpolatedZeroRate, expectedRate };
        cl::TapeFunction<double> f(rates, resultRates);
        tapeSize[rate_size - MIN_RATE_SIZE].memory_ = f.Memory();
        outPerform.log() << "Tape Memory in bytes: " << f.Memory() << std::endl;
        outPerform.log() << "Time for tape recording: " << pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ << std::endl;
        //Compute derivatives in forward and reverse mode
        std::vector<double> dy_forward(resultRates.size() * rates.size());
        std::vector<double> dy_reverce(resultRates.size() * rates.size());
        pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = std::min(gradReverse(f, dy_reverce, outPerform, false, false, repeatNumber)
            , gradForward(f, dy_forward, outPerform, false, false, repeatNumber));
        double h = 1.0e-6;
        std::vector<Real> dy_diff(resultRates.size() * rates.size());
        timer.restart();
        //Finite differences
        for (Size i = 0; i < rates.size(); i++)
        {
            rates[i] += h;
            CommonVars vars(rates, ts);

            std::vector<Date> spreadDates;
            spreadDates.push_back(vars.calendar.advance(vars.today, 8, Months));
            spreadDates.push_back(vars.calendar.advance(vars.today, 15, Months));
            spreadDates.push_back(vars.calendar.advance(vars.today, 25, Months));

            Date interpolationDate = vars.calendar.advance(vars.today, 11, Months);

            boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure;

            Frequency freq = NoFrequency;

            Cubic factory;
            factory = Cubic(CubicInterpolation::Spline, false);

            spreadedTermStructure =
                boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<Cubic> >(
                Handle<YieldTermStructure>(vars.termStructure),
                spreads, spreadDates, vars.compounding,
                freq, vars.dayCount, factory);

            Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
            Rate interpolatedZeroRate_h = spreadedTermStructure->zeroRate(t, vars.compounding);


            Real expectedRate_h = vars.termStructure->zeroRate(t, vars.compounding) +
                Real(0.026065770863);

            dy_diff[i] = (interpolatedZeroRate_h - interpolatedZeroRate) / h;
            dy_diff[rates.size() + i] = (expectedRate_h - expectedRate) / h;
            rates[i] -= h;
        }
        pTime[rate_size - MIN_RATE_SIZE].timeAnalytical_ = timer.elapsed();
        adjTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_;
        adjTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_;
        ok &= checkWithFiniteDiff(dy_forward, dy_reverce, dy_diff, outPerform, 1e-7, 1e-7);
    }
    outPerform << pTime;
    outPerform.log() << "Set Interpolation Factory Performance Time plot has been generated " << currentTime() << std::endl;
    outAdjoint << adjTime;
    outPerform.log() << "Set Interpolation Factory Adjoint Performance plot has been generated " << currentTime() << std::endl;
    outSize << tapeSize;
    outPerform.log() << "Tape size dependence on number of yield rates plot has been generated " << currentTime() << std::endl;
    outputSetInterpolationFactory(outPerform);
    return ok;
}

// Method testQuoteChanging()
// Testing quote update.
// Backward Flat Interpolation with quote changing.
// zeroRate is interpolated in a date point by the above mentioned method using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays) + spread quote.
// The output quantities (interpolatedZeroRate and expectedRate) are differentiated
// with respect to input yield rate vector values.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same plot.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it throws an exception and returns false.
bool AdjointPiecewiseZeroSpreadedTermStructureTest::testQuoteChanging() {

    BOOST_MESSAGE("Testing quote update...");
    std::vector<PerformanceTime> pTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<AdjointTime> adjTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<TapeSize> tapeSize(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    cl::AdjointTestOutput outPerform("AdjointPiecwZeroSprTermStruct//QuoteChanging", { { "filename", "AdjointPerformance" }
                                                                                  , { "not_clear", "Not" }
                                                                                  , { "title", "Interpolated zero rate(by Backward Flat Interp. with quote changing) diff. perf. with respect to yield rate" }
                                                                                  , { "ylabel", "Time (s)" }
                                                                                  , { "xlabel", "Number of yield rates" }
                                                                                  , { "line_box_width", "-5" } });
    cl::AdjointTestOutput outAdjoint("AdjointPiecwZeroSprTermStruct//QuoteChanging", { { "filename", "Adjoint" }
                                                                                 , { "not_clear", "Not" }
                                                                                 , { "title", "Interpolated zero rate (by Backward Flat Interp. with quote changing) adjoint differentiation with respect to yield rate" }
                                                                                 , { "xlabel", "Number of yield rates" }
                                                                                 , { "ylabel", "Adjoint calculations time (s)" } });
    cl::AdjointTestOutput outSize("AdjointPiecwZeroSprTermStruct//QuoteChanging", { { "filename", "TapeSize" }
                                                                                  , { "not_clear", "Not" }
                                                                                  , { "title", "Tape size dependence on number of yield rates" }
                                                                                  , { "xlabel", "Number of yield rates" }
                                                                                  , { "ylabel", "Memory(MB)" } });
    bool ok = true;

    for (Size rate_size = MIN_RATE_SIZE; rate_size <= MAX_RATE_SIZE; rate_size++)
    {
        pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        tapeSize[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        std::vector<double> rates_double = getDoubleVector(rate_size);
        std::vector<Integer> ts = getIntegerVector(rates_double.size() - 1);
        std::vector<Rate> rates(rates_double.size());

        for (Size i = 0; i < rates.size(); i++)
        {
            rates[i] = Real(rates_double[i]);
        }
        outPerform.log() << std::endl << "Start of tape recording " << currentTime() << std::endl;
        boost::timer timer;
        Independent(rates);
        CommonVars vars(rates, ts);

        std::vector<Handle<Quote> > spreads;
        boost::shared_ptr<SimpleQuote> spread1 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread2 = boost::make_shared<SimpleQuote>(0.03);
        spreads.push_back(Handle<Quote>(spread1));
        spreads.push_back(Handle<Quote>(spread2));

        std::vector<Date> spreadDates;
        spreadDates.push_back(vars.calendar.advance(vars.today, 100, Days));
        spreadDates.push_back(vars.calendar.advance(vars.today, 150, Days));

        Date interpolationDate = vars.calendar.advance(vars.today, 120, Days);

        boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
            boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<BackwardFlat> >(
            Handle<YieldTermStructure>(vars.termStructure),
            spreads, spreadDates);

        Time t = vars.dayCount.yearFraction(vars.settlementDate, interpolationDate);
        Rate interpolatedZeroRate = spreadedTermStructure->zeroRate(t, vars.compounding);

        Real expectedRate = vars.termStructure->zeroRate(t, vars.compounding) +
            Real(0.03);

        spread2->setValue(0.025);

        interpolatedZeroRate = spreadedTermStructure->zeroRate(t, vars.compounding);
        expectedRate = vars.termStructure->zeroRate(t, vars.compounding) +
            Real(0.025);
        outPerform.log() << "End of tape recording " << currentTime() << std::endl;
        pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ = timer.elapsed();
        std::vector<Rate> resultRates = { interpolatedZeroRate, expectedRate };
        cl::TapeFunction<double> f(rates, resultRates);
        tapeSize[rate_size - MIN_RATE_SIZE].memory_ = f.Memory();
        outPerform.log() << "Tape Memory in bytes: " << f.Memory() << std::endl;
        outPerform.log() << "Time for tape recording: " << pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ << std::endl;
        //Compute derivatives in forward and reverse mode
        std::vector<double> dy_forward(resultRates.size() * rates.size());
        std::vector<double> dy_reverce(resultRates.size() * rates.size());
        pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = std::min(gradReverse(f, dy_reverce, outPerform, false, false, repeatNumber)
            , gradForward(f, dy_forward, outPerform, false, false, repeatNumber));
        double h = 1.0e-6;
        std::vector<Real> dy_diff(resultRates.size() * rates.size());
        timer.restart();
        //Finite differences
        for (Size i = 0; i < rates.size(); i++)
        {
            rates[i] += h;
            CommonVars vars(rates, ts);

            std::vector<Date> spreadDates;
            spreadDates.push_back(vars.calendar.advance(vars.today, 100, Days));
            spreadDates.push_back(vars.calendar.advance(vars.today, 150, Days));

            Date interpolationDate = vars.calendar.advance(vars.today, 120, Days);

            boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
                boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<BackwardFlat> >(
                Handle<YieldTermStructure>(vars.termStructure),
                spreads, spreadDates);

            Time t = vars.dayCount.yearFraction(vars.settlementDate, interpolationDate);
            Rate interpolatedZeroRate_h = spreadedTermStructure->zeroRate(t, vars.compounding);

            Real expectedRate_h = vars.termStructure->zeroRate(t, vars.compounding) +
                Real(0.025);

            dy_diff[i] = (interpolatedZeroRate_h - interpolatedZeroRate) / h;
            dy_diff[rates.size() + i] = (expectedRate_h - expectedRate) / h;
            rates[i] -= h;
        }
        pTime[rate_size - MIN_RATE_SIZE].timeAnalytical_ = timer.elapsed();
        adjTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_;
        adjTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_;
        ok &= checkWithFiniteDiff(dy_forward, dy_reverce, dy_diff, outPerform, 1e-7, 1e-7);
    }
    outPerform << pTime;
    outPerform.log() << "Quote Changing Performance Time plot has been generated " << currentTime() << std::endl;
    outAdjoint << adjTime;
    outPerform.log() << "Quote Changing Adjoint Performance plot has been generated " << currentTime() << std::endl;
    outSize << tapeSize;
    outPerform.log() << "Tape size dependence on number of yield rates plot has been generated " << currentTime() << std::endl;
    outputQuoteChanging(outPerform);
    return ok;
}

// Function outputFlatInterpolationLeft(cl::AdjointTestOutput & logOut)
// Generating Zero Rate dependence on Yield Rate Dependence By Flat Interpolation Left plot.
// zeroRate is interpolated in a date point by the above mentioned method using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays) + spread quote.
// Parameter: logOut is reference to struct with log stream.
void outputFlatInterpolationLeft(cl::AdjointTestOutput & logOut)
{
    Size steps = Size((MAX_RATE_RANGE - MIN_RATE_RANGE) / RATE_STEP) + 1;
    std::vector<Rate> rates = { MIN_RATE_RANGE, 0.033};
    std::vector<RateDependence> rDep(steps);

    for (Size i = 0; i < steps; rates[0] += RATE_STEP, i++)
    {
        CommonVars vars(rates, getIntegerVector(rates.size()));
        std::vector<Handle<Quote> > spreads;
        boost::shared_ptr<SimpleQuote> spread1 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread2 = boost::make_shared<SimpleQuote>(0.03);
        spreads.push_back(Handle<Quote>(spread1));
        spreads.push_back(Handle<Quote>(spread2));

        std::vector<Date> spreadDates;
        spreadDates.push_back(vars.calendar.advance(vars.today, 8, Months));
        spreadDates.push_back(vars.calendar.advance(vars.today, 15, Months));

        Date interpolationDate = vars.today;

        boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
            boost::make_shared<PiecewiseZeroSpreadedTermStructure>(
            Handle<YieldTermStructure>(vars.termStructure),
            spreads, spreadDates);

        Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
        Rate interpolatedZeroRate = spreadedTermStructure->zeroRate(t, vars.compounding);


        Real expectedRate = vars.termStructure->zeroRate(t, vars.compounding) +
            spread1->value();

        rDep[i].inputRate_ = rates[0];
        rDep[i].interpolatedZeroRate_ = interpolatedZeroRate;
    }
    cl::AdjointTestOutput out("AdjointPiecwZeroSprTermStruct//FlatInterpLeft//output", { { "filename", "zeroRateOnRateByFlatInterpLeft" }
                                                                                       , { "not_clear", "Not" }
                                                                                       , { "title", "Zero Rate dependence on Yield Rate By Flat Interpolation Left" }
                                                                                       , { "ylabel", "Zero Rate" } });
    out << rDep;
    logOut.log() << "Zero Rate dependence on Yield Rate By Flat Interpolation Left plot has been generated " << currentTime() << std::endl;
}

// Function outputFlatInterpolationRight(cl::AdjointTestOutput & logOut)
// Generating Zero Rate dependence on Yield Rate Dependence By Flat Interpolation Right plot.
// zeroRate is interpolated in a date point by the above mentioned method using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays) + spread quote.
// Parameter: logOut is reference to struct with log stream.
void outputFlatInterpolationRight(cl::AdjointTestOutput & logOut)
{
    Size steps = Size((MAX_RATE_RANGE - MIN_RATE_RANGE) / RATE_STEP) + 1;
    std::vector<Rate> rates = { MIN_RATE_RANGE, 0.033 };
    std::vector<RateDependence> rDep(steps);

    for (Size i = 0; i < steps; rates[0] += RATE_STEP, i++)
    {
        CommonVars vars(rates, getIntegerVector(rates.size()));
        std::vector<Handle<Quote> > spreads;
        boost::shared_ptr<SimpleQuote> spread1 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread2 = boost::make_shared<SimpleQuote>(0.03);
        spreads.push_back(Handle<Quote>(spread1));
        spreads.push_back(Handle<Quote>(spread2));

        std::vector<Date> spreadDates;
        spreadDates.push_back(vars.calendar.advance(vars.today, 8, Months));
        spreadDates.push_back(vars.calendar.advance(vars.today, 15, Months));

        Date interpolationDate = vars.today;

        boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
            boost::make_shared<PiecewiseZeroSpreadedTermStructure>(
            Handle<YieldTermStructure>(vars.termStructure),
            spreads, spreadDates);
        spreadedTermStructure->enableExtrapolation();

        Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
        Rate interpolatedZeroRate = spreadedTermStructure->zeroRate(t, vars.compounding);


        Real expectedRate = vars.termStructure->zeroRate(t, vars.compounding) +
            spread2->value();

        rDep[i].inputRate_ = rates[0];
        rDep[i].interpolatedZeroRate_ = interpolatedZeroRate;
    }
    cl::AdjointTestOutput out("AdjointPiecwZeroSprTermStruct//FlatInterpRight//output", { { "filename", "zeroRateOnRateByFlatInterpRight" }
                                                                                        , { "not_clear", "Not" }
                                                                                        , { "title", "Zero Rate dependence on Yield Rate By Flat Interpolation Right" }
                                                                                        , { "ylabel", "Zero Rate" } });
    out << rDep;
    logOut.log() << "Zero Rate dependence on Yield Rate By Flat Interpolation Right plot has been generated " << currentTime() << std::endl;
}

// Function outputLinearInterpolationMultipleSpreads(cl::AdjointTestOutput & logOut)
// Generating Zero Rate dependence on Yield Rate Dependence By Linear Interpolation Multiple Spreads plot.
// zeroRate is interpolated in a date point by the above mentioned method using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays) + spread quote.
// Parameter: logOut is reference to struct with log stream.
void outputLinearInterpolationMultipleSpreads(cl::AdjointTestOutput & logOut)
{
    Size steps = Size((MAX_RATE_RANGE - MIN_RATE_RANGE) / RATE_STEP) + 1;
    std::vector<Rate> rates = { MIN_RATE_RANGE, 0.033 };
    std::vector<RateDependence> rDep(steps);

    for (Size i = 0; i < steps; rates[0] += RATE_STEP, i++)
    {
        CommonVars vars(rates, getIntegerVector(rates.size()));
        std::vector<Handle<Quote> > spreads;
        boost::shared_ptr<SimpleQuote> spread1 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread2 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread3 = boost::make_shared<SimpleQuote>(0.035);
        boost::shared_ptr<SimpleQuote> spread4 = boost::make_shared<SimpleQuote>(0.04);
        spreads.push_back(Handle<Quote>(spread1));
        spreads.push_back(Handle<Quote>(spread2));
        spreads.push_back(Handle<Quote>(spread3));
        spreads.push_back(Handle<Quote>(spread4));

        std::vector<Date> spreadDates;
        spreadDates.push_back(vars.calendar.advance(vars.today, 90, Days));
        spreadDates.push_back(vars.calendar.advance(vars.today, 150, Days));
        spreadDates.push_back(vars.calendar.advance(vars.today, 30, Months));
        spreadDates.push_back(vars.calendar.advance(vars.today, 40, Months));

        Date interpolationDate = vars.today;

        boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
            boost::make_shared<PiecewiseZeroSpreadedTermStructure>(
            Handle<YieldTermStructure>(vars.termStructure),
            spreads, spreadDates);

        Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
        Rate interpolatedZeroRate = spreadedTermStructure->zeroRate(t, vars.compounding);


        Real expectedRate = vars.termStructure->zeroRate(t, vars.compounding) +
            spread1->value();

        rDep[i].inputRate_ = rates[0];
        rDep[i].interpolatedZeroRate_ = interpolatedZeroRate;
    }
    cl::AdjointTestOutput out("AdjointPiecwZeroSprTermStruct//LinInterpMultSpreads//output", { { "filename", "zeroRateOnRateByLinInterpMultSpreads" }
                                                                                             , { "not_clear", "Not" }
                                                                                             , { "title", "Zero Rate dependence on Yield Rate By Linear Interpolation Multiple Spreads" }
                                                                                             , { "ylabel", "Zero Rate" } });
    out << rDep;
    logOut.log() << "Zero Rate dependence on Yield Rate By Linear Interpolation Multiple Spreads plot has been generated " << currentTime() << std::endl;
}

// Function outputLinearInterpolation(cl::AdjointTestOutput & logOut)
// Generating Zero Rate dependence on Yield Rate Dependence By Linear Interpolation plot.
// zeroRate is interpolated in a date point by the above mentioned method using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays) + spread quote.
// Parameter: logOut is reference to struct with log stream.
void outputLinearInterpolation(cl::AdjointTestOutput & logOut)
{
    Size steps = Size((MAX_RATE_RANGE - MIN_RATE_RANGE) / RATE_STEP) + 1;
    std::vector<Rate> rates = { MIN_RATE_RANGE, 0.033 };
    std::vector<RateDependence> rDep(steps);

    for (Size i = 0; i < steps; rates[0] += RATE_STEP, i++)
    {
        CommonVars vars(rates, getIntegerVector(rates.size()));
        std::vector<Handle<Quote> > spreads;
        boost::shared_ptr<SimpleQuote> spread1 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread2 = boost::make_shared<SimpleQuote>(0.03);
        spreads.push_back(Handle<Quote>(spread1));
        spreads.push_back(Handle<Quote>(spread2));

        std::vector<Date> spreadDates;
        spreadDates.push_back(vars.calendar.advance(vars.today, 100, Days));
        spreadDates.push_back(vars.calendar.advance(vars.today, 150, Days));

        Date interpolationDate = vars.calendar.advance(vars.today, 3, Days);

        boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
            boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<Linear> >(
            Handle<YieldTermStructure>(vars.termStructure),
            spreads, spreadDates);

        Date d0 = vars.calendar.advance(vars.today, 100, Days);
        Date d1 = vars.calendar.advance(vars.today, 150, Days);
        Date d2 = vars.calendar.advance(vars.today, 120, Days);

        Real m = (0.03 - 0.02) / vars.dayCount.yearFraction(d0, d1);
        Real expectedRate = m * vars.dayCount.yearFraction(d0, d2) + 0.054;

        Time t = vars.dayCount.yearFraction(vars.settlementDate, interpolationDate);
        Rate interpolatedZeroRate = spreadedTermStructure->zeroRate(t, vars.compounding);

        rDep[i].inputRate_ = rates[0];
        rDep[i].interpolatedZeroRate_ = interpolatedZeroRate;
    }
    cl::AdjointTestOutput out("AdjointPiecwZeroSprTermStruct//LinInterp//output", { { "filename", "zeroRateOnRateByLinInterp" }
                                                                                  , { "not_clear", "Not" }
                                                                                  , { "title", "Zero Rate dependence on Yield Rate By Linear Interpolation" }
                                                                                  , { "ylabel", "Zero Rate" } });
    out << rDep;
    logOut.log() << "Zero Rate dependence on Yield Rate By Linear Interpolation plot has been generated " << currentTime() << std::endl;
}

// Function outputForwardFlatInterpolation(cl::AdjointTestOutput & logOut)
// Generating Zero Rate dependence on Yield Rate Dependence By Forward Flat Interpolation plot.
// zeroRate is interpolated in a date point by the above mentioned method using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays) + spread quote.
// Parameter: logOut is reference to struct with log stream.
void outputForwardFlatInterpolation(cl::AdjointTestOutput & logOut)
{
    Size steps = Size((MAX_RATE_RANGE - MIN_RATE_RANGE) / RATE_STEP) + 1;
    std::vector<Rate> rates = { MIN_RATE_RANGE, 0.033 };
    std::vector<RateDependence> rDep(steps);

    for (Size i = 0; i < steps; rates[0] += RATE_STEP, i++)
    {
        CommonVars vars(rates, getIntegerVector(rates.size()));
        std::vector<Handle<Quote> > spreads;
        boost::shared_ptr<SimpleQuote> spread1 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread2 = boost::make_shared<SimpleQuote>(0.03);
        spreads.push_back(Handle<Quote>(spread1));
        spreads.push_back(Handle<Quote>(spread2));

        std::vector<Date> spreadDates;
        spreadDates.push_back(vars.calendar.advance(vars.today, 75, Days));
        spreadDates.push_back(vars.calendar.advance(vars.today, 260, Days));

        Date interpolationDate = vars.today;

        boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
            boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<ForwardFlat> >(
            Handle<YieldTermStructure>(vars.termStructure),
            spreads, spreadDates);

        Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
        Rate interpolatedZeroRate = spreadedTermStructure->zeroRate(t, vars.compounding);


        Real expectedRate = vars.termStructure->zeroRate(t, vars.compounding) +
            spread1->value();

        rDep[i].inputRate_ = rates[0];
        rDep[i].interpolatedZeroRate_ = interpolatedZeroRate;
    }
    cl::AdjointTestOutput out("AdjointPiecwZeroSprTermStruct//ForwardFlatInterp//output", { { "filename", "zeroRateOnRateByForwardFlatInterp" }
                                                                                          , { "not_clear", "Not" }
                                                                                          , { "title", "Zero Rate dependence on Yield Rate By Forward Flat Interpolation" }
                                                                                          , { "ylabel", "Zero Rate" } });
    out << rDep;
    logOut.log() << "Zero Rate dependence on Yield Rate By Forward Flat Interpolation plot has been generated " << currentTime() << std::endl;
}

// Function outputBackwardFlatInterpolation(cl::AdjointTestOutput & logOut)
// Generating Zero Rate dependence on Yield Rate Dependence By Backward Flat Interpolation plot.
// zeroRate is interpolated in a date point by the above mentioned method using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays) + spread quote.
// Parameter: logOut is reference to struct with log stream.
void outputBackwardFlatInterpolation(cl::AdjointTestOutput & logOut)
{
    Size steps = Size((MAX_RATE_RANGE - MIN_RATE_RANGE) / RATE_STEP) + 1;
    std::vector<Rate> rates = { MIN_RATE_RANGE, 0.033 };
    std::vector<RateDependence> rDep(steps);

    for (Size i = 0; i < steps; rates[0] += RATE_STEP, i++)
    {
        CommonVars vars(rates, getIntegerVector(rates.size()));
        std::vector<Handle<Quote> > spreads;
        boost::shared_ptr<SimpleQuote> spread1 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread2 = boost::make_shared<SimpleQuote>(0.03);
        boost::shared_ptr<SimpleQuote> spread3 = boost::make_shared<SimpleQuote>(0.04);
        spreads.push_back(Handle<Quote>(spread1));
        spreads.push_back(Handle<Quote>(spread2));
        spreads.push_back(Handle<Quote>(spread3));

        std::vector<Date> spreadDates;
        spreadDates.push_back(vars.calendar.advance(vars.today, 100, Days));
        spreadDates.push_back(vars.calendar.advance(vars.today, 200, Days));
        spreadDates.push_back(vars.calendar.advance(vars.today, 300, Days));

        Date interpolationDate = vars.today;

        boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
            boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<BackwardFlat> >(
            Handle<YieldTermStructure>(vars.termStructure),
            spreads, spreadDates);

        Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
        Rate interpolatedZeroRate = spreadedTermStructure->zeroRate(t, vars.compounding);


        Real expectedRate = vars.termStructure->zeroRate(t, vars.compounding) +
            spread2->value();

        rDep[i].inputRate_ = rates[0];
        rDep[i].interpolatedZeroRate_ = interpolatedZeroRate;
    }
    cl::AdjointTestOutput out("AdjointPiecwZeroSprTermStruct//BackwardFlatInterp//output", { { "filename", "zeroRateOnRateByBackwardFlatInterp" }
                                                                                           , { "not_clear", "Not" }
                                                                                           , { "title", "Zero Rate dependence on Yield Rate By Backward Flat Interpolation" }
                                                                                           , { "ylabel", "Zero Rate" } });
    out << rDep;
    logOut.log() << "Zero Rate dependence on Yield Rate By Backward Flat Interpolation plot has been generated " << currentTime() << std::endl;
}

// Function outputForwardFlatInterpolation(cl::AdjointTestOutput & logOut)
// Generating Zero Rate dependence on Yield Rate Dependence By Default (Linear) Interpolation plot.
// zeroRate is interpolated in a date point by the above mentioned method using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays) + spread quote.
// Parameter: logOut is reference to struct with log stream.
void outputDefaultInterpolation(cl::AdjointTestOutput & logOut)
{
    Size steps = Size((MAX_RATE_RANGE - MIN_RATE_RANGE) / RATE_STEP) + 1;
    std::vector<Rate> rates = { MIN_RATE_RANGE, 0.033 };
    std::vector<RateDependence> rDep(steps);

    for (Size i = 0; i < steps; rates[0] += RATE_STEP, i++)
    {
        CommonVars vars(rates, getIntegerVector(rates.size()));
        std::vector<Handle<Quote> > spreads;
        boost::shared_ptr<SimpleQuote> spread1 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread2 = boost::make_shared<SimpleQuote>(0.02);
        spreads.push_back(Handle<Quote>(spread1));
        spreads.push_back(Handle<Quote>(spread2));

        std::vector<Date> spreadDates;
        spreadDates.push_back(vars.calendar.advance(vars.today, 75, Days));
        spreadDates.push_back(vars.calendar.advance(vars.today, 160, Days));

        Date interpolationDate = vars.today;

        boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
            boost::make_shared<PiecewiseZeroSpreadedTermStructure>(
            Handle<YieldTermStructure>(vars.termStructure),
            spreads, spreadDates);

        Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
        Rate interpolatedZeroRate = spreadedTermStructure->zeroRate(t, vars.compounding);


        Real expectedRate = vars.termStructure->zeroRate(t, vars.compounding) +
            spread1->value();

        rDep[i].inputRate_ = rates[0];
        rDep[i].interpolatedZeroRate_ = interpolatedZeroRate;
    }
    cl::AdjointTestOutput out("AdjointPiecwZeroSprTermStruct//DefaultInterp//output", { { "filename", "zeroRateOnRateByDefaultInterp" }
                                                                                      , { "not_clear", "Not" }
                                                                                      , { "title", "Zero Rate dependence on Yield Rate By Default Interpolation" }
                                                                                      , { "ylabel", "Zero Rate" } });
    out << rDep;
    logOut.log() << "Zero Rate dependence on Yield Rate By Default Interpolation plot has been generated " << currentTime() << std::endl;
}

// Function outputForwardFlatInterpolation(cl::AdjointTestOutput & logOut)
// Generating Zero Rate dependence on Yield Rate Dependence using Interpolation Factory (Cubic Interpolation Spline) plot.
// zeroRate is interpolated in a date point by the above mentioned method using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays) + spread quote.
// Parameter: logOut is reference to struct with log stream.
void outputSetInterpolationFactory(cl::AdjointTestOutput & logOut)
{
    Size steps = Size((MAX_RATE_RANGE - MIN_RATE_RANGE) / RATE_STEP) + 1;
    std::vector<Rate> rates = { MIN_RATE_RANGE, 0.33 };
    std::vector<RateDependence> rDep(steps);

    for (Size i = 0; i < steps; rates[0] += RATE_STEP, i++)
    {
        CommonVars vars(rates, getIntegerVector(rates.size()));
        std::vector<Handle<Quote> > spreads;
        boost::shared_ptr<SimpleQuote> spread1 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread2 = boost::make_shared<SimpleQuote>(0.03);
        boost::shared_ptr<SimpleQuote> spread3 = boost::make_shared<SimpleQuote>(0.01);
        spreads.push_back(Handle<Quote>(spread1));
        spreads.push_back(Handle<Quote>(spread2));
        spreads.push_back(Handle<Quote>(spread3));

        std::vector<Date> spreadDates;
        spreadDates.push_back(vars.calendar.advance(vars.today, 8, Months));
        spreadDates.push_back(vars.calendar.advance(vars.today, 15, Months));
        spreadDates.push_back(vars.calendar.advance(vars.today, 25, Months));

        Date interpolationDate = vars.today;

        boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure;

        Frequency freq = NoFrequency;

        Cubic factory;
        factory = Cubic(CubicInterpolation::Spline, false);

        spreadedTermStructure =
            boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<Cubic> >(
            Handle<YieldTermStructure>(vars.termStructure),
            spreads, spreadDates, vars.compounding,
            freq, vars.dayCount, factory);

        Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
        Rate interpolatedZeroRate = spreadedTermStructure->zeroRate(t, vars.compounding);


        Real expectedRate = vars.termStructure->zeroRate(t, vars.compounding) +
            Real(0.026065770863);

        rDep[i].inputRate_ = rates[0];
        rDep[i].interpolatedZeroRate_ = interpolatedZeroRate;
    }
    cl::AdjointTestOutput out("AdjointPiecwZeroSprTermStruct//SetInterpFactory//output", { { "filename", "zeroRateOnRateByCubicInterpolationSpline" }
                                                                                         , { "not_clear", "Not" }
                                                                                         , { "title", "Zero Rate dependence on Yield Rate By Cubic Interpolation Spline" }
                                                                                         , { "ylabel", "Zero Rate" } });
    out << rDep;
    logOut.log() << "Zero Rate dependence on Yield Rate By Cubic Interpolation Spline plot has been generated " << currentTime() << std::endl;
}

// Function outputForwardFlatInterpolation(cl::AdjointTestOutput & logOut)
// Generating Zero Rate dependence on Yield Rate Dependence By Backward Flat Interpolation (with Quote Changing) plot.
// zeroRate is interpolated in a date point by the above mentioned method using the following formula:
// zeroRate = Interpolation(interpolationDate + settlementDays) + spread quote.
// Parameter: logOut is reference to struct with log stream.
void outputQuoteChanging(cl::AdjointTestOutput & logOut)
{
    Size steps = Size((MAX_RATE_RANGE - MIN_RATE_RANGE) / RATE_STEP) + 1;
    std::vector<Rate> rates = { MIN_RATE_RANGE, 0.33 };
    std::vector<RateDependence> rDep(steps);

    for (Size i = 0; i < steps; rates[0] += RATE_STEP, i++)
    {
        CommonVars vars(rates, getIntegerVector(rates.size()));
        std::vector<Handle<Quote> > spreads;
        boost::shared_ptr<SimpleQuote> spread1 = boost::make_shared<SimpleQuote>(0.02);
        boost::shared_ptr<SimpleQuote> spread2 = boost::make_shared<SimpleQuote>(0.03);
        spreads.push_back(Handle<Quote>(spread1));
        spreads.push_back(Handle<Quote>(spread2));

        std::vector<Date> spreadDates;
        spreadDates.push_back(vars.calendar.advance(vars.today, 75, Days));
        spreadDates.push_back(vars.calendar.advance(vars.today, 260, Days));

        Date interpolationDate = vars.today;

        boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
            boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<ForwardFlat> >(
            Handle<YieldTermStructure>(vars.termStructure),
            spreads, spreadDates);

        Time t = vars.dayCount.yearFraction(vars.today, interpolationDate);
        Rate interpolatedZeroRate = spreadedTermStructure->zeroRate(t, vars.compounding);


        Real expectedRate = vars.termStructure->zeroRate(t, vars.compounding) +
            spread1->value();

        rDep[i].inputRate_ = rates[0];
        rDep[i].interpolatedZeroRate_ = interpolatedZeroRate;
    }
    cl::AdjointTestOutput out("AdjointPiecwZeroSprTermStruct//QuoteChanging//output", { { "filename", "zeroRateOnRateWithQuoteChanging" }
                                                                                      , { "not_clear", "Not" }
                                                                                      , { "title", "Zero Rate dependence on Yield Rate with Quote Changing" }
                                                                                      , { "ylabel", "Zero Rate" } });
    out << rDep;
    logOut.log() << "Zero Rate dependence on Yield Rate With Quote Changing plot has been generated " << currentTime() << std::endl;
}

test_suite* AdjointPiecewiseZeroSpreadedTermStructureTest::suite() {
    test_suite* suite = BOOST_TEST_SUITE("Interpolated piecewise zero spreaded yield curve tests");
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointPiecewiseZeroSpreadedTermStructureTest::testFlatInterpolationLeft));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointPiecewiseZeroSpreadedTermStructureTest::testFlatInterpolationRight));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointPiecewiseZeroSpreadedTermStructureTest::testLinearInterpolationMultipleSpreads));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointPiecewiseZeroSpreadedTermStructureTest::testLinearInterpolation));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointPiecewiseZeroSpreadedTermStructureTest::testBackwardFlatInterpolation));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointPiecewiseZeroSpreadedTermStructureTest::testForwardFlatInterpolation));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointPiecewiseZeroSpreadedTermStructureTest::testDefaultInterpolation));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointPiecewiseZeroSpreadedTermStructureTest::testSetInterpolationFactory));
    suite->add(QUANTLIB_TEST_CASE(
        &AdjointPiecewiseZeroSpreadedTermStructureTest::testQuoteChanging));

    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_piecewise_zero_spreaded_term_structure)

BOOST_AUTO_TEST_CASE(testFlatInterpolationLeft)
{
    BOOST_CHECK(AdjointPiecewiseZeroSpreadedTermStructureTest::testFlatInterpolationLeft());
}

BOOST_AUTO_TEST_CASE(testFlatInterpolationRight)
{
    BOOST_CHECK(AdjointPiecewiseZeroSpreadedTermStructureTest::testFlatInterpolationRight());
}

BOOST_AUTO_TEST_CASE(testLinearInterpolationMultipleSpreads)
{
    BOOST_CHECK(AdjointPiecewiseZeroSpreadedTermStructureTest::testLinearInterpolationMultipleSpreads());
}

BOOST_AUTO_TEST_CASE(testLinearInterpolation)
{
    BOOST_CHECK(AdjointPiecewiseZeroSpreadedTermStructureTest::testLinearInterpolation());
}

BOOST_AUTO_TEST_CASE(testForwardFlatInterpolation)
{
    BOOST_CHECK(AdjointPiecewiseZeroSpreadedTermStructureTest::testForwardFlatInterpolation());
}

BOOST_AUTO_TEST_CASE(testBackwardFlatInterpolation)
{
    BOOST_CHECK(AdjointPiecewiseZeroSpreadedTermStructureTest::testBackwardFlatInterpolation());
}

BOOST_AUTO_TEST_CASE(testDefaultLinearInterpolation)
{
    BOOST_CHECK(AdjointPiecewiseZeroSpreadedTermStructureTest::testDefaultInterpolation());
}

BOOST_AUTO_TEST_CASE(testSetInterpolationFactory)
{
    BOOST_CHECK(AdjointPiecewiseZeroSpreadedTermStructureTest::testSetInterpolationFactory());
}

BOOST_AUTO_TEST_CASE(testBackwardFlatInterpQuoteDependence)
{
    BOOST_CHECK(AdjointPiecewiseZeroSpreadedTermStructureTest::testQuoteChanging());
}

BOOST_AUTO_TEST_SUITE_END()

#endif


