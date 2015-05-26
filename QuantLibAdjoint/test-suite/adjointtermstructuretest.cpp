/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2003 RiskMap srl
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

// Based on termstructure.cpp from Quantlib/test-suite.

#include "adjointtermstructuretest.hpp"
#include "utilities.hpp"
#include <ql/termstructures/yield/ratehelpers.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/yield/piecewiseyieldcurve.hpp>
#include <ql/termstructures/yield/impliedtermstructure.hpp>
#include <ql/termstructures/yield/forwardspreadedtermstructure.hpp>
#include <ql/termstructures/yield/zerospreadedtermstructure.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/calendars/nullcalendar.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/daycounters/thirty360.hpp>
#include <ql/math/comparison.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/currency.hpp>
#include <ql/utilities/dataformatters.hpp>
#include <boost/timer.hpp>
#include <iostream>
#include "adjointtestutilities.hpp"

using namespace QuantLib;
using namespace boost::unit_test_framework;

void outputImplied(cl::AdjointTestOutput & out);
void outputFSpreaded(cl::AdjointTestOutput & out);
void outputZSpreaded(cl::AdjointTestOutput & out);

namespace {

#if defined CL_GRAPH_GEN
    Size MIN_RATE_SIZE = 30;
    Size MAX_RATE_SIZE = 110;

    Real MIN_RATE_RANGE = 0.0;
#else
    Size MIN_RATE_SIZE = 30;
    Size MAX_RATE_SIZE = 30;
    Real MIN_RATE_RANGE = 1.0;
#endif

    

    Real MAX_RATE_RANGE = 80.0;
    Real RATE_STEP = 1;

    std::vector<double> getDoubleVector(Size size)
    {
        std::vector<double> vec =
        {
              0.50, 0.69, 0.57, 0.27, 0.48, 0.85, 0.07, 0.74, 0.47, 0.78, 0.07, 0.68
            , 0.92, 0.10, 0.93, 0.33, 0.20, 0.40, 0.04, 0.73, 0.91, 0.19, 0.50, 0.87
            , 0.81, 0.34, 0.10, 0.52, 0.43, 0.54, 0.34, 0.31, 0.75, 0.51, 0.30, 0.69
            , 0.36, 0.16, 0.90, 0.32, 0.40, 0.20, 0.64, 0.64, 0.72, 0.90, 0.78, 0.74
            , 0.27, 0.39, 0.24, 0.60, 0.88, 0.20, 0.97, 0.33, 0.59, 0.52, 0.91, 0.77
            , 0.73, 0.63, 0.89, 0.88, 0.54, 0.72, 0.92, 0.32, 0.09, 0.56, 0.68, 0.91
            , 0.82, 0.72, 0.35, 0.41, 0.42, 0.23, 0.92, 0.25, 0.67, 0.35, 0.58, 0.52
            , 0.65, 0.74, 0.19, 0.75, 0.06, 0.08, 0.51, 0.97, 0.01, 0.99, 0.32, 0.40
            , 0.69, 0.71, 0.83, 0.64, 0.08, 0.61, 0.56, 0.27, 0.00, 0.69, 0.70, 0.89
            , 0.98, 0.15, 0.06, 0.96, 0.18, 0.07, 0.05, 0.43, 0.58, 0.92, 0.00, 0.75
            , 0.38, 0.17, 0.31, 0.74, 0.99, 0.05, 0.59, 0.91, 0.94, 0.98, 0.44, 0.20
            , 0.36, 0.99, 0.83, 0.48, 0.28, 0.76, 0.41, 0.74, 0.88, 0.82, 0.08, 0.21
            , 0.91, 0.11, 0.26, 0.46, 0.59, 0.26, 0.74, 0.55, 0.39, 0.19, 0.38, 0.18
            , 0.59, 0.30, 0.10, 0.77, 0.46, 0.14, 0.90, 0.34, 0.08, 0.40, 0.20, 0.62
            , 0.32, 0.76, 0.50, 0.86, 0.48, 0.75, 0.59, 0.06, 0.36, 0.25, 0.39, 0.21
            , 0.79, 0.06, 0.55, 0.41, 0.60, 0.79, 0.36, 0.91, 0.93, 0.34, 0.06, 0.70
            , 0.41, 0.22, 0.30, 0.49, 0.59, 0.35, 0.51, 0.77, 0.53, 0.46, 0.63, 0.23
            , 0.98, 0.91, 0.65, 0.18, 0.71, 0.52, 0.98, 0.77, 0.80, 0.97, 0.31, 0.88
            , 0.26, 0.53, 0.85, 0.67
        };
        return std::vector<double>(vec.begin(), vec.begin() + size);
    }

    struct RateDiscount
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "Swap Rate", "Base Discount", "Discount", "Implied Discount"
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, RateDiscount& v)
        {
                stm << v.inputRate_
                    << ";" << v.baseDiscount_
                    << ";" << v.discount_
                    << ";" << v.impliedDiscount_ << std::endl;

                return stm;
            }

        Real inputRate_;
        Real baseDiscount_;
        Real discount_;
        Real impliedDiscount_;
    };

    struct ZeroRateConsistency
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
            operator << (stream_type& stm, ZeroRateConsistency& v)
        {
                stm << v.inputRate_
                    << ";" << v.functionValue_ << std::endl;

                return stm;
            }

        Real inputRate_;
        Real functionValue_;
    };

    struct ForwardRateConsistency
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
            operator << (stream_type& stm, ForwardRateConsistency& v)
        {
                stm << v.inputRate_
                    << ";" << v.functionValue_ << std::endl;

                return stm;
            }

        Real inputRate_;
        Real functionValue_;
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
        boost::shared_ptr<YieldTermStructure> termStructure;
        boost::shared_ptr<YieldTermStructure> dummyTermStructure;

        // cleanup
        SavedSettings backup;

        // setup
        CommonVars(std::vector<Datum> swapData)
        {
            std::vector<Datum> depositData = 
                {
                    { 1, Months, 4.581 },
                    { 2, Months, 4.573 },
                    { 3, Months, 4.557 },
                    { 6, Months, 4.496 },
                    { 9, Months, 4.490 }
                };
            calendar = TARGET();
            settlementDays = 2;
            Date today = calendar.adjust(Date::todaysDate());
            Settings::instance().evaluationDate() = today;
            Date settlement = calendar.advance(today, settlementDays, Days);
            Size deposits = depositData.size(),
                swaps = swapData.size();

            std::vector<boost::shared_ptr<RateHelper> > instruments(
                deposits + swaps);
            for (Size i = 0; i<deposits; i++) {
                instruments[i] = boost::shared_ptr<RateHelper>(new
                    DepositRateHelper(depositData[i].rate / 100,
                    depositData[i].n*depositData[i].units,
                    settlementDays, calendar,
                    ModifiedFollowing, true,
                    Actual360()));
            }
            boost::shared_ptr<IborIndex> index(new IborIndex("dummy",
                6 * Months,
                settlementDays,
                Currency(),
                calendar,
                ModifiedFollowing,
                false,
                Actual360()));
            for (Size i = 0; i<swaps; ++i) {
                instruments[i + deposits] = boost::shared_ptr<RateHelper>(new
                    SwapRateHelper(swapData[i].rate / 100,
                    swapData[i].n*swapData[i].units,
                    calendar,
                    Annual, Unadjusted, Thirty360(),
                    index));
            }
            termStructure = boost::shared_ptr<YieldTermStructure>(new
                PiecewiseYieldCurve<Discount, LogLinear>(settlement,
                instruments, Actual360()));
            dummyTermStructure = boost::shared_ptr<YieldTermStructure>(new
                PiecewiseYieldCurve<Discount, LogLinear>(settlement,
                instruments, Actual360()));
        }
    };

}

// Method testZSpreaded()
// Testing consistency of zero-spreaded term structure.
// Zero rate is interpolated in a date point by the log-linear method using the following formula:
// zero = Interpolation(testDate).
// The zero rate is differentiated
// with respect to input swap rate vector.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same plot.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it throws an exception and returns false.
bool AdjointTermStructureTest::testZSpreaded()
{
    BOOST_TEST_MESSAGE("Testing consistency of zero-spreaded term structure...");

    std::vector<PerformanceTime> pTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<AdjointTime> adjTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<TapeSize> tapeSize(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    cl::AdjointTestOutput outPerform("AdjointTermStructure//ZSpreaded", { { "filename", "AdjointPerformance" }
                                                                     , { "not_clear", "Not" }
                                                                     , { "title", "Zero rate differentiation performance with respect to swap rate" }
                                                                     , { "xlabel", "Number of swap rates" }
                                                                     , { "ylabel", "Time (s)" }
                                                                     , { "line_box_width", "-5" } });
    cl::AdjointTestOutput outAdjoint("AdjointTermStructure//ZSpreaded", { { "filename", "Adjoint" }
                                                                    , { "not_clear", "Not" }
                                                                    , { "title", " Zero rate adjoint differentiation with respect to swap rate" }
                                                                    , { "xlabel", "Number of swap rates" }
                                                                    , { "ylabel", "Adjoint calculations time" } });
    cl::AdjointTestOutput outSize("AdjointTermStructure//ZSpreaded", { { "filename", "TapeSize" }
                                                                     , { "not_clear", "Not" }
                                                                     , { "title", "Tape size dependence on number of swap rates" }
                                                                     , { "xlabel", "Number of swap rates" }
                                                                     , { "ylabel", "Memory(MB)" } });
    bool ok = true;

    for (Size rate_size = MIN_RATE_SIZE; rate_size <= MAX_RATE_SIZE; rate_size ++)
    {
        pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        tapeSize[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        std::vector<double> rate_values_double = getDoubleVector(rate_size);
        std::vector<cl::TapeDouble> rate_values(rate_values_double.size());
        for (Size i = 0; i < rate_values_double.size(); i++)
        {
            rate_values[i] = Real(rate_values_double[i]);
        }
        outPerform.log() << std::endl << "Start of tape recording " << currentTime() << std::endl;
        boost::timer timer;

        Independent(rate_values);

        std::vector<Datum> swapData(rate_size);
        
        for (Size i = 0; i < rate_size; i++)
        {
            swapData[i] = { i * 3 + 10, Months, rate_values[i] };
        }
        
        CommonVars vars(swapData);

        Real tolerance = 1.0e-10;
        boost::shared_ptr<Quote> me(new SimpleQuote(0.01));
        Handle<Quote> mh(me);
        boost::shared_ptr<YieldTermStructure> spreaded(
            new ZeroSpreadedTermStructure(
            Handle<YieldTermStructure>(vars.termStructure), mh));
        Date testDate = vars.termStructure->referenceDate() + 5 * Years;
        DayCounter rfdc = vars.termStructure->dayCounter();
        Rate zero = vars.termStructure->zeroRate(testDate, rfdc,
            Continuous, NoFrequency);
        Rate spreadedZero = spreaded->zeroRate(testDate, rfdc,
            Continuous, NoFrequency);
        pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ = timer.elapsed();
        outPerform.log() << "End of tape recording " << currentTime() << std::endl;
        std::vector<cl::TapeDouble> zero_value(1, zero);
        cl::TapeFunction<double> f(rate_values, zero_value);
        tapeSize[rate_size - MIN_RATE_SIZE].memory_ = f.Memory();
        outPerform.log() << "Tape Memory in bytes: " << f.Memory() << std::endl;
        outPerform.log() << "Time for tape recording: " << pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ << std::endl;
        std::vector<double> dy_forward(rate_values.size());
        std::vector<double> dy_reverse(rate_values.size());

        pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = std::min(
            // Start differentiation in Reverse mode
            gradReverse(f, dy_reverse, outPerform, false, false),
            //Start differentiation in Forward mode
            gradForward(f, dy_forward, outPerform, false, false));

        double h = 1.0e-6;
        std::vector<Real> dy_diff(rate_values.size());
        timer.restart();
        //Finite differences
        for (Size i = 0; i < rate_values.size(); i++)
        {
            swapData[i].rate += h;

            CommonVars vars(swapData);

            Real tolerance = 1.0e-10;
            boost::shared_ptr<Quote> me(new SimpleQuote(0.01));
            Handle<Quote> mh(me);
            boost::shared_ptr<YieldTermStructure> spreaded(
                new ZeroSpreadedTermStructure(
                Handle<YieldTermStructure>(vars.termStructure), mh));
            Date testDate = vars.termStructure->referenceDate() + 5 * Years;
            DayCounter rfdc = vars.termStructure->dayCounter();
            Rate zero_h = vars.termStructure->zeroRate(testDate, rfdc,
                Continuous, NoFrequency);

            dy_diff[i] = (zero_h - zero) / h;
            swapData[i].rate -= h;
        }
        pTime[rate_size - MIN_RATE_SIZE].timeAnalytical_ = timer.elapsed();
        adjTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_;
        adjTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_;
        ok &= checkWithFiniteDiff(dy_forward, dy_reverse, dy_diff, outPerform, 1e-2, 1e-5);
    }
    outPerform << pTime;
    outPerform.log() << "Zero-spreaded Performance Time plot has been generated " << currentTime() << std::endl;
    outAdjoint << adjTime;
    outPerform.log() << "Zero-spreaded Adjoint Performance plot has been generated " << currentTime() << std::endl;
    outSize << tapeSize;
    outPerform.log() << "Tape size dependence on number of independent variables plot has been generated " << currentTime() << std::endl;
    outputZSpreaded(outPerform);
    return ok;
}

// Method testFSpreaded()
// Testing consistency of forward-spreaded term structure.
// Forward rate is interpolated in a date point by the log-linear method using the following formula:
// forward = Interpolation(testDate).
// The forward rate is differentiated
// with respect to input depostit rate and swap rate vectors.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same plot.
// Method returns true if adjoint derivatives are close enough to finite difference derivative; otherwise, it throws an exception and returns false.
bool AdjointTermStructureTest::testFSpreaded() {

    BOOST_TEST_MESSAGE("Testing consistency of forward-spreaded term structure...");

    std::vector<PerformanceTime> pTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<AdjointTime> adjTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<TapeSize> tapeSize(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    cl::AdjointTestOutput outPerform("AdjointTermStructure//FSpreaded", { { "filename", "AdjointPerformance" }
                                                                     , { "not_clear", "Not" }
                                                                     , { "title", "Foward rate differentiation performance with respect to swap rate" }
                                                                     , { "xlabel", "Number of swap rates" }
                                                                     , { "ylabel", "Time (s)" }
                                                                     , { "line_box_width", "-5" } });
    cl::AdjointTestOutput outAdjoint("AdjointTermStructure//FSpreaded", { { "filename", "Adjoint" }
                                                                    , { "not_clear", "Not" }
                                                                    , { "title", "Forward rate adjoint differentiation with respect to swap rate" }
                                                                    , { "xlabel", "Number of swap rates" }
                                                                    , { "ylabel", "Adjoint calculations time" } });
    cl::AdjointTestOutput outSize("AdjointTermStructure//FSpreaded", { { "filename", "TapeSize" }
                                                                     , { "not_clear", "Not" }
                                                                     , { "title", "Tape size dependence on number of swap rates" }
                                                                     , { "xlabel", "Number of swap rates" }
                                                                     , { "ylabel", "Memory(MB)" } });
    bool ok = true;

    for (Size rate_size = MIN_RATE_SIZE; rate_size <= MAX_RATE_SIZE; rate_size ++)
    {
        pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        tapeSize[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        std::vector<double> rate_values_double = getDoubleVector(rate_size);
        std::vector<cl::TapeDouble> rate_values(rate_values_double.size());
        for (Size i = 0; i < rate_values_double.size(); i++)
        {
            rate_values[i] = Real(rate_values_double[i]);
        }
        outPerform.log() << std::endl << "Start of tape recording " << currentTime() << std::endl;
        boost::timer timer;
        Independent(rate_values);

        std::vector<Datum> swapData(rate_size);

        for (Size i = 0; i < rate_size; i++)
        {
            swapData[i] = { i * 3 + 10, Months, rate_values[i] };
        }

        CommonVars vars(swapData);

        Real tolerance = 1.0e-10;
        boost::shared_ptr<Quote> me(new SimpleQuote(0.01));
        Handle<Quote> mh(me);
        boost::shared_ptr<YieldTermStructure> spreaded(
            new ForwardSpreadedTermStructure(
            Handle<YieldTermStructure>(vars.termStructure), mh));
        Date testDate = vars.termStructure->referenceDate() + 5 * Years;
        DayCounter tsdc = vars.termStructure->dayCounter();
        DayCounter sprdc = spreaded->dayCounter();
        Rate forward = vars.termStructure->forwardRate(testDate, testDate, tsdc,
            Continuous, NoFrequency);
        Rate spreadedForward = spreaded->forwardRate(testDate, testDate, sprdc,
            Continuous, NoFrequency);
        pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ = timer.elapsed();
        outPerform.log() << "End of tape recording " << currentTime() << std::endl;
        std::vector<cl::TapeDouble> forward_value(1, forward);
        cl::TapeFunction<double> f(rate_values, forward_value);
        tapeSize[rate_size - MIN_RATE_SIZE].memory_ = f.Memory();
        outPerform.log() << "Tape Memory in bytes: " << f.Memory() << std::endl;
        outPerform.log() << "Time for tape recording: " << pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ << std::endl;
        std::vector<double> dy_forward(rate_values.size());
        std::vector<double> dy_reverse(rate_values.size());

        pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = std::min(
            // Start differentiation in Reverse mode
            gradReverse(f, dy_reverse, outPerform, false, false),
            //Start differentiation in Forward mode
            gradForward(f, dy_forward, outPerform, false, false));

        double h = 1.0e-6;
        std::vector<Real> dy_diff(rate_values.size());
        timer.restart();
        //Finite differences
        for (Size i = 0; i < rate_values.size(); i++)
        {
            swapData[i].rate += h;

            CommonVars vars(swapData);

            Real tolerance = 1.0e-10;
            boost::shared_ptr<Quote> me(new SimpleQuote(0.01));
            Handle<Quote> mh(me);
            boost::shared_ptr<YieldTermStructure> spreaded(
                new ForwardSpreadedTermStructure(
                Handle<YieldTermStructure>(vars.termStructure), mh));
            Date testDate = vars.termStructure->referenceDate() + 5 * Years;
            DayCounter tsdc = vars.termStructure->dayCounter();
            DayCounter sprdc = spreaded->dayCounter();
            Rate forward_h = vars.termStructure->forwardRate(testDate, testDate, tsdc,
                Continuous, NoFrequency);

            dy_diff[i] = (forward_h - forward) / h;

            swapData[i].rate -= h;
        }
        pTime[rate_size - MIN_RATE_SIZE].timeAnalytical_ = timer.elapsed();
        adjTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_;
        adjTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_;
        ok &= checkWithFiniteDiff(dy_forward, dy_reverse, dy_diff, outPerform, 1e-2, 1e-5);
    }
    outPerform << pTime;
    outPerform.log() << "Forward-spreaded Performance Time plot has been generated " << currentTime() << std::endl;
    outAdjoint << adjTime;
    outPerform.log() << "Forward-spreaded Adjoint Adjoint Performance plot has been generated " << currentTime() << std::endl;
    outSize << tapeSize;
    outPerform.log() << "Tape size dependence on number of independent variables plot has been generated " << currentTime() << std::endl;
    outputFSpreaded(outPerform);
    return ok;
}

// Method testImplied()
// Testing consistency of yield term structure.
// Discounts are interpolated in a date point by the log-linear method using the following formula:
// discount = exp( - InterpolatedInterestRate(testDate, newSettlementDate) * time).
// The output quantities (base discount, discount and implied discount) are differentiated
// with respect to input depostit rate and swap rate vectors.
// Dependencies of time for differentiation using adjoint mode and
// without it are plotted on the same plot.
// Method returns true if adjoint derivatives are close enough to finite difference derivatives; otherwise, it throws an exception and returns false.
bool AdjointTermStructureTest::testImplied() {
    
    BOOST_TEST_MESSAGE("Testing consistency of implied term structure...");
    std::vector<PerformanceTime> pTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<AdjointTime> adjTime(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    std::vector<TapeSize> tapeSize(MAX_RATE_SIZE - MIN_RATE_SIZE + 1);
    cl::AdjointTestOutput outPerform("AdjointTermStructure//Implied", { { "filename", "AdjointPerformance" }
                                                                   , { "not_clear", "Not" }
                                                                   , { "title", "Base discount, discount, implied discount differentiation performance with respect to swap rate" }
                                                                   , { "xlabel", "Number of swap rates" }
                                                                   , { "ylabel", "Time (s)" }
                                                                   , { "line_box_width", "-5" } });
    cl::AdjointTestOutput outAdjoint("AdjointTermStructure//Implied", { { "filename", "Adjoint" }
                                                                  , { "not_clear", "Not" }
                                                                  , { "title", "Implied Adjoint Performance" }
                                                                  , { "xlabel", "Number of swap rates" }
                                                                  , { "ylabel", "Base discount, discount, implied discount adjoint differentiation with respect to swap rate" } });
    cl::AdjointTestOutput outSize("AdjointTermStructure//Implied", { { "filename", "TapeSize" }
                                                                   , { "not_clear", "Not" }
                                                                   , { "title", "Tape size dependence on number of swap rates" }
                                                                   , { "xlabel", "Number of swap rates" }
                                                                   , { "ylabel", "Memory(MB)" } });
    bool ok = true;
    for (Size rate_size = MIN_RATE_SIZE; rate_size <= MAX_RATE_SIZE; rate_size ++)
    {
        pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        tapeSize[rate_size - MIN_RATE_SIZE].indepVarNumber_ = rate_size;
        std::vector<double> rate_values_double = getDoubleVector(rate_size);
        std::vector<cl::TapeDouble> rate_values(rate_values_double.size());
        for (Size i = 0; i < rate_values_double.size(); i++)
        {
            rate_values[i] = Real(rate_values_double[i]);
        }
        outPerform.log() << std::endl << "Start of tape recording " << currentTime() << std::endl;
        boost::timer timer;
        Independent(rate_values);
        std::vector<Datum> swapData(rate_size);

        for (Size i = 0; i < rate_size; i++)
        {
            swapData[i] = { i * 3 + 10, Months, rate_values[i] };
        }

        CommonVars vars(swapData);

        Real tolerance = 1.0e-10;
        Date today = Settings::instance().evaluationDate();
        Date newToday = today + 3 * Years;
        Date newSettlement = vars.calendar.advance(newToday,
            vars.settlementDays, Days);
        Date testDate = newSettlement + 5 * Years;
        boost::shared_ptr<YieldTermStructure> implied(
            new ImpliedTermStructure(Handle<YieldTermStructure>(vars.termStructure),
            newSettlement));
        DiscountFactor baseDiscount = vars.termStructure->discount(newSettlement);
        DiscountFactor discount = vars.termStructure->discount(testDate);
        DiscountFactor impliedDiscount = implied->discount(testDate);
        pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ = timer.elapsed();
        outPerform.log() << "End of tape recording " << currentTime() << std::endl;
        std::vector<cl::TapeDouble> discount_values = { baseDiscount, discount, impliedDiscount };

        cl::TapeFunction<double> f(rate_values, discount_values);
        tapeSize[rate_size - MIN_RATE_SIZE].memory_ = f.Memory();
        outPerform.log() << "Tape Memory in bytes: " << f.Memory() << std::endl;
        outPerform.log() << "Time for tape recording: " << pTime[rate_size - MIN_RATE_SIZE].timeTapeRecording_ << std::endl;
        timer.restart();
        //Compute derivatives using Jacobian
        pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = timer.elapsed();
        timer.restart();
        std::vector<double> dy_forward(discount_values.size() * rate_values.size());
        std::vector<double> dy_reverse(discount_values.size() * rate_values.size());
        pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = std::min(
            // Start differentiation in Reverse mode
            gradReverse(f, dy_reverse, outPerform, false, false),
            //Start differentiation in Forward mode
            gradForward(f, dy_forward, outPerform, false, false));
        double h = 1.0e-6;
        std::vector<Real> dy_diff(discount_values.size() * rate_values.size());
        timer.restart();
        //Finite differences
        for (Size i = 0; i < rate_values.size(); i++)
        {
            swapData[i].rate += h;

            CommonVars vars(swapData);

            Real tolerance = 1.0e-10;
            Date today = Settings::instance().evaluationDate();
            Date newToday = today + 3 * Years;
            Date newSettlement = vars.calendar.advance(newToday,
                vars.settlementDays, Days);
            Date testDate = newSettlement + 5 * Years;
            boost::shared_ptr<YieldTermStructure> implied(
                new ImpliedTermStructure(Handle<YieldTermStructure>(vars.termStructure),
                newSettlement));
            DiscountFactor baseDiscount_h = vars.termStructure->discount(newSettlement);
            DiscountFactor discount_h = vars.termStructure->discount(testDate);
            DiscountFactor impliedDiscount_h = implied->discount(testDate);
            dy_diff[i] = (baseDiscount_h - baseDiscount) / h;
            dy_diff[rate_values.size() + i] = (discount_h - discount) / h;
            dy_diff[2 * rate_values.size() + i] = (impliedDiscount_h - impliedDiscount) / h;

            swapData[i].rate -= h;
        }
        pTime[rate_size - MIN_RATE_SIZE].timeAnalytical_ = timer.elapsed();
        adjTime[rate_size - MIN_RATE_SIZE].indepVarNumber_ = pTime[rate_size - MIN_RATE_SIZE].indepVarNumber_;
        adjTime[rate_size - MIN_RATE_SIZE].timeAdjoint_ = pTime[rate_size - MIN_RATE_SIZE].timeAdjoint_;
        ok &= checkWithFiniteDiff(dy_forward, dy_reverse, dy_diff, outPerform, 1e-2, 1e-5);
    }
    outPerform << pTime;
    outPerform.log() << "Implied Performance Time plot has been generated " << currentTime() << std::endl;
    outAdjoint << adjTime;
    outPerform.log() << "Implied Adjoint Performance plot has been generated " << currentTime() << std::endl;
    outSize << tapeSize;
    outPerform.log() << "Tape size dependence on number of independent variables plot has been generated " << currentTime() << std::endl;
    outputImplied(outPerform);
    return ok;
}

// Function outputImplied(cl::AdjointTestOutput & logOut)
// Generating Discount On Swap Rate plot.
// Parameter: logOut is reference to struct with log stream.
void outputImplied(cl::AdjointTestOutput & logOut)
{
    Size steps = Size((MAX_RATE_RANGE - MIN_RATE_RANGE) / RATE_STEP) + 1;
    Rate rate = MIN_RATE_RANGE;
    std::vector<RateDiscount> rDep(steps);

    for (Size i = 0; i < steps; rate += RATE_STEP, i++)
    {
        CommonVars vars({ { 10, Years, rate } });
        Date today = Settings::instance().evaluationDate();
        Date newToday = today + 3 * Years;
        Date newSettlement = vars.calendar.advance(newToday,
            vars.settlementDays, Days);
        Date testDate = newSettlement + 5 * Years;
        boost::shared_ptr<YieldTermStructure> implied(
            new ImpliedTermStructure(Handle<YieldTermStructure>(vars.termStructure),
            newSettlement));
        DiscountFactor baseDiscount = vars.termStructure->discount(newSettlement);
        DiscountFactor discount = vars.termStructure->discount(testDate);
        DiscountFactor impliedDiscount = implied->discount(testDate);

        rDep[i].inputRate_ = rate / 100;
        rDep[i].baseDiscount_ = baseDiscount;
        rDep[i].discount_ = discount;
        rDep[i].impliedDiscount_ = impliedDiscount;
    }
    cl::AdjointTestOutput out("AdjointTermStructure//Implied//output", { { "filename", "discountOnRate" }
                                                                       , { "not_clear", "Not" }
                                                                       , { "line_box_width", "-3" }
                                                                       , { "title", "Discount values dependence on Swap Rate" }
                                                                       , { "ylabel", "Discount value" } });
    out << rDep;
    logOut.log() << "Discount On Swap Rate plot has been generated " << currentTime() << std::endl;
}

// Function outputZSpreaded(cl::AdjointTestOutput & logOut)
// Generating Zero Rate On Swap Rate plot.
// Parameter: logOut is reference to struct with log stream.
void outputZSpreaded(cl::AdjointTestOutput & logOut)
{
    Size steps = Size((MAX_RATE_RANGE - MIN_RATE_RANGE) / RATE_STEP) + 1;
    Rate rate = MIN_RATE_RANGE;
    std::vector<ZeroRateConsistency> rDep(steps);

    for (Size i = 0; i < steps; rate += RATE_STEP, i++)
    {
        CommonVars vars({ { 10, Years, rate } });
        boost::shared_ptr<Quote> me(new SimpleQuote(0.01));
        Handle<Quote> mh(me);
        boost::shared_ptr<YieldTermStructure> spreaded(
            new ZeroSpreadedTermStructure(
            Handle<YieldTermStructure>(vars.termStructure), mh));
        Date testDate = vars.termStructure->referenceDate() + 5 * Years;
        DayCounter rfdc = vars.termStructure->dayCounter();
        Rate zero = vars.termStructure->zeroRate(testDate, rfdc,
            Continuous, NoFrequency);
        Rate spreadedZero = spreaded->zeroRate(testDate, rfdc,
            Continuous, NoFrequency);
        rDep[i].inputRate_ = rate;
        rDep[i].functionValue_ = zero;
    }
    cl::AdjointTestOutput out("AdjointTermStructure//ZSpreaded//output", { { "filename", "zeroRateOnRate" }
                                                                         , { "not_clear", "Not" }
                                                                         , { "title", "Zero Rate dependence on Swap Rate" }
                                                                         , { "ylabel", "Zero Rate" } });
    out << rDep;
    logOut.log() << "Zero Rate On Swap Rate plot has been generated " << currentTime() << std::endl;
}

// Function outputFSpreaded(cl::AdjointTestOutput & logOut)
// Generating Forward Rate On Swap Rate plot.
// Parameter: logOut is reference to struct with log stream.
void outputFSpreaded(cl::AdjointTestOutput & logOut)
{
    Size steps = Size((MAX_RATE_RANGE - MIN_RATE_RANGE) / RATE_STEP) + 1;
    Rate rate = MIN_RATE_RANGE;
    std::vector<ForwardRateConsistency> rDep(steps);

    for (Size i = 0; i < steps; rate += RATE_STEP, i++)
    {
        CommonVars vars({ { 10, Years, rate } });
        boost::shared_ptr<Quote> me(new SimpleQuote(0.01));
        Handle<Quote> mh(me);
        boost::shared_ptr<YieldTermStructure> spreaded(
            new ForwardSpreadedTermStructure(
            Handle<YieldTermStructure>(vars.termStructure), mh));
        Date testDate = vars.termStructure->referenceDate() + 5 * Years;
        DayCounter tsdc = vars.termStructure->dayCounter();
        DayCounter sprdc = spreaded->dayCounter();
        Rate forward = vars.termStructure->forwardRate(testDate, testDate, tsdc,
            Continuous, NoFrequency);
        Rate spreadedForward = spreaded->forwardRate(testDate, testDate, sprdc,
            Continuous, NoFrequency);
        rDep[i].inputRate_ = rate;
        rDep[i].functionValue_ = forward;
    }
    cl::AdjointTestOutput out("AdjointTermStructure//FSpreaded//output", { { "filename", "forwardRateOnRate" }
                                                                         , { "not_clear", "Not" }
                                                                         , { "title", "Forward Rate dependence on Swap Rate" }
                                                                         , { "ylabel", "Forward Rate" } });
    out << rDep;
    logOut.log() << "Forward Rate On Swap Rate plot has been generated " << currentTime() << std::endl;
}

test_suite* AdjointTermStructureTest::suite() {
    test_suite* suite = BOOST_TEST_SUITE("Term structure tests");
    suite->add(QUANTLIB_TEST_CASE(&AdjointTermStructureTest::testZSpreaded));
    suite->add(QUANTLIB_TEST_CASE(&AdjointTermStructureTest::testFSpreaded));
    suite->add(QUANTLIB_TEST_CASE(&AdjointTermStructureTest::testImplied));

    return suite;
}

#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_termstructure)

BOOST_AUTO_TEST_CASE(testConsistencyZeroSpreadedTermStructure)
{
    BOOST_CHECK(AdjointTermStructureTest::testZSpreaded());
}

BOOST_AUTO_TEST_CASE(testConsistencyForwardSpreadedTermStructure)
{
    BOOST_CHECK(AdjointTermStructureTest::testFSpreaded());    
}

BOOST_AUTO_TEST_CASE(testConsistencyYieldTermStructure)
{
    BOOST_CHECK(AdjointTermStructureTest::testImplied());
}
BOOST_AUTO_TEST_SUITE_END()

#endif







