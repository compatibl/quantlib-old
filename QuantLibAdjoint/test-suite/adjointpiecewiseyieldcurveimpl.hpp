/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2005, 2006, 2007, 2008, 2009 StatPro Italia srl
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

// Based on piecewiseyieldcurve.cpp from test-suite.
#ifndef cl_adjoint_piecewice_yied_curve_impl
#define cl_adjoint_piecewice_yied_curve_impl
#pragma once

#include "adjointpiecewiseyieldcurvetest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include <ql/time/schedule.hpp>
#include <boost/make_shared.hpp>
#include <ql/quantlib.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;


#define OUTPUT_FOLDER_NAME "AdjointPiecewiseYieldCurve"


namespace
{
    enum
    {
#if defined CL_GRAPH_GEN
        // Number of points for dependency plots.
        pointNo = 45,
        // Number of points for performance plot.
        iterNo = 45,
        // Step for portfolio size for performance testing .
        step = 1,
#else
        // Number of points for dependency plots.
        pointNo = 20,
        // Number of points for performance plot.
        iterNo = 1,
        // Step for portfolio size for performance testing .
        step = 1,
#endif
        // Defines performance accuracy. Its value is a minimum number
        // of calling of O(1) complexity methods per one performance test.
        iterNumFactor = 10,
        maxSize = 50,
    };

    struct Datum
    {
        Integer n_;
        TimeUnit units_;
        Rate rate_;
    };


    struct SwapData
    {
        SwapData():
        calendar_(Japan())
        , settlementDays_(2)
        , today_(Date(4, October, 2007))
        , settlement_(calendar_.advance(today_, settlementDays_, Days))
        , fixedLegConvention_(Unadjusted)
        , fixedLegFrequency_ (Annual)
        , fixedLegDayCounter_( Thirty360())
        , swapData_(maxSize)
        , rates_(1)
        , instruments_(1)
        , termStructure_()
        {
            Settings::instance().evaluationDate() = today_;

            for (Integer i = 0; i < maxSize; i++)
            {
                swapData_[i] = { i + 1, Years, (4.5 + i*0.03) / 100 };
            }
        }

        // Calculate total NPV of swaps based on input rates using Japanese Yen LIBOR rate.
        // Japanese Yen LIBOR rate can be considered as the interbank cost of borrowing funds in Japanese yens.
        Real calculateSwapNPV(Real swapRate, Integer id)
        {
            // Create ibor index using 6 Month JPY Libor.
            boost::shared_ptr<IborIndex> index = boost::make_shared<JPYLibor>(6 * Months);
            
            rates_[0] = boost::make_shared<SimpleQuote>(swapRate);
            instruments_[0] = boost::make_shared<SwapRateHelper>(Handle<Quote>(rates_[0]), swapData_[id].n_*swapData_[id].units_,
                                                                 calendar_,
                                                                 fixedLegFrequency_, fixedLegConvention_,
                                                                 fixedLegDayCounter_, index);

            // Initialize yield term structure.
            termStructure_ = boost::make_shared<PiecewiseYieldCurve<Discount, LogLinear>>(
                settlement_, instruments_,
                Actual360(),
                1.0e-12);


            RelinkableHandle<YieldTermStructure> curveHandle;
            curveHandle.linkTo(termStructure_);

            // Create ibor index using 6 Month JPY Libor.
            boost::shared_ptr<IborIndex> jpylibor6m = boost::make_shared<JPYLibor>(6 * Months, curveHandle);
            
                Period tenor = swapData_[id].n_*swapData_[id].units_;

            // Create vanilla swap
            VanillaSwap swap = MakeVanillaSwap(tenor, jpylibor6m, 0.0)
                .withEffectiveDate(settlement_)
                .withFixedLegDayCount(fixedLegDayCounter_)
                .withFixedLegTenor(Period(fixedLegFrequency_))
                .withFixedLegConvention(fixedLegConvention_)
                .withFixedLegTerminationDateConvention(fixedLegConvention_)
                .withFixedLegCalendar(calendar_)
                .withFloatingLegCalendar(calendar_);

            return swap.NPV();
        }

        // Global variables.
        Calendar calendar_;
        Natural settlementDays_;
        Date today_, settlement_;
        BusinessDayConvention fixedLegConvention_;
        Frequency fixedLegFrequency_;
        DayCounter fixedLegDayCounter_;

        std::vector<boost::shared_ptr<SimpleQuote> > rates_;
        std::vector<boost::shared_ptr<RateHelper> > instruments_;
        boost::shared_ptr<YieldTermStructure> termStructure_;

        std::vector<Datum> swapData_;

        // Cleanup.
        SavedSettings backup_;
        IndexHistoryCleaner cleaner_;
    };

    // Struct for plots recording.
    struct RateSens
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
            operator << (stream_type& stm, RateSens& v)
        {
                stm << v.rate_
                    << ";" << v.sensitivity_
                    << std::endl;
                return stm;
            }

        Real rate_;
        Real sensitivity_;
    };


    struct TestData
        : public SwapData
    {
        struct Test
        : public cl::AdjointTest<Test>
        {
            Test(Size size, TestData* data)
            : size_(size)
            , data_(data)
            , rates_(size)
            , totalNPV_()
            {
                setLogger(&data_->outPerform_);

                for (Size i = 0; i < size_; i++)
                {
                    rates_[i] = data_->swapData_[i].rate_;
                }

            }

            Size indepVarNumber() { return size_; }

            Size depVarNumber() { return 1; }

            Size minPerfIteration() { return iterNumFactor; }

            void recordTape()
            {
                cl::Independent(rates_);
                calculateNPV();
                f_ = std::make_unique<cl::tape_function<double>>(rates_, totalNPV_);
            }

            // Calculates price of portfolio and each option.
            void calculateNPV()
            {
                totalNPV_.resize(1, 0);
                for (Size i = 0; i < size_; i++)
                    totalNPV_[0] += data_->calculateSwapNPV(rates_[i], i);
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                double h = 1.0e-4;
                analyticalResults_.resize(size_);
                std::vector<Real> npvs(size_);
                for (Size i = 0; i < size_; i++)
                {
                    npvs[i] = data_->calculateSwapNPV(rates_[i], i);
                }
                for (Size i = 0; i < size_; i++)
                {
                    //Evaluate derivative using finite difference.
                    analyticalResults_[i] = (data_->calculateSwapNPV(rates_[i] + h, i) - npvs[i]) / h;
                }
            }

            double relativeTol() const { return 1e-2; }

            double absTol() const { return 1e-4; }

            Size size_;
            TestData* data_;
            std::vector<cl::tape_double> rates_;
            std::vector<cl::tape_double> totalNPV_;
        };

        TestData()
            : SwapData()
            , outPerform_(OUTPUT_FOLDER_NAME
            , { { "filename", "AdjointPerformance" }
            , { "not_clear", "Not" }
            , { "line_box_width", "-5" }
            , { "title", "Swap NPV differentiation performance with respect to rate" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of rates" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        })
            , outAdjoint_(OUTPUT_FOLDER_NAME
            , { { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "title", "Swap NPV adjoint differentiation performance with respect to rate" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of rates" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        })
            , outSize_(OUTPUT_FOLDER_NAME
            , { { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "title", "Tape size dependence on number of rates" }
            , { "ylabel", "Size (MB)" }
            , { "xlabel", "Number of rates" }
            , { "cleanlog", "false" }
        })
            , out_(OUTPUT_FOLDER_NAME "//output", {
            { "filename", "StrikeDependence" }
            , { "not_clear", "Not" }
            , { "title", "Swap NPV sensitivity dependence on rate" }
            , { "ylabel", "Swap NPV" }
            , { "xlabel", "Rate" }
            , { "cleanlog", "false" }
        })
        {
        }

        bool makeOutput()
        {
            bool ok = true;
            if (pointNo > 0)
            {
                ok &= recordDependencePlot();
            }
            ok &= cl::recordPerformance(*this, iterNo, step);
            return ok;
        }

        std::shared_ptr<Test> getTest(size_t size)
        {
            return std::make_shared<Test>(size, this);
        }

        // Makes plots for strike sensitivity dependence.
        bool recordDependencePlot()
        {
            std::vector<RateSens> outData(pointNo);
            auto test = getTest(pointNo);
            bool ok = test->testReverse();
            for (Size i = 0; i < pointNo; i++)
            {
                outData[i] = { test->rates_[i], test->data_->calculateSwapNPV(test->rates_[i], i) };
            }
            out_ << outData;
            return ok;
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output out_;
    };

    typedef TestData::Test PiecewiceYieldCurveTest;
}
#endif