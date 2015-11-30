/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2005, 2007 StatPro Italia srl
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

// Based on bermudanswaption.cpp file from test-suite.


#ifndef cl_adjoint_berbudan_swaption_impl_hpp
#define cl_adjoint_berbudan_swaption_impl_hpp
#pragma once


#include "adjointbermudanswaptiontest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include <ql/instruments/swaption.hpp>
#include <ql/pricingengines/swaption/treeswaptionengine.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>
#include <ql/pricingengines/swaption/fdhullwhiteswaptionengine.hpp>
#include <ql/models/shortrate/onefactormodels/hullwhite.hpp>
#include <ql/cashflows/coupon.hpp>
#include <ql/time/daycounters/thirty360.hpp>
#include <ql/indexes/ibor/euribor.hpp>
#include <ql/time/schedule.hpp>
#include <boost/make_shared.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;


#define OUTPUT_FOLDER_NAME "AdjointBermudanSwaption"


namespace
{
    enum
    {
#if defined CL_GRAPH_GEN
        // Number of points for dependency plots.
        pointNo = 100,
        // Number of points for performance plot.
        iterNo = 30,
        // Step for portfolio size for performance testing .
        step = 1,
#else
        // Number of points for dependency plots.
        pointNo = 0,
        // Number of points for performance plot.
        iterNo = 0,
        // Step for portfolio size for performance testing .
        step = 1,
#endif
        // Defines performance accuracy. Its value is a minimum number
        // of calling of O(1) complexity methods per one performance test.
        iterNumFactor = 0,
    };


    enum EngineType { tree, fdHullWhite };

    std::string toString(EngineType engine)
    {
        switch (engine)
        {
        case tree:
            return "TreeSwaptionEngine";
        case fdHullWhite:
            return "FdHullWhiteSwaptionEngine";
        default:
            return "";
        }
    }


    struct CommonVars
    {
        CommonVars()
        : backup_()
        , startYears_(1)
        , length_(5)
        , type_(VanillaSwap::Payer)
        , nominal_(1000.0)
        , fixedConvention_(Unadjusted)
        , floatingConvention_(ModifiedFollowing)
        , fixedFrequency_(Annual)
        , floatingFrequency_(Semiannual)
        , fixedDayCount_(Thirty360())
        , termStructure_()
        , index_(boost::shared_ptr<IborIndex>(new Euribor6M(termStructure_)))
        , settlementDays_(2)
        , calendar_(index_->fixingCalendar())
        , today_(calendar_.adjust(Date::todaysDate()))
        , settlement_(calendar_.advance(today_, settlementDays_, Days))
        , fixedRate_(0.4)
        {
        }

        // Provides swaps which are different only by their fixed rate strike.
        boost::shared_ptr<VanillaSwap> makeSwap(Spread spread)
        {
            Date start = calendar_.advance(settlement_, startYears_, Years);
            Date maturity = calendar_.advance(start, length_, Years);
            Schedule fixedSchedule(start, maturity,
                                   Period(fixedFrequency_),
                                   calendar_,
                                   fixedConvention_,
                                   fixedConvention_,
                                   DateGeneration::Forward, false);
            Schedule floatSchedule(start, maturity,
                                   Period(floatingFrequency_),
                                   calendar_,
                                   floatingConvention_,
                                   floatingConvention_,
                                   DateGeneration::Forward, false);
            boost::shared_ptr<VanillaSwap> swap(
                new VanillaSwap(type_, nominal_,
                fixedSchedule, fixedRate_, fixedDayCount_,
                floatSchedule, index_, spread,
                index_->dayCounter()));
            swap->setPricingEngine(boost::shared_ptr<PricingEngine>(
                new DiscountingSwapEngine(termStructure_)));
            return swap;
        }

        // cleanup
        SavedSettings backup_;

        // underlying swap parameters
        Integer startYears_;
        Integer length_;
        VanillaSwap::Type type_;
        Real nominal_;
        BusinessDayConvention fixedConvention_;
        BusinessDayConvention floatingConvention_;
        Frequency fixedFrequency_;
        Frequency floatingFrequency_;
        DayCounter fixedDayCount_;
        RelinkableHandle<YieldTermStructure> termStructure_;
        boost::shared_ptr<IborIndex> index_;
        Natural settlementDays_;

        // global data
        Calendar calendar_;
        Date today_;
        Date settlement_;
        Rate fixedRate_;
    };

    struct SwaptionsData
        : public CommonVars
    {
        explicit SwaptionsData(EngineType engine)
        : CommonVars()
        , engineId_(engine)
        {
            today_ = Date(15, February, 2002);
            Settings::instance().evaluationDate() = today_;
            settlement_ = Date(19, February, 2002);

            // flat yield term structure impling 1x5 swap at 5%
            termStructure_.linkTo(flatRate(settlement_,
                0.04875825,
                Actual365Fixed()));
            atmRate_ = makeSwap(0.0)->fairRate();
            Real a = 0.048696;
            Real sigma = 0.0058904;
            boost::shared_ptr<VanillaSwap> atmSwap = makeSwap(atmRate_);
            const Leg& leg = atmSwap->fixedLeg();
            model_ = boost::make_shared<HullWhite>(termStructure_, a, sigma);
            std::vector<Date> exerciseDates(leg.size());
            for (Size i = 0; i < leg.size(); i++)
            {
                boost::shared_ptr<Coupon> coupon =
                    boost::dynamic_pointer_cast<Coupon>(leg[i]);
                exerciseDates[i] = coupon->accrualStartDate();
            }
            exercise_ = boost::make_shared<BermudanExercise>(exerciseDates);
            switch (engine)
            {
                case tree:
                    engine_ = boost::make_shared<TreeSwaptionEngine>(model_, 200);
                    break;
                case fdHullWhite:
                    engine_ = boost::make_shared<FdHullWhiteSwaptionEngine>(model_, 100, 100, 0, 1e-7);
                    break;
            }
        }

        // Returns net present value (NPV) of swaption with given fixed rate strike and other
        // parameters determined from class instance.
        Real swaptionNVP(Spread spread)
        {
            auto swap = makeSwap(spread);
            Swaption swaption(swap, exercise_);
            swaption.setPricingEngine(engine_);
            return swaption.NPV();
        }

        Rate atmRate_;
        EngineType engineId_;
        boost::shared_ptr<HullWhite> model_;
        boost::shared_ptr<Exercise> exercise_;
        boost::shared_ptr<PricingEngine> engine_;
    };


    // Struct for plots recording.
    struct StrikeSens
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "Spread", ""
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, StrikeSens& v)
        {
                stm << v.strike_
                    << ";" << v.sensitivity_
                    << std::endl;
                return stm;
            }

        Real strike_;
        Real sensitivity_;
    };


    struct SwaptionTest
        : public cl::AdjointTest<SwaptionTest>
    {
        SwaptionTest(Size size, SwaptionsData* data, cl::tape_empty_test_output* logger = nullptr)
        : AdjointTest()
        , size_(size)
        , data_(data)
        , spread_(size)
        , prices_()
        , totalPrice_()
        {
            setLogger(logger);

            for (Size i = 0; i < size_; i++)
            {
                spread_[i] = 0.300909 + 0.1 * i / (size_ - 1);
            }
        }

        Size indepVarNumber() { return size_; }

        Size minPerfIteration() { return iterNumFactor; }

        void recordTape()
        {
            cl::Independent(spread_);
            calculatePrices();
            f_ = std::make_unique<cl::tape_function<double>>(spread_, totalPrice_);
        }

        // Calculates price of portfolio and each option.
        void calculatePrices()
        {
            totalPrice_.resize(1, 0);
            prices_.resize(size_);
            for (Size i = 0; i < size_; i++)
            {
                prices_[i] = data_->swaptionNVP(spread_[i]);
                totalPrice_.front() += prices_[i];
            }
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            Real h = 1e-6 ;  // shift for finite diff. method
            analyticalResults_.resize(size_);
            for (Size i = 0; i < size_; i++)
            {
                analyticalResults_[i] = (data_->swaptionNVP(spread_[i] + h)
                                         - data_->swaptionNVP(spread_[i] - h)) / (2 * h);
            }
        }

        double relativeTol() const { return 1e-2; }

        double absTol() const { return 1e-10; }

        Size size_;
        SwaptionsData* data_;
        std::vector<cl::tape_double> spread_;
        std::vector<Real> prices_;
        std::vector<cl::tape_double> totalPrice_;
    };


    struct TestData
        : public SwaptionsData
    {
        explicit TestData(EngineType engine)
        : SwaptionsData(engine)
        , outPerform_(OUTPUT_FOLDER_NAME "//" + toString(engine), {
            { "title", "Bermudan swaption NPV differentiation performance with respect to spread" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of spreads" }
            , { "line_box_width", "-5" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        })
            , outAdjoint_(OUTPUT_FOLDER_NAME "//" + toString(engine), {
                { "title", "Swaption NPV adjoint differentiation with respect to spread" }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of spreads" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        })
            , outSize_(OUTPUT_FOLDER_NAME "//" + toString(engine), {
                { "title", "Tape size dependence on number of spreads" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
            , { "cleanlog", "false" }
        })
            , out_(OUTPUT_FOLDER_NAME "//" + toString(engine) + "//output", {
                { "filename", "NVPonSpread" }
            , { "not_clear", "Not" }
            , { "title", "Swaption NPV dependence on spread" }
            , { "ylabel", "Swaption NPV" }
            , { "xlabel", "Spread" }
            , { "cleanlog", "false" }
        })
        {
        }

        bool makeOutput()
        {
            bool ok = true;
            if (pointNo > 0)
            {
                recordDependencePlot(pointNo);
            }
            ok &= cl::recordPerformance(*this, iterNo, step);
            return ok;
        }

        std::shared_ptr<SwaptionTest> getTest(size_t size)
        {
            return std::make_shared<SwaptionTest>(size, this, &outPerform_);
        }

        // Makes plots for strike sensitivity dependence.
        void recordDependencePlot(size_t pointNo)
        {
            std::vector<StrikeSens> outData(pointNo);
            auto test = getTest(pointNo);
            test->recordTape();
            for (Size i = 0; i < pointNo; i++)
            {
                outData[i] = { test->spread_[i], test->prices_[i] };
            }
            out_ << outData;
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output out_;
    };
}

#endif