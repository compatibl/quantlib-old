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

// Based on swaptionvolatilitymatrix.cpp project from Quantlib.

#ifndef cl_adjoint_default_probability_curve_impl_hpp
#define cl_adjoint_default_probability_curve_impl_hpp
#pragma once

#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include <boost/make_shared.hpp>
#include "adjointdefaultprobabilitycurvetest.hpp"
#include "utilities.hpp"
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
using namespace boost::unit_test_framework;


#define OUTPUT_FOLDER_NAME "AdjointSwaptionVolatilityMatrix"


namespace
{
    struct CommonVars
    {
        CommonVars() :
        dayCounter_(Actual360())
        , calendar_(TARGET())
        , today_(Settings::instance().evaluationDate())
        , startDate_(today_)
        , endDate_(startDate_)
        , notional_(1.0)
        , frequency_(Quarterly)
        , convention_(Following)
        , rule_(DateGeneration::TwentiethIMM)
        , settlementDays_(1)
        , recoveryRate_(0.4)
        , fixedRate_(0.05)
        {
            discountCurve_.linkTo(boost::shared_ptr<YieldTermStructure>(
                new FlatForward(today_, 0.06, Actual360())));
        }

        DayCounter dayCounter_;
        Calendar calendar_;

        Date today_;
        Date startDate_;
        Date endDate_;

        Real notional_;
        Frequency frequency_;
        BusinessDayConvention convention_;
        DateGeneration::Rule rule_;
        Integer settlementDays_;
        Real recoveryRate_;
        Real fixedRate_;
        RelinkableHandle<YieldTermStructure> discountCurve_;

    };

    template <class T, class I>
    struct FairRate : CommonVars
    {
        FairRate() :
        CommonVars()
        {
        }

        template <class T, class I>
        Real computeFairRateFromSpread(Real quote, Integer n)
        {
            settlementDays_ = 1;
            std::vector<boost::shared_ptr<DefaultProbabilityHelper> > helpers = {
                boost::make_shared<SpreadCdsHelper>(quote, Period(n, Months),
                settlementDays_, calendar_,
                frequency_, convention_, rule_,
                dayCounter_, recoveryRate_,
                discountCurve_) };


            RelinkableHandle<DefaultProbabilityTermStructure> piecewiseCurve;
            piecewiseCurve.linkTo(
                boost::make_shared<PiecewiseDefaultCurve<T, I>>(
                today_, helpers,
                Thirty360()));

            // ensure apple-to-apple comparison
            SavedSettings backup;
            Settings::instance().includeTodaysCashFlows() = true;


            Date protectionStart = today_ + settlementDays_;
            startDate_ = calendar_.adjust(protectionStart, convention_);
            endDate_ = today_ + n * Months;

            Schedule schedule(startDate_, endDate_, Period(frequency_), calendar_,
                              convention_, Unadjusted, rule_, false);

            CreditDefaultSwap cds(Protection::Buyer, notional_, quote,
                                  schedule, convention_, dayCounter_,
                                  true, true, protectionStart);
            cds.setPricingEngine(boost::shared_ptr<PricingEngine>(
                new MidPointCdsEngine(piecewiseCurve, recoveryRate_,
                discountCurve_)));

            return cds.fairSpread();
        }

        template <class T, class I>
        Real computeFairRateFromUpfront(Real quote, Integer n)
        {
            settlementDays_ = 0;
            std::vector<boost::shared_ptr<DefaultProbabilityHelper> > helpers = {
                boost::make_shared<UpfrontCdsHelper>(quote, fixedRate_, Period(n, Months),
                settlementDays_, calendar_,
                frequency_, convention_, rule_,
                dayCounter_, recoveryRate_,
                discountCurve_) };

            RelinkableHandle<DefaultProbabilityTermStructure> piecewiseCurve;
            piecewiseCurve.linkTo(
                boost::make_shared<PiecewiseDefaultCurve<T, I>>(today_, helpers,
                Thirty360()));

            // ensure apple-to-apple comparison
            SavedSettings backup;
            Settings::instance().includeTodaysCashFlows() = true;


            Date protectionStart = today_ + settlementDays_;
            startDate_ = calendar_.adjust(protectionStart, convention_);
            endDate_ = today_ + n * Months;

            Schedule schedule(startDate_, endDate_, Period(frequency_), calendar_,
                              convention_, Unadjusted, rule_, false);

            CreditDefaultSwap cds(Protection::Buyer, notional_, quote, fixedRate_,
                                  schedule, convention_, dayCounter_,
                                  true, true, protectionStart);
            cds.setPricingEngine(boost::shared_ptr<PricingEngine>(
                new MidPointCdsEngine(piecewiseCurve, recoveryRate_,
                discountCurve_, true)));

            return cds.fairUpfront();
        }
    };

    struct TestFlatHazardRate
        : public cl::AdjointTest<TestFlatHazardRate>
    {
        TestFlatHazardRate(Size size)
            : size_(size)
            , hazardRate_(1, 0.0100)
            , t_(size)
            , computedProbability_(size)
            , var_()
        {
            for (Size i = 0; i < size_; i++)
            {
                var_.endDate_ = var_.calendar_.advance(var_.endDate_, 1, Months);
                t_[i] = var_.dayCounter_.yearFraction(var_.startDate_, var_.endDate_);
            }
        }

        Size indepVarNumber() { return 1; }

        Size depVarNumber() { return size_; }

        void computeProbability()
        {
            Handle<Quote> hazardRateQuote = Handle<Quote>(
                boost::make_shared<SimpleQuote>(hazardRate_[0]));
            FlatHazardRate flatHazardRate(var_.today_, hazardRateQuote, var_.dayCounter_);

            for (Size i = 0; i < size_; i++)
            {
                computedProbability_[i] = flatHazardRate.defaultProbability(t_[i]);
            }
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            analyticalResults_.resize(size_, 0);
            for (Size i = 0; i < size_; i++)
            {
                analyticalResults_[i] = t_[i] * std::exp(-hazardRate_[0] * t_[i]);
            }
        }

        double relativeTol() const { return 1e-10; }

        double absTol() const { return 1e-10; }

        Size size_;
        std::vector<cl::tape_double> hazardRate_;
        std::vector<Time> t_;
        std::vector<cl::tape_double> computedProbability_;

        CommonVars var_;
    };

    struct TestFlatHazardRateTime
        : public cl::AdjointTest<TestFlatHazardRateTime>
    {
        static const CalcMethod default_method = other;

        TestFlatHazardRateTime(Size size)
            : size_(size)
            , hazardRate_(1, 0.0100)
            , t_(size)
            , doubleT_(size)
            , computedProbability_(size)
            , var_()
        {
            for (Size i = 0; i < size_; i++)
            {
                var_.endDate_ = var_.calendar_.advance(var_.endDate_, 1, Months);
                t_[i] = var_.dayCounter_.yearFraction(var_.startDate_, var_.endDate_);
                doubleT_[i] = double(var_.dayCounter_.yearFraction(var_.startDate_, var_.endDate_));
            }
        }

        Size indepVarNumber() { return size_; }

        Size depVarNumber() { return size_; }

        void computeProbability()
        {
            Handle<Quote> hazardRateQuote = Handle<Quote>(
                boost::make_shared<SimpleQuote>(hazardRate_[0]));
            FlatHazardRate flatHazardRate(var_.today_, hazardRateQuote, var_.dayCounter_);

            for (Size i = 0; i < size_; i++)
            {
                computedProbability_[i] = flatHazardRate.defaultProbability(t_[i]);
            }
        }

        // Calculates derivatives using finite difference method.

        void calcAnalytical()
        {
            double h = 1e-4;  // shift for finite diff. method
            analyticalResults_.resize(size_*size_, 0);

            Handle<Quote> hazardRateQuote = Handle<Quote>(
                boost::make_shared<SimpleQuote>(0.0100));
            FlatHazardRate flatHazardRate(var_.today_, hazardRateQuote, var_.dayCounter_);

            for (Size i = 0; i < size_; i++)
            {
                Real rightValue = flatHazardRate.defaultProbability(t_[i] + h);
                Real leftValue = flatHazardRate.defaultProbability(t_[i] - h);
                analyticalResults_[i*size_ + i] = (rightValue - leftValue) / (2 * h);
            }
        }

        double relativeTol() const { return 1e-4; }

        double absTol() const { return 1e-4; }

        Size size_;
        std::vector<cl::tape_double> hazardRate_;
        std::vector<cl::tape_double> t_;
        std::vector<double> doubleT_;
        std::vector<cl::tape_double> computedProbability_;


        CommonVars var_;
    };

    template <class T, class I>
    struct TestBootstrapFromSpread
        : public cl::AdjointTest<TestBootstrapFromSpread<T,I>>
        {
        static const CalcMethod default_method = other;

        TestBootstrapFromSpread(Size size)
            : size_(size)
            , quote_(size)
            , doubleQuote_(size)
            , computedRate_(size)
            , var_()
        {
            for (Size i = 0; i < size_; i++)
            {
                doubleQuote_[i] = (i + 1)*0.001;
                quote_[i] = (i + 1)*0.001;
            }
        }

            Size indepVarNumber() { return size_; }

            Size depVarNumber() { return size_; }

            void computeFairRate()
            {
                for (Size i = 0; i < size_; i++)
                    computedRate_[i] = var_.computeFairRateFromSpread<T, I>(quote_[i], i+1);
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                double h = 1e-4;  // shift for finite diff. method
                analyticalResults_.resize(size_*size_, 0);
                for (Size i = 0; i < size_; i++)
                {
                    Real rightValue = var_.computeFairRateFromSpread<T, I>(quote_[i] + h, i + 1);
                    Real leftValue = var_.computeFairRateFromSpread<T, I>(quote_[i] - h, i + 1);
                    analyticalResults_[i*size_ + i] = (rightValue - leftValue) / (2 * h);
                }
            }

            double relativeTol() const { return 1e-4; }

            double absTol() const { return 1e-4; }

            Size size_;
            std::vector<cl::tape_double> quote_;
            std::vector<double> doubleQuote_;
            std::vector<cl::tape_double> computedRate_;


            FairRate<T, I> var_;
        };


    template <class T, class I>
    struct TestBootstrapFromUpfront
        : public cl::AdjointTest<TestBootstrapFromUpfront<T, I>>
    {
        static const CalcMethod default_method = other;

        TestBootstrapFromUpfront(Size size)
            : size_(size)
            , quote_(size)
            , doubleQuote_(size)
            , computedRate_(size)
            , var_()
        {
            for (Size i = 0; i < size_; i++)
            {
                doubleQuote_[i] = (i + 1)*0.001;
                quote_[i] = (i + 1)*0.001;
            }
        }

        Size indepVarNumber() { return size_; }

        Size depVarNumber() { return size_; }

        void computeFairRate()
        {
            for (Size i = 0; i < size_; i++)
                computedRate_[i] = var_.computeFairRateFromUpfront<T, I>(quote_[i], i + 2);
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            double h = 1e-4;  // shift for finite diff. method
            analyticalResults_.resize(size_*size_, 0);
            for (Size i = 0; i < size_; i++)
            {
                Real rightValue = var_.computeFairRateFromUpfront<T, I>(quote_[i] + h, i + 2);
                Real leftValue = var_.computeFairRateFromUpfront<T, I>(quote_[i] - h, i + 2);
                analyticalResults_[i*size_ + i] = (rightValue - leftValue) / (2 * h);
            }
        }

        double relativeTol() const { return 1e-4; }

        double absTol() const { return 1e-4; }

        Size size_;
        std::vector<cl::tape_double> quote_;
        std::vector<double> doubleQuote_;
        std::vector<cl::tape_double> computedRate_;


        FairRate<T, I> var_;
    };
}
#endif