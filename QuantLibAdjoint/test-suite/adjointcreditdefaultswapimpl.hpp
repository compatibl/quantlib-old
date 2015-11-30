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

#ifndef cl_adjoint_creditdefaultswap_impl_hpp
#define cl_adjoint_creditdefaultswap_impl_hpp
#pragma once

#include <ql/quantlib.hpp>
#include "utilities.hpp"
#include "adjointcreditdefaultswaptest.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include <boost/shared_ptr.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;


#define OUTPUT_FOLDER_NAME "AdjointCreditDefaultSwap"

namespace
{
    // Support structure for plotting
    struct RealPair
    {
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

        Real x, y;
    };

    struct SwapData
    {
        SwapData()
        : calendar_()
        , discountCurve_()
        , curveDayCounter_()
        , backup_()
        , discountDates_()
        , discountFactors_()
        , defaultProbabilities_()
        {
            Settings::instance().evaluationDate() = Date(9, June, 2006);
            evalDate_ = Settings::instance().evaluationDate();
        }

        void setDefaultProbabilitiesData()
        {
            calendar_ = UnitedStates();

            discountDates_ = {
                evalDate_,
                calendar_.advance(evalDate_, 1, Weeks, ModifiedFollowing),
                calendar_.advance(evalDate_, 1, Months, ModifiedFollowing),
                calendar_.advance(evalDate_, 2, Months, ModifiedFollowing),
                calendar_.advance(evalDate_, 3, Months, ModifiedFollowing),
                calendar_.advance(evalDate_, 6, Months, ModifiedFollowing),
                calendar_.advance(evalDate_, 1, Years, ModifiedFollowing),
                calendar_.advance(evalDate_, 2, Years, ModifiedFollowing),
                calendar_.advance(evalDate_, 3, Years, ModifiedFollowing),
                calendar_.advance(evalDate_, 4, Years, ModifiedFollowing),
                calendar_.advance(evalDate_, 5, Years, ModifiedFollowing),
                calendar_.advance(evalDate_, 6, Years, ModifiedFollowing),
                calendar_.advance(evalDate_, 7, Years, ModifiedFollowing),
                calendar_.advance(evalDate_, 8, Years, ModifiedFollowing),
                calendar_.advance(evalDate_, 9, Years, ModifiedFollowing),
                calendar_.advance(evalDate_, 10, Years, ModifiedFollowing),
                calendar_.advance(evalDate_, 15, Years, ModifiedFollowing)
            };

            discountFactors_ = {
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

            // Build the discount curve
            curveDayCounter_ = Actual360();
            discountCurve_.linkTo(boost::make_shared<DiscountCurve>(discountDates_
                                                                    , discountFactors_
                                                                    , curveDayCounter_));

        }

        std::vector<Real> calculateDefaultProbabilitiesFunction(std::vector<Real> defaultProbabilities)
        {
            Size size = defaultProbabilities.size();

            std::vector<Date> dates(size);
            for (Size i = 0; i < size; i++)
                dates[i] = evalDate_ + i * Weeks;

            std::vector<Real> hazardRates;
            DayCounter dayCounter = Thirty360();
            hazardRates.push_back(0.0);

            for (Size i = 1; i < size; ++i)
            {
                Time t1 = dayCounter.yearFraction(dates[0], dates[i - 1]);
                Time t2 = dayCounter.yearFraction(dates[0], dates[i]);
                Probability S1 = 1.0 - defaultProbabilities[i - 1];
                Probability S2 = 1.0 - defaultProbabilities[i];
                hazardRates.push_back(std::log(S1 / S2) / (t2 - t1));
            }

            // Build the hazard rates structure
            RelinkableHandle<DefaultProbabilityTermStructure> piecewiseFlatHazardRate;
            piecewiseFlatHazardRate.linkTo(
                boost::make_shared<InterpolatedHazardRateCurve<BackwardFlat>>(dates, hazardRates, Thirty360()));

            // Build the schedule
            Date issueDate(20, March, 2006);
            Date maturity = issueDate + size * Weeks;
            Frequency cdsFrequency = Weekly;
            BusinessDayConvention cdsConvention = ModifiedFollowing;
            Schedule schedule(issueDate
                              , maturity
                              , Period(cdsFrequency)
                              , calendar_
                              , cdsConvention
                              , cdsConvention
                              , DateGeneration::Forward
                              , false);

            // Build the credit default swap
            Real recoveryRate = 0.25;
            Rate fixedRate = 0.0224;
            DayCounter dayCount = Actual360();
            Real cdsNotional = 100.0;
            CreditDefaultSwap cds(Protection::Seller
                                  , cdsNotional
                                  , fixedRate
                                  , schedule
                                  , cdsConvention
                                  , dayCount
                                  , true
                                  , true);

            cds.setPricingEngine(
                boost::make_shared<MidPointCdsEngine>(piecewiseFlatHazardRate, recoveryRate, discountCurve_));

            // Compute NPV and fair spread of the credit default swap
            return std::vector<Real>{cds.NPV(), cds.fairSpread() };
        }

        void setDiscountFactorData()
        {
            calendar_ = UnitedStates();

            // Set dates for default probabilities
            discountDates_ = {
                evalDate_,
                calendar_.advance(evalDate_, 6, Months, ModifiedFollowing),
                calendar_.advance(evalDate_, 1, Years, ModifiedFollowing),
                calendar_.advance(evalDate_, 2, Years, ModifiedFollowing),
                calendar_.advance(evalDate_, 3, Years, ModifiedFollowing),
                calendar_.advance(evalDate_, 4, Years, ModifiedFollowing),
                calendar_.advance(evalDate_, 5, Years, ModifiedFollowing),
                calendar_.advance(evalDate_, 7, Years, ModifiedFollowing),
                calendar_.advance(evalDate_, 10, Years, ModifiedFollowing),
                calendar_.advance(evalDate_, 12, Years, ModifiedFollowing)
            };

            defaultProbabilities_ = {
                0.0000, 0.0047, 0.0093, 0.0286, 0.0619, 0.0953, 0.1508, 0.2288, 0.3666, 0.5 };

            curveDayCounter_ = Actual360();
        }

        std::vector<Real> getDiscountFactors(std::vector<Real> dfs)
        {
            int size = dfs.size() + 1;
            std::vector<Date> discountDates(size);
            for (Size i = 0; i < size; i++)
                discountDates[i] = evalDate_ + i * Weeks;

            std::vector<DiscountFactor> discountFactors;
            discountFactors.push_back(1.0);
            discountFactors.insert(discountFactors.end(), dfs.begin(), dfs.end());
            return discountFactors;
        }

        std::vector<Date> getDiscountDates(int size)
        {
            std::vector<Date> discountDates(size);
            for (Size i = 0; i < size; i++)
                discountDates[i] = evalDate_ + i * Weeks;
            return discountDates;
        }

        RelinkableHandle<DefaultProbabilityTermStructure> getHazardRates(std::vector<Date> discountDates)
        {
            std::vector<Real> hazardRates;
            hazardRates.push_back(0.0);
            DayCounter dayCounter = Thirty360();
            for (Size i = 1; i < discountDates_.size(); ++i)
            {
                Time t1 = dayCounter.yearFraction(discountDates_[0], discountDates_[i - 1]);
                Time t2 = dayCounter.yearFraction(discountDates_[0], discountDates_[i]);
                Probability S1 = 1.0 - defaultProbabilities_[i - 1];
                Probability S2 = 1.0 - defaultProbabilities_[i];
                hazardRates.push_back(std::log(S1 / S2) / (t2 - t1));
            }

            RelinkableHandle<DefaultProbabilityTermStructure> piecewiseFlatHazardRate;
            piecewiseFlatHazardRate.linkTo(
                boost::make_shared<InterpolatedHazardRateCurve<BackwardFlat>>(discountDates_, hazardRates, Thirty360()));
            return piecewiseFlatHazardRate;
        }

        Schedule getSchedule(int size)
        {
            Date issueDate(20, March, 2006);
            Date maturity = issueDate + size * Weeks;
            Frequency cdsFrequency = Weekly;
            BusinessDayConvention cdsConvention = ModifiedFollowing;
            return Schedule(issueDate
                              , maturity
                              , Period(cdsFrequency)
                              , calendar_
                              , cdsConvention
                              , cdsConvention
                              , DateGeneration::Forward
                              , false);
        }

        std::vector<Real> calculateDiscountFactorFunction(std::vector<Real> dfs)
        {
            Size size = dfs.size() + 1;

            std::vector<Date> discountDates = getDiscountDates(size);

            // Build the discount curve
            discountCurve_.linkTo(boost::make_shared<DiscountCurve>(discountDates, 
                getDiscountFactors(dfs), curveDayCounter_));

            // Set hazard rates
            RelinkableHandle<DefaultProbabilityTermStructure> piecewiseFlatHazardRate = getHazardRates(discountDates);

            // Build the credit default swap
            Real recoveryRate = 0.25;
            Rate fixedRate = 0.0224;
            DayCounter dayCount = Actual360();
            Real cdsNotional = 100.0;
            CreditDefaultSwap cds(Protection::Seller
                                  , cdsNotional
                                  , fixedRate
                                  , getSchedule(size)
                                  , BusinessDayConvention::ModifiedFollowing
                                  , dayCount
                                  , true
                                  , true);

            cds.setPricingEngine(
                boost::make_shared<MidPointCdsEngine>(piecewiseFlatHazardRate, recoveryRate, discountCurve_));

            // Compute NPV and fair spread of the credit default swap
            return std::vector<Real>{cds.NPV(), cds.fairSpread() };
        }

        std::vector<Real> calculateNotionalSpreadRateFunction(std::vector<Real> parameters)
        {
            Date today = Settings::instance().evaluationDate();
            calendar_ = TARGET();

            // Build the default probability curve
            RelinkableHandle<DefaultProbabilityTermStructure> probabilityCurve;
            Handle<Quote> hazardRate = Handle<Quote>(boost::make_shared<SimpleQuote>(0.01234));
            probabilityCurve.linkTo(
                boost::make_shared<FlatHazardRate>(0, calendar_, hazardRate, Actual360()));

            // Build the discount curve
            discountCurve_.linkTo(boost::make_shared<FlatForward>(today, 0.06, Actual360()));

            // Build the schedule
            Date issueDate = calendar_.advance(today, -1, Years);
            Date maturity = calendar_.advance(issueDate, 10, Years);
            Frequency frequency = Semiannual;
            BusinessDayConvention convention = ModifiedFollowing;
            Schedule schedule(issueDate
                              , maturity
                              , Period(frequency)
                              , calendar_
                              , convention
                              , convention
                              , DateGeneration::Forward
                              , false);


            // Build the credit default swap
            Real notional = parameters[0];
            Rate spread = parameters[1];
            Real recoveryRate = parameters[2];
            DayCounter dayCount = Actual360();
            CreditDefaultSwap cds(Protection::Seller
                                  , notional
                                  , spread
                                  , schedule
                                  , convention
                                  , dayCount
                                  , true
                                  , true);

            // Compute NPV and fair spread of the credit default swap for the first pricing engine
            auto midPointCdsEngine = boost::make_shared<MidPointCdsEngine>(probabilityCurve, recoveryRate, discountCurve_);
            cds.setPricingEngine(midPointCdsEngine);
            std::vector<Real> function(2);
            function[0] = cds.NPV();

            // Compute NPV and fair spread of the credit default swap for the second pricing engine
            auto integralCdsEngine = boost::make_shared<IntegralCdsEngine>(1 * Days, probabilityCurve, recoveryRate, discountCurve_);
            cds.setPricingEngine(integralCdsEngine);
            function[1] = cds.NPV();
            return function;
        }

        template <class Func>
        void calculateFinDiff(std::vector<Real>& X
                              , Real h
                              , std::vector<Real>& sf_Finite
                              , Func function)
        {
            Size sizeX = X.size();

            std::vector<Real> totalfunction = function(X);

            Size sizeY = totalfunction.size();

            sf_Finite.resize(sizeX * sizeY);

            for (Size i = 0; i < sizeX; i++)
            {
                X[i] += h;
                std::vector<Real> finiteTotalfunction = function(X);
                X[i] -= h;

                for (Size j = 0; j < sizeY; j++)
                {
                    sf_Finite[sizeX*j + i] = (finiteTotalfunction[j] - totalfunction[j]) / h;
                }
            }
        }

        SavedSettings backup_;

        Date evalDate_;
        Calendar calendar_;

        std::vector<Date> discountDates_;
        std::vector<DiscountFactor> discountFactors_;
        std::vector<Probability> defaultProbabilities_;

        RelinkableHandle<YieldTermStructure> discountCurve_;
        DayCounter curveDayCounter_;
    };

    struct DefaultProbabilitiesTestData
        : public SwapData
    {
        struct Test
        : public cl::AdjointTest<Test>
        {
            static const CalcMethod default_method = other;

            Test(Size size, DefaultProbabilitiesTestData* data)
            : size_(size)
            , data_(data)
            , probability_(size)
            , doubleProbability_(size)
            , iterNumFactor_(1000)
            {
                setLogger(&data_->outPerform_);

                for (Size i = 0; i < size_; i++)
                    doubleProbability_[i] = (double)(0.4 * i / size_);

                std::copy(doubleProbability_.begin(), doubleProbability_.end(), probability_.begin());

                data->setDefaultProbabilitiesData();
            }

            Size indepVarNumber() { return size_; }

            Size depVarNumber() { return 2; }

            Size minPerfIteration() { return iterNumFactor_; }

            void recordTape()
            {
                cl::Independent(probability_);
                calculateFunction();
                f_ = std::make_unique<cl::tape_function<double>>(probability_, calculatedFunction_);
            }

            void calculateFunction()
            {
                calculatedFunction_ = data_->calculateDefaultProbabilitiesFunction(probability_);
            }

            // Calculates derivatives using adjoint.
            void calcAdjoint()
            {
                adjointResults_ = f_->Jacobian(doubleProbability_);
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                double h = 1.0e-6;  // shift for finite diff. method
                analyticalResults_.resize(size_);
                data_->calculateFinDiff(probability_
                                        , h
                                        , analyticalResults_
                                        , [this] (std::vector<Real>& v) -> std::vector<Real>
                {
                    return data_->calculateDefaultProbabilitiesFunction(v);
                });
            }

            double relativeTol() const { return 1e-4; }

            double absTol() const { return 1e-5; }

            Size size_;
            Size iterNumFactor_;
            DefaultProbabilitiesTestData* data_;
            std::vector<cl::tape_double> probability_;
            std::vector<double> doubleProbability_;
            std::vector<cl::tape_double> calculatedFunction_;
        };

        DefaultProbabilitiesTestData()
            : SwapData()

            , outPerform_(OUTPUT_FOLDER_NAME "//DefaultProbabilities"
            , { { "filename", "AdjointPerformance" }
              , { "not_clear", "Not" }
              , { "cleanlog", "true" }
              , { "title", "Credit default swap NPV and fair spread differentiation performance with respect to default probability" }
              , { "xlabel", "Number of default probabilities" }
              , { "ylabel", "Time (s)" }
              , { "smooth", "12" }
              , { "line_box_width", "-5" }})

            , outAdjoint_(OUTPUT_FOLDER_NAME "//DefaultProbabilities"
            , { { "filename", "Adjoint" }
              , { "not_clear", "Not" }
              , { "cleanlog", "false" }
              , { "title", "Credit default swap NPV and fair spread adjoint differentiation performance with respect to default probability" }
              , { "xlabel", "Number of default probabilities" }
              , { "ylabel", "Time (s)" }
              , { "smooth", "12" }
              })

            , outSize_(OUTPUT_FOLDER_NAME "//DefaultProbabilities"
            , { { "filename", "TapeSize" }
              , { "not_clear", "Not" }
              , { "cleanlog", "false" }
              , { "title", "Tape size dependence on number of default probabilities" }
              , { "xlabel", "Number of default probabilities" }
              , { "ylabel", "Size (MB)" } })

            , outNpv_(OUTPUT_FOLDER_NAME "//DefaultProbabilities//output"
            , { { "filename", "NPV on default probability" }
              , { "not_clear", "Not" }
              , { "title", "Credit default swap NPV dependence on default probability" }
              , { "xlabel", "Default probability" }
              , { "ylabel", "NPV" } })

            , outFairSpread_(OUTPUT_FOLDER_NAME "//DefaultProbabilities//output"
            , { { "filename", "Fair spread on default probability" }
              , { "not_clear", "Not" }
              , { "title", "Credit default swap fair spread dependence on default probability" }
              , { "xlabel", "Default probability" }
              , { "ylabel", "Fair spread" }
              })

#if defined CL_GRAPH_GEN
            , pointNo_(50)
            , iterNo_(50)
            , step_(12)
#else
            , pointNo_(1)
            , iterNo_(1)
            , step_(1)
#endif
        {
        }

        bool makeOutput()
        {
            bool ok = true;
            if (pointNo_ > 0)
            {
                ok &= recordDependencePlot();
            }
            ok &= cl::recordPerformance(*this, iterNo_, step_);
            return ok;
        }

        std::shared_ptr<Test> getTest(Size size)
        {
            return std::make_shared<Test>(size + 12, this);
        }

        // Makes plots for strike sensitivity dependence.
        bool recordDependencePlot()
        {
            std::vector<RealPair> outNpv(pointNo_);
            std::vector<RealPair> outFairSpread(pointNo_);
            auto test = getTest(0);
            Real dpDelta = (test->probability_[2] - test->probability_[0]) / (pointNo_ + 1);
            test->probability_[1] = test->probability_[0];

            for (Size i = 0; i < pointNo_; i++)
            {
                test->probability_[1] += dpDelta;
                std::vector<Real> out = test->data_->calculateDefaultProbabilitiesFunction(test->probability_);
                outNpv[i] = { test->probability_[1], out[0] };
                outFairSpread[i] = { test->probability_[1], out[1] };
            }
            outNpv_ << outNpv;
            outFairSpread_ << outFairSpread;
            return true;
        }

        Size pointNo_;
        Size iterNo_;
        Size step_;

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output outNpv_;
        cl::tape_empty_test_output outFairSpread_;
    };
    typedef DefaultProbabilitiesTestData::Test DefaultProbabilitiesTest;

    struct DiscountFactorTestData
        : public SwapData
    {
        struct Test
        : public cl::AdjointTest<Test>
        {
            static const CalcMethod default_method = other;

            Test(Size size, DiscountFactorTestData* data)
                : size_(size)
                , data_(data)
                , discountFactor_(size)
                , doubleDiscountFactor_(size)
                , iterNumFactor_(1000)
            {
                setLogger(&data_->outPerform_);

                double coef = -1.0 / size_;
                for (Size i = 0; i < size_; i++)
                    doubleDiscountFactor_[i] = (double)(std::exp(coef * (i + 1)));

                std::copy(doubleDiscountFactor_.begin(), doubleDiscountFactor_.end(), discountFactor_.begin());

                data->setDiscountFactorData();
            }

            Size indepVarNumber() { return size_; }

            Size depVarNumber() { return 2; }

            Size minPerfIteration() { return iterNumFactor_; }

            void recordTape()
            {
                cl::Independent(discountFactor_);
                calculateFunction();
                f_ = std::make_unique<cl::tape_function<double>>(discountFactor_, calculatedFunction_);
            }

            void calculateFunction()
            {
                calculatedFunction_ = data_->calculateDiscountFactorFunction(discountFactor_);
            }

            // Calculates derivatives using adjoint.
            void calcAdjoint()
            {
                adjointResults_ = f_->Jacobian(doubleDiscountFactor_);
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                double h = 1.0e-6;  // shift for finite diff. method
                analyticalResults_.resize(size_);
                data_->calculateFinDiff(discountFactor_
                                        , h
                                        , analyticalResults_
                                        , [this] (std::vector<Real>& v) -> std::vector<Real>
                {
                    return data_->calculateDiscountFactorFunction(v);
                });
            }

            double relativeTol() const { return 1e-4; }

            double absTol() const { return 1e-5; }

            Size size_;
            Size iterNumFactor_;
            DiscountFactorTestData* data_;
            std::vector<cl::tape_double> discountFactor_;
            std::vector<double> doubleDiscountFactor_;
            std::vector<cl::tape_double> calculatedFunction_;
        };

        DiscountFactorTestData()
            : SwapData()

            , outPerform_(OUTPUT_FOLDER_NAME "//DiscountFactor"
            , { { "filename", "AdjointPerformance" }
              , { "not_clear", "Not" }
              , { "cleanlog", "true"}
              , { "title", "Credit default swap NPV and fair spread differentiation performance with respect to discount factors" }
              , { "xlabel", "Number of discount factors" }
              , { "ylabel", "Time (s)" }
              , { "smooth", "12" }
              , { "line_box_width", "-5" } })

            , outAdjoint_(OUTPUT_FOLDER_NAME "//DiscountFactor"
            , { { "filename", "Adjoint" }
              , { "not_clear", "Not" }
              , { "cleanlog", "false" }
              , { "title", "Credit default swap NPV and fair spread adjoint differentiation performance with respect to discount factors" }
              , { "xlabel", "Number of discount factors" }
              , { "ylabel", "Time (s)" }
              , { "smooth", "12" } })

            , outSize_(OUTPUT_FOLDER_NAME "//DiscountFactor"
            , { { "filename", "TapeSize" }
              , { "not_clear", "Not" }
              , { "cleanlog", "false" }
              , { "title", "Tape size dependence on number of discount factors" }
              , { "xlabel", "Number of discount factors" }
              , { "ylabel", "Size (MB)" } })

            , outNpv_(OUTPUT_FOLDER_NAME "//DiscountFactor//output"
            , { { "filename", "NPV on discount factor" }
              , { "not_clear", "Not" }
              , { "cleanlog", "false" }
              , { "title", "Credit default swap NPV dependence on discount factor" }
              , { "xlabel", "Discount factor" }
              , { "ylabel", "NPV" } })

           , outFairSpread_(OUTPUT_FOLDER_NAME "//DiscountFactor//output"
           , { { "filename", "Fair spread on discount factor" }
             , { "not_clear", "Not" }
             , { "cleanlog", "false" }
             , { "title", "Credit default swap fair spread dependence on discount factor" }
             , { "xlabel", "Discount factor" }
             , { "ylabel", "Fair spread" } })

#if defined CL_GRAPH_GEN
            , pointNo_(50)
            , iterNo_(50)
            , step_(12)
#else
            , pointNo_(1)
            , iterNo_(1)
            , step_(1)
#endif
        {
        }

        bool makeOutput()
        {
            bool ok = true;
            if (pointNo_ > 0)
            {
                ok &= recordDependencePlot();
            }
            ok &= cl::recordPerformance(*this, iterNo_, step_);
            return ok;
        }

        std::shared_ptr<Test> getTest(Size size)
        {
            return std::make_shared<Test>(size + 12, this);
        }

        // Makes plots for strike sensitivity dependence.
        bool recordDependencePlot()
        {
            std::vector<RealPair> outNpv(pointNo_);
            std::vector<RealPair> outFairSpread(pointNo_);
            auto test = getTest(0);

            Real dfsDelta = (test->discountFactor_[0] - test->discountFactor_[2]) / (pointNo_ + 1);
            test->discountFactor_[1] = test->discountFactor_[2];
            for (Size i = 0; i < pointNo_; i++)
            {
                test->discountFactor_[1] += dfsDelta;
                std::vector<Real> out = test->data_->calculateDiscountFactorFunction(test->discountFactor_);
                outNpv[i] = { test->discountFactor_[1], out[0] };
                outFairSpread[i] = { test->discountFactor_[1], out[1] };
            }
            outNpv_ << outNpv;
            outFairSpread_ << outFairSpread;
            return true;
        }

        Size pointNo_;
        Size iterNo_;
        Size step_;

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output outNpv_;
        cl::tape_empty_test_output outFairSpread_;
    };
    typedef DiscountFactorTestData::Test DiscountFactorTest;

    struct NotionalSpreadRateTest
        : public cl::AdjointTest<NotionalSpreadRateTest>
    {
        static const CalcMethod default_method = other;

        NotionalSpreadRateTest()
            : iterNumFactor_(1)
            , parameters_(3)
        {
            doubleParameters_ = { 10000.0, 0.012, 0.4 };
            std::copy(doubleParameters_.begin(), doubleParameters_.end(), parameters_.begin());
        }

        Size indepVarNumber() { return 3; }

        Size depVarNumber() { return 3; }

        Size minPerfIteration() { return iterNumFactor_; }

        void recordTape()
        {
            cl::Independent(parameters_);
            calculateFunction();
            f_ = std::make_unique<cl::tape_function<double>>(parameters_, calculatedFunction_);
        }

        void calculateFunction()
        {
            calculatedFunction_ = swapData_.calculateNotionalSpreadRateFunction(parameters_);
        }

        // Calculates derivatives using adjoint.
        void calcAdjoint()
        {
            adjointResults_ = f_->Jacobian(doubleParameters_);
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            double h = 1e-4;
            analyticalResults_.resize(indepVarNumber() * depVarNumber());
            swapData_.calculateFinDiff(parameters_
                                , h
                                , analyticalResults_
                                , [this] (std::vector<Real>& v) -> std::vector<Real>
            {
                return swapData_.calculateNotionalSpreadRateFunction(v);
            });

        }

        double relativeTol() const { return 1e-4; }

        double absTol() const { return 1e-5; }

        Size size_;
        Size iterNumFactor_;
        std::vector<cl::tape_double> parameters_;
        std::vector<double> doubleParameters_;
        std::vector<cl::tape_double> calculatedFunction_;
        SwapData swapData_;
    };
}

#endif