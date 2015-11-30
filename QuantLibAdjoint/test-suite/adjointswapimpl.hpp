/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2003, 2004, 2007 StatPro Italia srl
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

// Based on swap.cpp file from test-suite.

#ifndef cl_adjoint_swap_impl_hpp
#define cl_adjoint_swap_impl_hpp
#pragma once

#include <ql/quantlib.hpp>
#include "utilities.hpp"
#include "adjointswaptest.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include <boost/shared_ptr.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;

#define OUTPUT_FOLDER_NAME "AdjointSwap"

namespace
{
    struct Variation
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                " ", ""
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, Variation& v)
        {
                stm << v.param_
                    << ";" << v.swapNpv_
                    << std::endl;
                return stm;
            }

        Real param_;
        Real swapNpv_;
    };

    struct SwapData
    {
        SwapData()
        : type_(VanillaSwap::Payer)
        , settlementDays_(2)
        , nominal_(100.0)
        , fixedConvention_(Unadjusted)
        , floatingConvention_(ModifiedFollowing)
        , fixedFrequency_(Annual)
        , floatingFrequency_(Semiannual)
        , fixedDayCount_(Thirty360())
        , lengths_({ 1, 2, 5, 10, 20 })
        , rates_({ 0.04, 0.05, 0.06, 0.07 })
        , spreads_({ -0.01, -0.001, 0.001, 0.01 })
        , backup_()
        {
            index_ = boost::make_shared<Euribor>(Period(floatingFrequency_), termStructure_);
            calendar_ = index_->fixingCalendar();
            today_ = calendar_.adjust(Date::todaysDate());
            Settings::instance().evaluationDate() = today_;
            settlement_ = calendar_.advance(today_, settlementDays_, Days);
            termStructure_.linkTo(flatRate(settlement_, 0.05, Actual365Fixed()));
        }

        boost::shared_ptr<VanillaSwap> makeSwap(Integer length
                                                , Rate fixedRate
                                                , Spread floatingSpread)
        {
            Date maturity = calendar_.advance(settlement_
                                              , length
                                              , Years
                                              , floatingConvention_);

            Schedule fixedSchedule(settlement_
                                   , maturity
                                   , Period(fixedFrequency_)
                                   , calendar_
                                   , fixedConvention_
                                   , fixedConvention_
                                   , DateGeneration::Forward
                                   , false);

            Schedule floatSchedule(settlement_
                                   , maturity
                                   , Period(floatingFrequency_)
                                   , calendar_
                                   , floatingConvention_
                                   , floatingConvention_
                                   , DateGeneration::Forward
                                   , false);

            boost::shared_ptr<VanillaSwap> swap =
                boost::make_shared<VanillaSwap>(type_
                                , nominal_
                                , fixedSchedule
                                , fixedRate
                                , fixedDayCount_
                                , floatSchedule
                                , index_
                                , floatingSpread
                                , index_->dayCounter());

            swap->setPricingEngine(boost::make_shared<DiscountingSwapEngine>(termStructure_));

            return swap;
        }


        Real calculateRateDependencyNPV(Real rate)
        {
            Real Y = 0;
            for (Size i = 0; i < lengths_.size(); i++)
            {
                for (Size j = 0; j < spreads_.size(); j++)
                {
                    boost::shared_ptr<VanillaSwap> swap =
                        makeSwap(lengths_[i], rate, spreads_[j]);
                    Y += swap->NPV();
                }
            }
            return Y;
        }

        Real calculateSpreadDependencyNPV(Real spread)
        {
            Real Y = 0;

            for (Size i = 0; i < lengths_.size(); i++)
            {
                for (Size j = 0; j < rates_.size(); j++)
                {
                    boost::shared_ptr<VanillaSwap> swap =
                        makeSwap(lengths_[i], rates_[j], spread);
                    Y += swap->NPV();
                }
            }
            return Y;
        }

        Real calculateInArrearsNPV(Real capletVolatility)
        {
            Date maturity = today_ + 5 * Years;
            Calendar calendar = NullCalendar();

            Schedule schedule(today_
                              , maturity
                              , Period(Annual)
                              , calendar
                              , Following
                              , Following
                              , DateGeneration::Forward
                              , false);

            DayCounter dayCounter = SimpleDayCounter();
            std::vector<Real> nominals(1, 100000000.0);

            boost::shared_ptr<IborIndex> index_
                = boost::make_shared<IborIndex>("dummy"
                , 1 * Years
                , 0
                , EURCurrency()
                , calendar
                , Following
                , false
                , dayCounter
                , termStructure_);
            Rate oneYear = 0.05;
            Rate r = std::log(1.0 + oneYear);
            termStructure_.linkTo(flatRate(today_, r, dayCounter));


            std::vector<Rate> coupons(1, oneYear);
            Leg fixedLeg = FixedRateLeg(schedule)
                .withNotionals(nominals)
                .withCouponRates(coupons, dayCounter);

            std::vector<Real> gearings;
            std::vector<Rate> spreads;
            Natural fixingDays = 0;

            Real Y = 0;

            Handle<OptionletVolatilityStructure> vol(
                boost::make_shared<ConstantOptionletVolatility>(today_
                , NullCalendar()
                , Following
                , capletVolatility
                , dayCounter));

            boost::shared_ptr<IborCouponPricer> pricer
                = boost::make_shared<BlackIborCouponPricer>(vol);

            Leg floatingLeg = IborLeg(schedule, index_)
                .withNotionals(nominals)
                .withPaymentDayCounter(dayCounter)
                .withFixingDays(fixingDays)
                .withGearings(gearings)
                .withSpreads(spreads)
                .inArrears();
            setCouponPricer(floatingLeg, pricer);

            Swap swap(floatingLeg, fixedLeg);
            swap.setPricingEngine(boost::make_shared<DiscountingSwapEngine>(termStructure_));

            return swap.NPV();
        }

        template <class Func>
        void calculateFinDiff(std::vector<Real>& X
                              , Real h
                              , std::vector<Real>& sf_Finite
                              , Func swapNpv)
        {
            Size sizeX = X.size();

            std::vector<Real> swapPortfolioNpv(sizeX);
            for (int i = 0; i < sizeX; i++)
            {
                swapPortfolioNpv[i] = swapNpv(X[i]);
            }

            sf_Finite.resize(sizeX);

            for (int i = 0; i < sizeX; i++)
            {
                sf_Finite[i] = (swapNpv(X[i] + h) - swapPortfolioNpv[i]) / h;
            }
        }

        // Global data.
        Calendar calendar_;
        Date today_, settlement_;
        VanillaSwap::Type type_;
        Real nominal_;
        BusinessDayConvention fixedConvention_, floatingConvention_;
        Frequency fixedFrequency_, floatingFrequency_;
        DayCounter fixedDayCount_;
        Natural settlementDays_;
        RelinkableHandle<YieldTermStructure> termStructure_;
        boost::shared_ptr<IborIndex> index_;

        std::vector<Integer> lengths_;
        std::vector<Rate> rates_;
        std::vector<Spread> spreads_;

        // Cleanup.
        SavedSettings backup_;
    };

    struct RateDependencyTestData
        : public SwapData
    {
        struct Test
        : public cl::AdjointTest<Test>
        {
            explicit Test(Size size, RateDependencyTestData* data)
            : size_(size)
            , data_(data)
            , rate_(size)
            , iterNumFactor_(1000)
            {
                setLogger(&data_->outPerform_);

                Rate startRate = 0.01;
                Rate maxRate = 0.06;
                Real step = (maxRate - startRate) / size_;

                for (Size i = 0; i < size_; i++)
                {
                    rate_[i] = startRate + i*step;
                }

            }

            Size indepVarNumber() { return size_; }

            Size depVarNumber() { return 1; }

            Size minPerfIteration() { return iterNumFactor_; }

            void recordTape()
            {
                cl::Independent(rate_);
                calculateTotalSwapNpv();
                f_ = std::make_unique<cl::tape_function<double>>(rate_, swapNpv_);
            }

            // Calculates total calibration error.
            void calculateTotalSwapNpv()
            {
                swapNpv_.resize(depVarNumber());
                for (Size i = 0; i < size_; i++)
                    swapNpv_[0] += data_->calculateRateDependencyNPV(rate_[i]);
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                double h = 1.0e-12;  // shift for finite diff. method
                analyticalResults_.resize(size_);
                data_->calculateFinDiff(rate_
                                        , h
                                        , analyticalResults_
                                        , [this] (Real v) -> Real
                                          {
                                              return data_->calculateRateDependencyNPV(v);
                                          });

            }

            double relativeTol() const { return 1e-2; }

            double absTol() const { return 1e-10; }

            Size size_;
            Size iterNumFactor_;
            RateDependencyTestData* data_;
            std::vector<cl::tape_double> rate_;
            std::vector<cl::tape_double> swapNpv_;
        };

        RateDependencyTestData()
            : SwapData()

            , outPerform_(OUTPUT_FOLDER_NAME "//RateDependency"
            , { { "filename", "AdjointPerformance" }
              , { "not_clear", "Not" }
              , { "line_box_width", "-5" }
              , { "title", "Swap NPV differentiation performance with respect to rate" }
              , { "ylabel", "Time (s)" }
              , { "xlabel", "Number of rates" }
              , { "smooth", "default" }
              , { "cleanlog", "true" } })

            , outAdjoint_(OUTPUT_FOLDER_NAME "//RateDependency"
            , { { "filename", "Adjoint" }
              , { "not_clear", "Not" }
              , { "smooth", "default" }
              , { "title", "Swap NPV adjoint differentiation performance with respect to rate" }
              , { "cleanlog", "false" }
              , { "ylabel", "Time (s)" }
              , { "xlabel", "Number of rates" } })

            , outSize_(OUTPUT_FOLDER_NAME "//RateDependency"
            , { { "filename", "TapeSize" }
              , { "not_clear", "Not" }
              , { "title", "Tape size dependence on number of rates" }
              , { "cleanlog", "false" }
              , { "smooth", "default" }
              , { "ylabel", "Size (MB)" }
              , { "xlabel", "Number of rates" } })

            , out_(OUTPUT_FOLDER_NAME "//RateDependency//output"
            , { { "filename", "SwapNPVonRates" }
              , { "ylabel", "Swap NPV" }
              , { "not_clear", "Not" }
              , { "title", "Swap NPV dependence on rate" }
              , { "cleanlog", "false" }
              , { "xlabel", "Rate" } })

#if defined CL_GRAPH_GEN
            , pointNo_(10)
            , iterNo_(10)
            , step_(1)
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

        std::shared_ptr<Test> getTest(size_t size)
        {
            return std::make_shared<Test>(size, this);
        }

        // Makes plots for rate sensitivity dependence.
        bool recordDependencePlot()
        {
            std::vector<Variation> outData(pointNo_);
            auto test = getTest(pointNo_);
            for (Size i = 0; i < pointNo_; i++)
            {
                outData[i].param_ = test->rate_[i];
                outData[i].swapNpv_ = test->data_->calculateRateDependencyNPV(test->rate_[i]);
            }
            out_ << outData;
            return true;
        }

        Size pointNo_;
        Size iterNo_;
        Size step_;

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output out_;
    };
    typedef RateDependencyTestData::Test RateDependencyTest;

    struct SpreadDependencyTestData
        : public SwapData
    {
        struct Test
        : public cl::AdjointTest<Test>
        {
            explicit Test(Size size, SpreadDependencyTestData* data)
            : size_(size)
            , data_(data)
            , spread_(size)
            , iterNumFactor_(1000)
            {
                setLogger(&data_->outPerform_);

                Spread startSpread = -0.01;
                Spread maxSpread = 0.01;
                Real step = (maxSpread - startSpread) / size_;

                for (Size i = 0; i < size_; i++)
                {
                    spread_[i] = (startSpread + i*step == 0) ? 0.00001 : startSpread + i*step;
                }

            }

            Size indepVarNumber() { return size_; }

            Size depVarNumber() { return 1; }

            Size minPerfIteration() { return iterNumFactor_; }

            void recordTape()
            {
                cl::Independent(spread_);
                calculateTotalSwapNpv();
                f_ = std::make_unique<cl::tape_function<double>>(spread_, swapNpv_);
            }

            // Calculates total calibration error.
            void calculateTotalSwapNpv()
            {
                swapNpv_.resize(depVarNumber());
                for (Size i = 0; i < size_; i++)
                    swapNpv_[0] += data_->calculateSpreadDependencyNPV(spread_[i]);
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                double h = 1.0e-12;  // shift for finite diff. method
                analyticalResults_.resize(size_);
                data_->calculateFinDiff(spread_
                                        , h
                                        , analyticalResults_
                                        , [this] (Real v) -> Real
                                          {
                                              return data_->calculateSpreadDependencyNPV(v);
                                          });
            }

            double relativeTol() const { return 1e-2; }

            double absTol() const { return 1e-10; }

            Size size_;
            Size iterNumFactor_;
            SpreadDependencyTestData* data_;
            std::vector<cl::tape_double> spread_;
            std::vector<cl::tape_double> swapNpv_;
        };

        SpreadDependencyTestData()
            : SwapData()

            , outPerform_(OUTPUT_FOLDER_NAME "//SpreadDependency"
            , { { "filename", "AdjointPerformance" }
              , { "not_clear", "Not" }
              , { "line_box_width", "-5" }
              , { "title", "Swap NPV differentiation performance with respect to spread" }
              , { "ylabel", "Time (s)" }
              , { "xlabel", "Number of spreads" }
              , { "smooth", "default" }
              , { "cleanlog", "true" } })

            , outAdjoint_(OUTPUT_FOLDER_NAME "//SpreadDependency"
            , { { "filename", "Adjoint" }
              , { "not_clear", "Not" }
              , { "smooth", "default" }
              , { "title", "Swap NPV adjoint differentiation performance with respect to spread" }
              , { "cleanlog", "false" }
              , { "ylabel", "Time (s)" }
              , { "xlabel", "Number of spreads" } })

            , outSize_(OUTPUT_FOLDER_NAME "//SpreadDependency"
            , { { "filename", "TapeSize" }
              , { "not_clear", "Not" }
              , { "title", "Tape size dependence on number of spreads" }
              , { "cleanlog", "false" }
              , { "smooth", "default" }
              , { "ylabel", "Size (MB)" }
              , { "xlabel", "Number of spreads" } })

            , out_(OUTPUT_FOLDER_NAME "//SpreadDependency//output"
            , { { "filename", "SwapNPVonSpreads" }
              , { "ylabel", "Swap NPV" }
              , { "not_clear", "Not" }
              , { "title", "Swap NPV dependence on spread" }
              , { "cleanlog", "false" }
              , { "xlabel", "Spread" } })

#if defined CL_GRAPH_GEN
            , pointNo_(10)
            , iterNo_(10)
            , step_(1)
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

        std::shared_ptr<Test> getTest(size_t size)
        {
            return std::make_shared<Test>(size, this);
        }

        // Makes plots for spread sensitivity dependence.
        bool recordDependencePlot()
        {
            std::vector<Variation> outData(pointNo_);
            auto test = getTest(pointNo_);
            for (Size i = 0; i < pointNo_; i++)
            {
                outData[i].param_ = test->spread_[i];
                outData[i].swapNpv_ = test->data_->calculateSpreadDependencyNPV(test->spread_[i]);
            }
            out_ << outData;
            return true;
        }

        Size pointNo_;
        Size iterNo_;
        Size step_;

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output out_;
    };
    typedef SpreadDependencyTestData::Test SpreadDependencyTest;

    struct InArrearsTestData
        : public SwapData
    {
        struct Test
        : public cl::AdjointTest<Test>
        {
            explicit Test(Size size, InArrearsTestData* data)
            : size_(size)
            , data_(data)
            , volatility_(size)
            , iterNumFactor_(1000)
            {
                setLogger(&data_->outPerform_);

                Volatility startVol = 0.18;
                Volatility maxVol = 0.22;
                Real step = (maxVol - startVol) / size_;

                for (Size i = 0; i < size_; i++)
                {
                    volatility_[i] = startVol + i * step;
                }
            }

            Size indepVarNumber() { return size_; }

            Size depVarNumber() { return 1; }

            Size minPerfIteration() { return iterNumFactor_; }

            void recordTape()
            {
                cl::Independent(volatility_);
                calculateTotalSwapNpv();
                f_ = std::make_unique<cl::tape_function<double>>(volatility_, swapNpv_);
            }

            // Calculates total calibration error.
            void calculateTotalSwapNpv()
            {
                swapNpv_.resize(depVarNumber());
                for (Size i = 0; i < size_; i++)
                    swapNpv_[0] += data_->calculateInArrearsNPV(volatility_[i]);
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                double h = 1.0e-12;  // shift for finite diff. method
                analyticalResults_.resize(size_);
                data_->calculateFinDiff(volatility_
                                        , h
                                        , analyticalResults_
                                        , [this] (Real v) -> Real
                                          {
                                              return data_->calculateInArrearsNPV(v);
                                          });

            }

            double relativeTol() const { return 1e-2; }

            double absTol() const { return 1e-10; }

            Size size_;
            Size iterNumFactor_;
            InArrearsTestData* data_;
            std::vector<cl::tape_double> volatility_;
            std::vector<cl::tape_double> swapNpv_;
        };

        InArrearsTestData()
            : SwapData()

            , outPerform_(OUTPUT_FOLDER_NAME "//InArrears"
            , { { "filename", "AdjointPerformance" }
              , { "not_clear", "Not" }
              , { "line_box_width", "-5" }
              , { "title", "Swap NPV differentiation performance with respect to volatility" }
              , { "ylabel", "Time (s)" }
              , { "xlabel", "Number of volatilities" }
              , { "smooth", "default" }
              , { "cleanlog", "true" } })

            , outAdjoint_(OUTPUT_FOLDER_NAME "//InArrears"
            , { { "filename", "Adjoint" }
              , { "not_clear", "Not" }
              , { "smooth", "default" }
              , { "title", "Swap NPV adjoint differentiation performance with respect to volatility" }
              , { "cleanlog", "false" }
              , { "ylabel", "Time (s)" }
              , { "xlabel", "Number of volatilities" } })

            , outSize_(OUTPUT_FOLDER_NAME "//InArrears"
            , { { "filename", "TapeSize" }
              , { "not_clear", "Not" }
              , { "title", "Tape size dependence on number of volatilities" }
              , { "cleanlog", "false" }
              , { "smooth", "default" }
              , { "ylabel", "Size (MB)" }
              , { "xlabel", "Number of volatilities" } })

            , out_(OUTPUT_FOLDER_NAME "//InArrears//output"
            , { { "filename", "SwapNPVonVolatilities" }
              , { "ylabel", "Swap NPV" }
              , { "not_clear", "Not" }
              , { "title", "Swap NPV dependence on volatility" }
              , { "cleanlog", "false" }
              , { "xlabel", "Spread" } })

#if defined CL_GRAPH_GEN
            , pointNo_(20)
            , iterNo_(10)
            , step_(20)
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

        std::shared_ptr<Test> getTest(size_t size)
        {
            return std::make_shared<Test>(size, this);
        }

        // Makes plots for volatility sensitivity dependence.
        bool recordDependencePlot()
        {
            std::vector<Variation> outData(pointNo_);
            auto test = getTest(pointNo_);
            for (Size i = 0; i < pointNo_; i++)
            {
                outData[i].param_ = test->volatility_[i];
                outData[i].swapNpv_ = test->data_->calculateInArrearsNPV(test->volatility_[i]);
            }
            out_ << outData;
            return true;
        }

        Size pointNo_;
        Size iterNo_;
        Size step_;

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output out_;
    };
    typedef InArrearsTestData::Test InArrearsTest;
}

#endif