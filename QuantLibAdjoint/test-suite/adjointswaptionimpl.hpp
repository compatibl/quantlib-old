/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2003 RiskMap srl
Copyright (C) 2006, 2007 Ferdinando Ametrano
Copyright (C) 2006 Marco Bianchetti
Copyright (C) 2006 Cristina Duminuco
Copyright (C) 2007, 2008 StatPro Italia srl
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

// based on swaption.cpp file from test-suite

#ifndef cl_adjoint_swaption_impl_hpp
#define cl_adjoint_swaption_impl_hpp
#pragma once

#include <ql/quantlib.hpp>
#include "utilities.hpp"
#include "adjointswaptiontest.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include <boost/shared_ptr.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;


#define OUTPUT_FOLDER_NAME "AdjointSwaption"

namespace
{
    struct Variation
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "param", ""
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type& operator << (stream_type& stm, Variation& v)
        {
            stm << v.param_
                << ";" << v.swaptionNpv_
                << std::endl;
            return stm;
        }

        Real param_;
        Real swaptionNpv_;
    };

    struct SwaptionData
    {
        SwaptionData()
        : settlementDays_(2)
        , nominal_(1000000.0)
        , fixedConvention_(Unadjusted)
        , fixedFrequency_(Annual)
        , fixedDayCount_(Thirty360())
        , type_(VanillaSwap::Payer)
        {
            index_ = boost::make_shared<Euribor6M>(termStructure_);
            floatingTenor_ = index_->tenor();
            floatingConvention_ = index_->businessDayConvention();
            calendar_ = index_->fixingCalendar();
            today_ = calendar_.adjust(Date::todaysDate());
            Settings::instance().evaluationDate() = today_;
            settlement_ = calendar_.advance(today_, settlementDays_, Days);
            termStructure_.linkTo(flatRate(settlement_, 0.05, Actual365Fixed()));

            exercises_ = { 1 * Years, 3 * Years, 7 * Years, 10 * Years };
            lengths_ = { 1 * Years, 2 * Years, 5 * Years, 10 * Years, 20 * Years };
        }

        // Make swaption with given parameters.
        boost::shared_ptr<Swaption> makeSwaption(const boost::shared_ptr<VanillaSwap>& swap
                                                 , const Date& exercise
                                                 , Volatility volatility
                                                 , Settlement::Type settlementType = Settlement::Physical)
        {
            // Set volatility from simple quote.
            Handle<Quote> vol(boost::make_shared<SimpleQuote>(volatility));

            // Create pricing engine.
            boost::shared_ptr<PricingEngine> engine(new BlackSwaptionEngine(termStructure_, vol));

            // Create swaption using vanilla swap.
            boost::shared_ptr<Swaption> result(new Swaption(swap
                , boost::make_shared<EuropeanExercise>(exercise)
                , settlementType));

            // Set pricing engine in swaption.
            result->setPricingEngine(engine);
            return result;
        }


        // Calculate sum of swaption NPV for given spreads.
        std::vector<Real> swaptionNpvOnSpreads(Real spread)
        {
            std::vector<Real> Y(1);
            for (Size i = 0; i < exercises_.size(); i++)
            {
                for (Size j = 0; j < lengths_.size(); j++)
                {
                    // Set date of exercize.
                    Date exerciseDate = calendar_.advance(today_, exercises_[i]);

                    // Set start date.
                    Date startDate = calendar_.advance(exerciseDate, settlementDays_, Days);

                    // Create vanilla swap.
                    boost::shared_ptr<VanillaSwap> swap =
                        MakeVanillaSwap(lengths_[j], index_, 0.06)
                        .withFixedLegTenor(1 * Years)
                        .withFixedLegDayCount(fixedDayCount_)
                        .withEffectiveDate(startDate)
                        .withFloatingLegSpread(spread)
                        .withType(type_);

                    // Create swaption using vanilla swap.
                    boost::shared_ptr<Swaption> swaption =
                        makeSwaption(swap, exerciseDate, 0.20);

                    Y[0] += swaption->NPV();
                }
            }
            return Y;
        }

        // Calculate sum of swaption NPV for given spreads.
        std::vector<Real> swaptionNpvWithSpreadCorrection(Real spread)
        {
            std::vector<Real> Y(4, 0);

            for (Size i = 0; i < exercises_.size(); i++)
            {
                for (Size j = 0; j < lengths_.size(); j++)
                {
                    // Set date of exercize.
                    Date exerciseDate = calendar_.advance(today_, exercises_[i]);

                    // Set start date.
                    Date startDate = calendar_.advance(exerciseDate, settlementDays_, Days);

                    // Create vanilla swap.
                    boost::shared_ptr<VanillaSwap> swap =
                        MakeVanillaSwap(lengths_[j], index_, 0.06)
                        .withFixedLegTenor(1 * Years)
                        .withFixedLegDayCount(fixedDayCount_)
                        .withEffectiveDate(startDate)
                        .withFloatingLegSpread(spread)
                        .withType(type_);

                    // Calculate spread correction.
                    Spread correction = spread *
                        swap->floatingLegBPS() /
                        swap->fixedLegBPS();

                    // Create equivalent vanilla swap.
                    boost::shared_ptr<VanillaSwap> equivalentSwap =
                        MakeVanillaSwap(lengths_[j], index_, 0.06 + correction)
                        .withFixedLegTenor(1 * Years)
                        .withFixedLegDayCount(fixedDayCount_)
                        .withEffectiveDate(startDate)
                        .withFloatingLegSpread(0.0)
                        .withType(type_);

                    // Create swaption using vanilla swap.
                    boost::shared_ptr<Swaption> swaption1 =
                        makeSwaption(swap, exerciseDate, 0.20);

                    Y[0] += swaption1->NPV();

                    // Create swaption using equivalent vanilla swap.
                    boost::shared_ptr<Swaption> swaption2 =
                        makeSwaption(equivalentSwap, exerciseDate, 0.20);

                    Y[1] += swaption2->NPV();

                    // Create swaption with cash settlement using vanilla swap.
                    boost::shared_ptr<Swaption> swaption1_cash =
                        makeSwaption(swap, exerciseDate, 0.20,
                        Settlement::Cash);

                    Y[2] += swaption1_cash->NPV();

                    // Create swaption with cash settlement using equivalent vanilla swap.
                    boost::shared_ptr<Swaption> swaption2_cash =
                        makeSwaption(equivalentSwap, exerciseDate, 0.20,
                        Settlement::Cash);

                    Y[3] += swaption2_cash->NPV();

                    // Check swaption NPV and equivalent swaption NPV
                    if (std::fabs(swaption1->NPV() - swaption2->NPV()) > 1.0e-6)
                        BOOST_ERROR("wrong spread treatment:" <<
                        "\nexercise: " << exerciseDate <<
                        "\nlength:   " << lengths_[j] <<
                        "\ntype      " << type_ <<
                        "\nspread:   " << io::rate(spread) <<
                        "\noriginal swaption value:   " << swaption1->NPV() <<
                        "\nequivalent swaption value: " << swaption2->NPV());

                    if (std::fabs(swaption1_cash->NPV() - swaption2_cash->NPV()) > 1.0e-6)
                        BOOST_ERROR("wrong spread treatment:" <<
                        "\nexercise date: " << exerciseDate <<
                        "\nlength: " << lengths_[j] <<
                        "\npay " << (type_ ? "fixed" : "floating") <<
                        "\nspread: " << io::rate(spread) <<
                        "\nvalue of original swaption:   " << swaption1_cash->NPV() <<
                        "\nvalue of equivalent swaption: " << swaption2_cash->NPV());

                }
            }
            return Y;
        }

        // Calculate sum of swaption NPV for given volatilities.
        std::vector<Real> swaptionNpvOnCachedValues(Real volatility)
        {
            today_ = Date(13, March, 2002);
            settlement_ = Date(15, March, 2002);
            Settings::instance().evaluationDate() = today_;
            termStructure_.linkTo(flatRate(settlement_, 0.05, Actual365Fixed()));

            // Set date of exercize.
            Date exerciseDate = calendar_.advance(settlement_, 5 * Years);

            // Set start date.
            Date startDate = calendar_.advance(exerciseDate,
                                               settlementDays_, Days);

            // Create vanilla swap.
            boost::shared_ptr<VanillaSwap> swap =
                MakeVanillaSwap(10 * Years, index_, 0.06)
                .withEffectiveDate(startDate)
                .withFixedLegTenor(1 * Years)
                .withFixedLegDayCount(fixedDayCount_);

            std::vector<Real> Y(1);
            // Create swaption using vanilla swap.
            boost::shared_ptr<Swaption> swaption =
                makeSwaption(swap, exerciseDate, volatility);
            Y[0] += swaption->NPV();

            return Y;
        }

        template <class Func>
        void calculateFinDiff(std::vector<Real>& X
                              , Size dep_number
                              , Real h
                              , std::vector<Real>& sf_Finite
                              , Func swaptionNpv)
        {
            Size sizeX = X.size();

            std::vector<Real> swaptionPortfolioNpv(sizeX * dep_number);
            for (int i = 0; i < sizeX; i++)
            {
                std::vector<Real> temp = swaptionNpv(X[i]);
                for (int j = 0; j < dep_number; j++)
                    swaptionPortfolioNpv[sizeX*j + i] = temp[j];
            }

            sf_Finite.resize(swaptionPortfolioNpv.size());

            for (int i = 0; i < sizeX; i++)
            {
                std::vector<Real> temp = swaptionNpv(X[i] + h);
                for (int j = 0; j < dep_number; j++)
                {
                    int curIndex = sizeX*j + i;
                    sf_Finite[curIndex] = (temp[j] - swaptionPortfolioNpv[curIndex]) / h;
                }
            }
        }

        // Global data.
        Date today_;
        Date settlement_;
        Real nominal_;
        Calendar calendar_;

        BusinessDayConvention fixedConvention_;
        Frequency fixedFrequency_;
        DayCounter fixedDayCount_;

        BusinessDayConvention floatingConvention_;
        Period floatingTenor_;
        boost::shared_ptr<IborIndex> index_;

        Natural settlementDays_;
        RelinkableHandle<YieldTermStructure> termStructure_;

        std::vector<Period> exercises_;
        std::vector<Period> lengths_;
        VanillaSwap::Type type_;

        // Cleanup
        SavedSettings backup_;

    };

    struct SpreadDependencyTestData
        : public SwaptionData
    {
        struct Test
        : public cl::AdjointTest<Test>
        {
            Test(Size size, SpreadDependencyTestData* data)
            : size_(size)
            , data_(data)
            , swaptionNpv_(1)
            , spread_(size)
            , iterNumFactor_(1)
            {
                setLogger(&data_->outPerform_);

                Real startSpread = -0.01;
                Real maxSpread = 0.01;
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
                calculateTotalSwaptionNpv();
                f_ = std::make_unique<cl::tape_function<double>>(spread_, swaptionNpv_);
            }

            // Calculates total calibration error.
            void calculateTotalSwaptionNpv()
            {
                for (int i = 0; i < size_; i++)
                    swaptionNpv_[0] += data_->swaptionNpvOnSpreads(spread_[i])[0];
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                double h = 1.0e-10;  // shift for finite diff. method
                analyticalResults_.resize(size_);
                data_->calculateFinDiff(spread_
                                        , depVarNumber()
                                        , h
                                        , analyticalResults_
                                        , [this] (Real v) -> std::vector<Real>
                {
                    return data_->swaptionNpvOnSpreads(v);
                });
            }

            double relativeTol() const { return 1e-2; }

            double absTol() const { return 1e-10; }

            Size size_;
            Size iterNumFactor_;
            SpreadDependencyTestData* data_;
            std::vector<cl::tape_double> spread_;
            std::vector<cl::tape_double> swaptionNpv_;
        };

        SpreadDependencyTestData()
            : SwaptionData()

            , outPerform_(OUTPUT_FOLDER_NAME "//SpreadDependency"
            , { { "filename", "AdjointPerformance" }
        , { "not_clear", "Not" }
        , { "line_box_width", "-5" }
        , { "title", "Swaption NPV differentiation performance with respect to spread" }
        , { "ylabel", "Time (s)" }
        , { "xlabel", "Number of spreads" }
        , { "smooth", "default" }
        , { "cleanlog", "true" } })

            , outAdjoint_(OUTPUT_FOLDER_NAME "//SpreadDependency"
            , { { "filename", "Adjoint" }
        , { "not_clear", "Not" }
        , { "smooth", "default" }
        , { "title", "Swaption NPV adjoint differentiation performance with respect to spread" }
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
            , { { "filename", "SwaptNPVonSpreads" }
        , { "ylabel", "Swaption NPV" }
        , { "not_clear", "Not" }
        , { "title", "Swaption NPV dependence on spread" }
        , { "cleanlog", "false" }
        , { "smooth", "default" }
        , { "xlabel", "Spread" } })

#if defined CL_GRAPH_GEN
            , pointNo_(50)
            , iterNo_(10)
            , step_(5)
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

        // Makes plots for strike sensitivity dependence.
        bool recordDependencePlot()
        {
            std::vector<Variation> outData(pointNo_);
            auto test = getTest(pointNo_);
            for (Size i = 0; i < pointNo_; i++)
            {
                outData[i].param_ = test->spread_[i];
                outData[i].swaptionNpv_ = test->data_->swaptionNpvOnSpreads(test->spread_[i])[0];
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

    struct SpreadCorrectionTestData
        : public SwaptionData
    {
        struct Test
        : public cl::AdjointTest<Test>
        {
            static const CalcMethod default_method = other;

            Test(Size size, SpreadCorrectionTestData* data)
                : size_(size)
                , data_(data)
                , spread_(size)
                , doubleSpread_(size)
                , iterNumFactor_(1)
            {
                setLogger(&data_->outPerform_);

                double startSpread = -0.01;
                double maxSpread = 0.01;
                double step = (maxSpread - startSpread) / size_;

                for (Size i = 0; i < size_; i++)
                {
                    doubleSpread_[i] = (startSpread + i*step == 0) ? 0.00001 : startSpread + i*step;
                }

                std::copy(doubleSpread_.begin(), doubleSpread_.end(), spread_.begin());
            }

            Size indepVarNumber() { return size_; }

            Size depVarNumber() { return 4; }

            Size minPerfIteration() { return iterNumFactor_; }

            void recordTape()
            {
                cl::Independent(spread_);
                calculateTotalSwaptionNpv();
                f_ = std::make_unique<cl::tape_function<double>>(spread_, swaptionNpv_);
            }

            // Calculates total calibration error.
            void calculateTotalSwaptionNpv()
            {
                swaptionNpv_.resize(depVarNumber());
                for (int i = 0; i < size_; i++)
                {
                    std::vector<cl::tape_double> temp = data_->swaptionNpvWithSpreadCorrection(spread_[i]);
                    for (int j = 0; j < depVarNumber(); j++)
                        swaptionNpv_[j] += temp[j];
                }
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                double h = 1.0e-10;  // shift for finite diff. method
                analyticalResults_.resize(size_ * depVarNumber());
                data_->calculateFinDiff(spread_
                                        , depVarNumber()
                                        , h
                                        , analyticalResults_
                                        , [this] (Real v) -> std::vector<Real>
                {
                    return data_->swaptionNpvWithSpreadCorrection(v);
                });
            }

            // Calculates derivatives using adjoint.
            void calcAdjoint()
            {
                adjointResults_ = f_->Jacobian(doubleSpread_);
            }


            double relativeTol() const { return 1e-2; }

            double absTol() const { return 1e-10; }

            Size size_;
            Size iterNumFactor_;
            SpreadCorrectionTestData* data_;
            std::vector<double> doubleSpread_;
            std::vector<cl::tape_double> spread_;
            std::vector<cl::tape_double> swaptionNpv_;
        };

        SpreadCorrectionTestData()
            : SwaptionData()

            , outPerform_(OUTPUT_FOLDER_NAME "//SpreadCorrection"
            , { { "filename", "AdjointPerformance" }
        , { "not_clear", "Not" }
        , { "line_box_width", "-5" }
        , { "title", "Swaption NPV differentiation performance with respect to spread" }
        , { "ylabel", "Time (s)" }
        , { "xlabel", "Number of spreads" }
        , { "smooth", "default" }
        , { "cleanlog", "true" } })

            , outAdjoint_(OUTPUT_FOLDER_NAME "//SpreadCorrection"
            , { { "filename", "Adjoint" }
        , { "not_clear", "Not" }
        , { "smooth", "default" }
        , { "title", "Swaption NPV adjoint differentiation performance with respect to spread" }
        , { "cleanlog", "false" }
        , { "ylabel", "Time (s)" }
        , { "xlabel", "Number of spreads" } })

            , outSize_(OUTPUT_FOLDER_NAME "//SpreadCorrection"
            , { { "filename", "TapeSize" }
        , { "not_clear", "Not" }
        , { "title", "Tape size dependence on number of spreads" }
        , { "cleanlog", "false" }
        , { "smooth", "default" }
        , { "ylabel", "Size (MB)" }
        , { "xlabel", "Number of spreads" } })

            , out_(OUTPUT_FOLDER_NAME "//SpreadCorrection//output"
            , { { "filename", "SwaptNPVonSpreads" }
        , { "ylabel", "Swaption NPV" }
        , { "not_clear", "Not" }
        , { "title", "Swaption NPV dependence on spread" }
        , { "cleanlog", "false" }
        , { "smooth", "default" }
        , { "xlabel", "Spread" } })

#if defined CL_GRAPH_GEN
            , pointNo_(20)
            , iterNo_(20)
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

        // Makes plots for strike sensitivity dependence.
        bool recordDependencePlot()
        {
            std::vector<Variation> outData(pointNo_);
            auto test = getTest(pointNo_);
            for (Size i = 0; i < pointNo_; i++)
            {
                outData[i].param_ = test->spread_[i];
                outData[i].swaptionNpv_ = test->data_->swaptionNpvWithSpreadCorrection(test->spread_[i])[0];
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
    typedef SpreadCorrectionTestData::Test SpreadCorrectionTest;

    struct CachedValueTestData
        : public SwaptionData
    {
        struct Test
        : public cl::AdjointTest<Test>
        {

            Test(Size size, CachedValueTestData* data)
            : size_(size)
            , data_(data)
            , volatility_(size)
            , swaptionNpv_(1)
            , iterNumFactor_(1)
            {
                setLogger(&data_->outPerform_);

                Real startVol = 0.15;
                Real maxVol = 0.25;
                Real step = (maxVol - startVol) / size_;

                for (Size i = 0; i < size_; i++)
                {
                    volatility_[i] = startVol + i*step;
                }


            }

            Size indepVarNumber() { return size_; }

            Size depVarNumber() { return 1; }

            Size minPerfIteration() { return iterNumFactor_; }

            void recordTape()
            {
                cl::Independent(volatility_);
                calculateTotalSwaptionNpv();
                f_ = std::make_unique<cl::tape_function<double>>(volatility_, swaptionNpv_);
            }

            // Calculates total calibration error.
            void calculateTotalSwaptionNpv()
            {
                for (int i = 0; i < size_; i++)
                    swaptionNpv_[0] += data_->swaptionNpvOnCachedValues(volatility_[i])[0];
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                double h = 1.0e-10;  // shift for finite diff. method
                analyticalResults_.resize(size_ * depVarNumber());
                data_->calculateFinDiff(volatility_
                                        , depVarNumber()
                                        , h
                                        , analyticalResults_
                                        , [this] (Real v) -> std::vector<Real>
                {
                    return data_->swaptionNpvOnCachedValues(v);
                });
            }

            double relativeTol() const { return 1e-2; }

            double absTol() const { return 1e-10; }

            Size size_;
            Size iterNumFactor_;
            CachedValueTestData* data_;
            std::vector<cl::tape_double> volatility_;
            std::vector<cl::tape_double> swaptionNpv_;
        };

        CachedValueTestData()
            : SwaptionData()

            , outPerform_(OUTPUT_FOLDER_NAME "//CachedValue"
            , { { "filename", "AdjointPerformance" }
        , { "not_clear", "Not" }
        , { "line_box_width", "-5" }
        , { "title", "Swaption NPV differentiation performance with respect to volatility" }
        , { "ylabel", "Time (s)" }
        , { "xlabel", "Number of volatilities" }
        , { "smooth", "default" }
        , { "cleanlog", "true" } })

            , outAdjoint_(OUTPUT_FOLDER_NAME "//CachedValue"
            , { { "filename", "Adjoint" }
        , { "not_clear", "Not" }
        , { "smooth", "default" }
        , { "title", "Swaption NPV adjoint differentiation performance with respect to volatility" }
        , { "cleanlog", "false" }
        , { "ylabel", "Time (s)" }
        , { "xlabel", "Number of volatilities" } })

            , outSize_(OUTPUT_FOLDER_NAME "//CachedValue"
            , { { "filename", "TapeSize" }
        , { "not_clear", "Not" }
        , { "title", "Tape size dependence on number of volatilities" }
        , { "cleanlog", "false" }
        , { "smooth", "default" }
        , { "ylabel", "Size (MB)" }
        , { "xlabel", "Number of volatilities" } })

            , out_(OUTPUT_FOLDER_NAME "//CachedValue//output"
            , { { "filename", "SwaptNPVonVol" }
        , { "ylabel", "Swaption NPV" }
        , { "not_clear", "Not" }
        , { "title", "Swaption NPV dependence on volatilities" }
        , { "cleanlog", "false" }
        , { "smooth", "default" }
        , { "xlabel", "Volatility" } })

#if defined CL_GRAPH_GEN
            , pointNo_(400)
            , iterNo_(16)
            , step_(25)
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

        // Makes plots for strike sensitivity dependence.
        bool recordDependencePlot()
        {
            std::vector<Variation> outData(pointNo_);
            auto test = getTest(pointNo_);
            for (Size i = 0; i < pointNo_; i++)
            {
                outData[i].param_ = test->volatility_[i];
                outData[i].swaptionNpv_ = test->data_->swaptionNpvOnCachedValues(test->volatility_[i])[0];
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
    typedef CachedValueTestData::Test CachedValueTest;
}

#endif