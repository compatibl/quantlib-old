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

#ifndef cl_adjoint_term_structure_impl_hpp
#define cl_adjoint_term_structure_impl_hpp
#pragma once

#include <ql/quantlib.hpp>
#include "utilities.hpp"
#include "adjointtermstructuretest.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include <boost/make_shared.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;

#define OUTPUT_FOLDER_NAME "AdjointTermStructure"

namespace
{
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

    struct RateConsistency
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "Swap Rate", ""
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, RateConsistency& v)
        {
                stm << v.inputRate_
                    << ";" << v.functionValue_ << std::endl;

                return stm;
            }

        Real inputRate_;
        Real functionValue_;
    };

    struct Datum
    {
        Integer n;
        TimeUnit units;
        Rate rate;
    };

    struct TermStructureData
    {
        TermStructureData()
        : calendar_(TARGET())
        , settlementDays_(2)
        {
        }

        std::vector<Datum> createSwapData(std::vector<Real> rateValues)
        {
            Size rateSize = rateValues.size();
            std::vector<Datum> swapData(rateSize);
            for (Integer i = 0; i < rateSize; i++)
            {
                swapData[i] = { i * 3 + 10, Months, rateValues[i] };
            }
            return swapData;
        }

        void setTermStructure(std::vector<Datum> swapData)
        {
            std::vector<Datum> depositData =
            {
                { 1, Months, 4.581 },
                { 2, Months, 4.573 },
                { 3, Months, 4.557 },
                { 6, Months, 4.496 },
                { 9, Months, 4.490 }
            };

            Date today = calendar_.adjust(Date::todaysDate());
            Settings::instance().evaluationDate() = today;
            Date settlement = calendar_.advance(today, settlementDays_, Days);
            Size deposits = depositData.size(),
                swaps = swapData.size();

            std::vector<boost::shared_ptr<RateHelper> > instruments(deposits + swaps);
            for (Size i = 0; i < deposits; i++)
            {
                instruments[i] = boost::make_shared<DepositRateHelper>(
                    depositData[i].rate / 100
                    , depositData[i].n*depositData[i].units
                    , settlementDays_
                    , calendar_
                    , ModifiedFollowing
                    , true
                    , Actual360());
            }
            boost::shared_ptr<IborIndex> index(new IborIndex("dummy"
                , 6 * Months
                , settlementDays_
                , Currency()
                , calendar_
                , ModifiedFollowing
                , false
                , Actual360()));
            for (Size i = 0; i < swaps; ++i)
            {
                instruments[i + deposits] = boost::make_shared<SwapRateHelper>(
                    swapData[i].rate / 100
                    , swapData[i].n*swapData[i].units
                    , calendar_
                    , Annual
                    , Unadjusted
                    , Thirty360()
                    , index);
            }
            termStructure_ = boost::make_shared<PiecewiseYieldCurve<Discount, LogLinear>>(settlement,
                                                                                          instruments, Actual360());

        }

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
                , 0.65, 0.74, 0.19, 0.75, 0.06, 0.08, 0.51, 0.97, 0.10, 0.99, 0.32, 0.40
                , 0.69, 0.71, 0.83, 0.64, 0.08, 0.61, 0.56, 0.27, 0.15, 0.69, 0.70, 0.89
                , 0.52, 0.15
            };
            return std::vector<double>(vec.begin(), vec.begin() + size);
        }

        std::vector<Real> calculateZeroValue(std::vector<Datum> swapData)
        {
            setTermStructure(swapData);

            Handle<Quote> spread(boost::make_shared<SimpleQuote>(0.01));

            boost::shared_ptr<YieldTermStructure> spreaded(
                new ZeroSpreadedTermStructure(Handle<YieldTermStructure>(termStructure_), spread));

            Date testDate = termStructure_->referenceDate() + 5 * Years;
            DayCounter resultDayCounter = termStructure_->dayCounter();

            Rate zero = termStructure_->zeroRate(testDate, resultDayCounter, Continuous, NoFrequency);

            return std::vector<Real>(1, zero);
        }

        std::vector<Real> calculateForwardValue(std::vector<Datum> swapData)
        {
            setTermStructure(swapData);

            Handle<Quote> spread(boost::shared_ptr<Quote>(new SimpleQuote(0.01)));

            boost::shared_ptr<YieldTermStructure> spreaded(
                new ForwardSpreadedTermStructure(Handle<YieldTermStructure>(termStructure_), spread));

            Date testDate = termStructure_->referenceDate() + 5 * Years;
            DayCounter resultDayCounter = termStructure_->dayCounter();

            Rate forward = termStructure_->forwardRate(testDate, testDate, resultDayCounter, Continuous, NoFrequency);

            return std::vector<Real>(1, forward);
        }

        std::vector<Real> calculateDiscountValues(std::vector<Datum> swapData)
        {
            setTermStructure(swapData);

            Date today = Settings::instance().evaluationDate();
            Date newToday = today + 3 * Years;
            Date newSettlement = calendar_.advance(newToday, settlementDays_, Days);
            Date testDate = newSettlement + 5 * Years;

            boost::shared_ptr<YieldTermStructure> implied(
                new ImpliedTermStructure(Handle<YieldTermStructure>(termStructure_), newSettlement));

            DiscountFactor baseDiscount = termStructure_->discount(newSettlement);
            DiscountFactor discount = termStructure_->discount(testDate);
            DiscountFactor impliedDiscount = implied->discount(testDate);

            return std::vector<Real>{ baseDiscount, discount, impliedDiscount };
        }

        template <class Func>
        void calculateFinDiff(std::vector<Real>& rateValues
                              , Real h
                              , std::vector<Real>& sf_Finite
                              , Func getValue)
        {
            Size sizeX = rateValues.size();

            std::vector<Real> value = getValue(createSwapData(rateValues));

            Size sizeY = value.size();

            sf_Finite.resize(sizeX * sizeY);

            for (Size i = 0; i < sizeX; i++)
            {
                rateValues[i] += h;
                std::vector<Real> finiteValue = getValue(createSwapData(rateValues));
                rateValues[i] -= h;

                for (Size j = 0; j < sizeY; j++)
                {
                    sf_Finite[sizeX*j + i] = (finiteValue[j] - value[j]) / h;
                }
            }
        }

        // Global data.
        Calendar calendar_;
        Natural settlementDays_;
        boost::shared_ptr<YieldTermStructure> termStructure_;

        // Cleanup.
        SavedSettings backup_;

    };

    struct ZSpreadedTestData
        : public TermStructureData
    {
        struct Test
        : public cl::AdjointTest<Test>
        {
            Test(Size size, ZSpreadedTestData* data)
            : size_(size)
            , data_(data)
            , rate_(size)
            , iterNumFactor_(10)
            {
                setLogger(&data_->outPerform_);

                std::vector<double> doubleRate = data_->getDoubleVector(size);
                for (Size i = 0; i < size; i++)
                    rate_[i] = doubleRate[i];
            }

            Size indepVarNumber() { return size_; }

            Size depVarNumber() { return 1; }

            Size minPerfIteration() { return iterNumFactor_; }

            void recordTape()
            {
                cl::Independent(rate_);
                calculateZeroValue();
                f_ = std::make_unique<cl::tape_function<double>>(rate_, zeroValue_);
            }

            // Calculates total calibration error.
            void calculateZeroValue()
            {
                zeroValue_ = data_->calculateZeroValue(data_->createSwapData(rate_));
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                double h = 1.0e-6;  // shift for finite diff. method
                analyticalResults_.resize(size_);
                data_->calculateFinDiff(rate_
                                        , h
                                        , analyticalResults_
                                        , [this] (std::vector<Datum>& v) -> std::vector<Real>
                {
                    return data_->calculateZeroValue(v);
                });
            }

            double relativeTol() const { return 1e-2; }

            double absTol() const { return 1e-5; }

            Size size_;
            Size iterNumFactor_;
            ZSpreadedTestData* data_;
            std::vector<cl::tape_double> rate_;
            std::vector<cl::tape_double> zeroValue_;
        };

        ZSpreadedTestData()
            : TermStructureData()

            , outPerform_(OUTPUT_FOLDER_NAME "//ConsistencyZeroSpreadedTermStructure"
            , { { "filename", "AdjointPerformance" }
        , { "not_clear", "Not" }
        , { "cleanlog", "true" }
        , { "title", "Zero rate differentiation performance with respect to swap rate" }
        , { "xlabel", "Number of swap rates" }
        , { "ylabel", "Time (s)" }
        , { "smooth", "default" }
        , { "line_box_width", "-5" } })

            , outAdjoint_(OUTPUT_FOLDER_NAME "//ConsistencyZeroSpreadedTermStructure"
            , { { "filename", "Adjoint" }
        , { "not_clear", "Not" }
        , { "cleanlog", "false" }
        , { "title", " Zero rate adjoint differentiation performance with respect to swap rate" }
        , { "xlabel", "Number of swap rates" }
        , { "smooth", "default" }
        , { "ylabel", "Time (s)" } })

            , outSize_(OUTPUT_FOLDER_NAME "//ConsistencyZeroSpreadedTermStructure"
            , { { "filename", "TapeSize" }
        , { "not_clear", "Not" }
        , { "cleanlog", "false" }
        , { "title", "Tape size dependence on number of swap rates" }
        , { "xlabel", "Number of swap rates" }
        , { "ylabel", "Memory(MB)" } })

            , out_(OUTPUT_FOLDER_NAME "//ConsistencyZeroSpreadedTermStructure//output"
            , { { "filename", "zeroRateOnSwapRate" }
        , { "not_clear", "Not" }
        , { "cleanlog", "false" }
        , { "title", "Zero Rate dependence on Swap Rate" }
        , { "ylabel", "Zero Rate" } })

#if defined CL_GRAPH_GEN
            , pointNo_(80)
            , iterNo_(80)
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
            return std::make_shared<Test>(size + 30, this);
        }

        // Makes plots for strike sensitivity dependence.
        bool recordDependencePlot()
        {
            std::vector<RateConsistency> outData(pointNo_);
            auto test = getTest(pointNo_);
            Rate rate = 0.0;
            for (Size i = 0; i < pointNo_; rate += 1.0, i++)
            {
                outData[i] = { rate / 100, test->data_->calculateZeroValue({ { 10, Years, rate } })[0] };
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
    typedef ZSpreadedTestData::Test ZSpreadedTest;

    struct FSpreadedTestData
        : public TermStructureData
    {
        struct Test
        : public cl::AdjointTest<Test>
        {
            Test(Size size, FSpreadedTestData* data)
            : size_(size)
            , data_(data)
            , rate_(size)
            , iterNumFactor_(10)
            {
                setLogger(&data_->outPerform_);

                std::vector<double> doubleRate = data_->getDoubleVector(size);
                for (Size i = 0; i < size; i++)
                    rate_[i] = doubleRate[i];
            }

            Size indepVarNumber() { return size_; }

            Size depVarNumber() { return 1; }

            Size minPerfIteration() { return iterNumFactor_; }

            void recordTape()
            {
                cl::Independent(rate_);
                calculateForwardValue();
                f_ = std::make_unique<cl::tape_function<double>>(rate_, forwardValue_);
            }

            // Calculates total calibration error.
            void calculateForwardValue()
            {
                forwardValue_ = data_->calculateForwardValue(data_->createSwapData(rate_));
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                double h = 1.0e-6;  // shift for finite diff. method
                analyticalResults_.resize(size_);
                data_->calculateFinDiff(rate_
                                        , h
                                        , analyticalResults_
                                        , [this] (std::vector<Datum>& v) -> std::vector<Real>
                {
                    return data_->calculateForwardValue(v);
                });
            }

            double relativeTol() const { return 1e-2; }

            double absTol() const { return 1e-5; }

            Size size_;
            Size iterNumFactor_;
            FSpreadedTestData* data_;
            std::vector<cl::tape_double> rate_;
            std::vector<cl::tape_double> forwardValue_;
        };

        FSpreadedTestData()
            : TermStructureData()

            , outPerform_(OUTPUT_FOLDER_NAME "//ConsistencyForwardSpreadedTermStructure"
            , { { "filename", "AdjointPerformance" }
        , { "not_clear", "Not" }
        , { "cleanlog", "true" }
        , { "title", "Forward rate differentiation performance with respect to swap rate" }
        , { "xlabel", "Number of swap rates" }
        , { "ylabel", "Time (s)" }
        , { "smooth", "default" }
        , { "line_box_width", "-5" } })

            , outAdjoint_(OUTPUT_FOLDER_NAME "//ConsistencyForwardSpreadedTermStructure"
            , { { "filename", "Adjoint" }
        , { "not_clear", "Not" }
        , { "cleanlog", "false" }
        , { "title", " Forward rate adjoint differentiation performance with respect to swap rate" }
        , { "xlabel", "Number of swap rates" }
        , { "smooth", "default" }
        , { "ylabel", "Time (s)" } })

            , outSize_(OUTPUT_FOLDER_NAME "//ConsistencyForwardSpreadedTermStructure"
            , { { "filename", "TapeSize" }
        , { "not_clear", "Not" }
        , { "cleanlog", "false" }
        , { "title", "Tape size dependence on number of swap rates" }
        , { "xlabel", "Number of swap rates" }
        , { "ylabel", "Memory(MB)" } })

            , out_(OUTPUT_FOLDER_NAME "//ConsistencyForwardSpreadedTermStructure//output"
            , { { "filename", "forwardRateOnSwapRate" }
        , { "not_clear", "Not" }
        , { "cleanlog", "false" }
        , { "title", "Forward Rate dependence on Swap Rate" }
        , { "ylabel", "Forward Rate" } })

#if defined CL_GRAPH_GEN
            , pointNo_(80)
            , iterNo_(80)
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
            return std::make_shared<Test>(size + 30, this);
        }

        // Makes plots for strike sensitivity dependence.
        bool recordDependencePlot()
        {
            std::vector<RateConsistency> outData(pointNo_);
            auto test = getTest(pointNo_);
            Rate rate = 0.0;
            for (Size i = 0; i < pointNo_; rate += 1.0, i++)
            {
                outData[i] = { rate / 100, test->data_->calculateForwardValue({ { 10, Years, rate } })[0] };
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
    typedef FSpreadedTestData::Test FSpreadedTest;

    struct ImpliedTestData
        : public TermStructureData
    {
        struct Test
        : public cl::AdjointTest<Test>
        {
            static const CalcMethod default_method = other;

            Test(Size size, ImpliedTestData* data)
                : size_(size)
                , data_(data)
                , rate_(size)
                , iterNumFactor_(10)
            {
                setLogger(&data_->outPerform_);

                doubleRate_ = data_->getDoubleVector(size);

                std::copy(doubleRate_.begin(), doubleRate_.end(), rate_.begin());
            }

            Size indepVarNumber() { return size_; }

            Size depVarNumber() { return 3; }

            Size minPerfIteration() { return iterNumFactor_; }

            void recordTape()
            {
                cl::Independent(rate_);
                calculateDiscountValues();
                f_ = std::make_unique<cl::tape_function<double>>(rate_, discountValues_);
            }

            // Calculates total calibration error.
            void calculateDiscountValues()
            {
                discountValues_ = data_->calculateDiscountValues(data_->createSwapData(rate_));
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                double h = 1.0e-6;  // shift for finite diff. method
                analyticalResults_.resize(size_ * depVarNumber());
                data_->calculateFinDiff(rate_
                                        , h
                                        , analyticalResults_
                                        , [this] (std::vector<Datum>& v) -> std::vector<Real>
                {
                    return data_->calculateDiscountValues(v);
                });
            }

            // Calculates derivatives using adjoint.
            void calcAdjoint()
            {
                adjointResults_ = f_->Jacobian(doubleRate_);
            }

            double relativeTol() const { return 1e-2; }

            double absTol() const { return 1e-5; }

            Size size_;
            Size iterNumFactor_;
            ImpliedTestData* data_;
            std::vector<double> doubleRate_;
            std::vector<cl::tape_double> rate_;
            std::vector<cl::tape_double> discountValues_;
        };

        ImpliedTestData()
            : TermStructureData()

            , outPerform_(OUTPUT_FOLDER_NAME "//ConsistencyYieldTermStructure"
            , { { "filename", "AdjointPerformance" }
        , { "not_clear", "Not" }
        , { "cleanlog", "true" }
        , { "title", "Base discount, discount, implied discount\\ndifferentiation performance with respect to swap rate" }
        , { "xlabel", "Number of swap rates" }
        , { "ylabel", "Time (s)" }
        , { "smooth", "default" }
        , { "line_box_width", "-5" } })

            , outAdjoint_(OUTPUT_FOLDER_NAME "//ConsistencyYieldTermStructure"
            , { { "filename", "Adjoint" }
        , { "not_clear", "Not" }
        , { "cleanlog", "false" }
        , { "title", "Base discount, discount, implied discount\\nadjoint differentiation performance with respect to swap rate" }
        , { "xlabel", "Number of swap rates" }
        , { "smooth", "default" }
        , { "ylabel", "Time (s)" } })

            , outSize_(OUTPUT_FOLDER_NAME "//ConsistencyYieldTermStructure"
            , { { "filename", "TapeSize" }
        , { "not_clear", "Not" }
        , { "cleanlog", "false" }
        , { "title", "Tape size dependence on number of swap rates" }
        , { "xlabel", "Number of swap rates" }
        , { "ylabel", "Memory(MB)" } })

            , out_(OUTPUT_FOLDER_NAME "//ConsistencyYieldTermStructure//output"
            , { { "filename", "discountOnSwapRate" }
        , { "not_clear", "Not" }
        , { "cleanlog", "false" }
        , { "title", "Discount values dependence on Swap Rate" }
        , { "ylabel", "Forward Rate" } })

#if defined CL_GRAPH_GEN
            , pointNo_(80)
            , iterNo_(80)
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
            return std::make_shared<Test>(size + 30, this);
        }

        // Makes plots for strike sensitivity dependence.
        bool recordDependencePlot()
        {
            std::vector<RateDiscount> outData(pointNo_);
            auto test = getTest(pointNo_);
            Rate rate = 0.0;
            for (Size i = 0; i < pointNo_; rate += 1.0, i++)
            {
                std::vector<Real> out = test->data_->calculateDiscountValues({ { 10, Years, rate } });
                outData[i] = { rate / 100
                    , out[0]
                    , out[1]
                    , out[2] };
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
    typedef ImpliedTestData::Test ImpliedTest;

}

#endif