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


#ifndef cl_adjoint_test_zero_spreaded_term_structure_impl_hpp
#define cl_adjoint_test_zero_spreaded_term_structure_impl_hpp
#pragma once

#include "adjointzerospreadedtermstructuretest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include <ql/termstructures/yield/piecewisezerospreadedtermstructure.hpp>
#include <ql/termstructures/yield/zerocurve.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/termstructures/yield/ratehelpers.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/thirty360.hpp>
#include <ql/math/interpolations/all.hpp>
#include <boost/make_shared.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;

#define OUTPUT_FOLDER_NAME "AdjointPiecwZeroSprTermStruct"

namespace
{
    enum
    {
#if defined CL_GRAPH_GEN
        // Number of points for dependency plots.
        pointNo = 100,
        // Number of points for performance plot.
        iterNo = 100,
        // Step for portfolio size for performance testing .
        step = 1,
        // Defines performance accuracy. Its value is a minimum number
        // of calling of O(1) complexity methods per one performance test.
        iterNumFactor = 1000,
#else
        // Number of points for dependency plots.
        pointNo = 10,
        // Number of points for performance plot.
        iterNo = 10,
        // Step for portfolio size for performance testing .
        step = 1,
        // Defines performance accuracy. Its value is a minimum number
        // of calling of O(1) complexity methods per one performance test.
        iterNumFactor = 1
#endif
    };

    std::vector<double> getDoubleVector(Size size)
    {
        std::vector<double> vec =
        {
            0.04, 0.09, 0.12, 0.06, 0.05, 0.19, 0.52, 0.66, 0.45, 0.82, 0.73, 0.33, 0.82
            , 0.66, 0.39, 0.74, 0.58, 0.55, 0.36, 0.91, 0.11, 0.24, 0.5, 0.22, 0.43, 0.5
            , 0.99, 0.64, 0.58, 0.35, 0.46, 0.86, 0.71, 0.23, 0.3, 0.07, 0.03, 0.82, 0.12
            , 0.55, 0.82, 0.95, 0.11, 0.91, 0.16, 0.73, 0.24, 0.96, 0.31, 0.32, 0.46, 0.88
            , 0.71, 0.03, 0.55, 0.58, 0.1, 0.63, 0.23, 0.02, 0.62, 0.7, 0.37, 0.21, 0.19
            , 0.13, 0.2, 0.46, 0.41, 0.9, 0.99, 0.01, 0.73, 0.19, 0.2, 0.58, 0.59, 0.83
            , 0.1, 0.85, 0.76, 0.9, 0.49, 0.01, 0.65, 0.99, 0.5, 0.63, 0.52, 0.43, 0.79
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

    struct Interpol
    {
        Interpol() :
        calendar_(TARGET())
        , settlementDays_(2)
        , today_(Date(9, June, 2009))
        , compounding_(Continuous)
        , dayCount_(Actual360())
        , settlementDate_(calendar_.advance(today_, settlementDays_, Days))
        {
            Settings::instance().evaluationDate() = today_;
        }

        virtual Date setInterpolationDate() = 0;

        virtual void setSpreads(std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates) = 0;

        virtual boost::shared_ptr<ZeroYieldStructure> createSpreadedTermStructure(boost::shared_ptr<YieldTermStructure>& termStructure
                                                                                  , std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates)
                                                                                  = 0;
        virtual std::string getName() = 0;

        Calendar calendar_;
        Natural settlementDays_;
        DayCounter dayCount_;
        Compounding compounding_;
        Date today_;
        Date settlementDate_;

        SavedSettings backup_;

    };

    struct FlatInterpolLeft : Interpol
    {
        FlatInterpolLeft()
        {
        }

        Date setInterpolationDate() { return calendar_.advance(today_, 6, Months); }

        void setSpreads(std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates)
        {
            spreads.resize(2);
            spreadDates.resize(2);
            for (Size i = 0; i < 2; i++)
            {
                spreads[i] = Handle<Quote>(boost::make_shared<SimpleQuote>(0.01*i + 0.02));
                spreadDates[i] = calendar_.advance(today_, 7 * i + 8, Months);
            }
        }

        boost::shared_ptr<ZeroYieldStructure> createSpreadedTermStructure(boost::shared_ptr<YieldTermStructure>& termStructure
                                                                          , std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates)
        {
            return boost::make_shared<PiecewiseZeroSpreadedTermStructure>(
                Handle<YieldTermStructure>(termStructure),
                spreads, spreadDates);
        }

        std::string getName() { return "FlatInterpolationLeft"; }

    };

    struct FlatInterpolRight : Interpol
    {
        FlatInterpolRight() {}

        Date setInterpolationDate() { return calendar_.advance(today_, 6, Months); }

        void setSpreads(std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates)
        {
            spreads.resize(2);
            spreadDates.resize(2);
            for (Size i = 0; i < 2; i++)
            {
                spreads[i] = Handle<Quote>(boost::make_shared<SimpleQuote>(0.01*i + 0.02));
                spreadDates[i] = calendar_.advance(today_, 7 * i + 8, Months);
            }
        }

        boost::shared_ptr<ZeroYieldStructure> createSpreadedTermStructure(boost::shared_ptr<YieldTermStructure>& termStructure
                                                                          , std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates)
        {
            boost::shared_ptr<ZeroYieldStructure> sprTermStructure =
                boost::make_shared<PiecewiseZeroSpreadedTermStructure>(
                Handle<YieldTermStructure>(termStructure),
                spreads, spreadDates);
            sprTermStructure->enableExtrapolation();
            return sprTermStructure;
        }

        std::string getName() { return "FlatInterpolationRight"; }
    };

    struct LinearInterpolMultipleSpreads : Interpol
    {
        LinearInterpolMultipleSpreads() {}

        Date setInterpolationDate() { return  calendar_.advance(today_, 120, Days); }

        void setSpreads(std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates)
        {
            spreads.resize(4);
            spreadDates.resize(4);
            for (Size i = 0; i < 4; i++)
            {
                spreads[i] = Handle<Quote>(boost::make_shared<SimpleQuote>(0.005*i + 0.02));
            }
            for (Size i = 0; i < 2; i++)
            {
                spreadDates[i] = calendar_.advance(today_, 3 * i + 3, Months);
                spreadDates[i + 2] = calendar_.advance(today_, 10 * i + 30, Months);
            }
        }

        boost::shared_ptr<ZeroYieldStructure> createSpreadedTermStructure(boost::shared_ptr<YieldTermStructure>& termStructure
                                                                          , std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates)
        {
            return boost::make_shared<PiecewiseZeroSpreadedTermStructure>(
                Handle<YieldTermStructure>(termStructure),
                spreads, spreadDates);
        }

        std::string getName() { return "LinearInterpolationMultipleSpreads"; }
    };

    struct LinearInterpol : Interpol
    {
        LinearInterpol() {}

        Date setInterpolationDate() { return  calendar_.advance(today_, 120, Days); }

        void setSpreads(std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates)
        {
            spreads.resize(2);
            spreadDates.resize(2);
            for (Size i = 0; i < 2; i++)
            {
                spreads[i] = Handle<Quote>(boost::make_shared<SimpleQuote>(0.01*i + 0.02));
                spreadDates[i] = calendar_.advance(today_, 50 * i + 100, Days);
            }
        }

        boost::shared_ptr<ZeroYieldStructure> createSpreadedTermStructure(boost::shared_ptr<YieldTermStructure>& termStructure
                                                                          , std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates)
        {
            return boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<Linear> >(
                Handle<YieldTermStructure>(termStructure),
                spreads, spreadDates);
        }

        std::string getName() { return "LinearInterpolation"; }
    };

    struct ForwardFlatInterpol : Interpol
    {
        ForwardFlatInterpol() {}

        Date setInterpolationDate() { return  calendar_.advance(today_, 100, Days); }

        void setSpreads(std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates)
        {
            spreads.resize(2);
            spreadDates.resize(2);
            for (Size i = 0; i < 2; i++)
            {
                spreads[i] = Handle<Quote>(boost::make_shared<SimpleQuote>(0.01*i + 0.02));
                spreadDates[i] = calendar_.advance(today_, 185 * i + 75, Days);
            }
        }

        boost::shared_ptr<ZeroYieldStructure> createSpreadedTermStructure(boost::shared_ptr<YieldTermStructure>& termStructure
                                                                          , std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates)
        {
            return boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<ForwardFlat> >(
                Handle<YieldTermStructure>(termStructure),
                spreads, spreadDates);
        }

        std::string getName() { return "ForwardFlatInterpolation"; }
    };

    struct BackwardFlatInterpol : Interpol
    {
        BackwardFlatInterpol() {}

        Date setInterpolationDate() { return  calendar_.advance(today_, 110, Days); }

        void setSpreads(std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates)
        {
            spreads.resize(3);
            spreadDates.resize(3);
            for (Size i = 0; i < 3; i++)
            {
                spreads[i] = Handle<Quote>(boost::make_shared<SimpleQuote>(0.01*i + 0.02));
                spreadDates[i] = calendar_.advance(today_, 100 * (i + 1), Days);
            }
        }

        boost::shared_ptr<ZeroYieldStructure> createSpreadedTermStructure(boost::shared_ptr<YieldTermStructure>& termStructure
                                                                          , std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates)
        {
            return boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<BackwardFlat> >(
                Handle<YieldTermStructure>(termStructure),
                spreads, spreadDates);
        }

        std::string getName() { return "BackwardFlatInterpolation"; }
    };

    struct DefaultInterpol : Interpol
    {
        DefaultInterpol() {}

        Date setInterpolationDate() { return  calendar_.advance(today_, 100, Days); }

        void setSpreads(std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates)
        {
            spreads.resize(2, Handle<Quote>(boost::make_shared<SimpleQuote>(0.02)));
            spreadDates.resize(2);
            for (Size i = 0; i < 2; i++)
            {
                spreadDates[i] = calendar_.advance(today_, 75 * (i + 1), Days);
            }
        }

        boost::shared_ptr<ZeroYieldStructure> createSpreadedTermStructure(boost::shared_ptr<YieldTermStructure>& termStructure
                                                                          , std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates)
        {
            return boost::make_shared<PiecewiseZeroSpreadedTermStructure>(
                Handle<YieldTermStructure>(termStructure),
                spreads, spreadDates);
        }

        std::string getName() { return "DefaultInterpolation"; }
    };

    struct SetInterpolFactory : Interpol
    {
        SetInterpolFactory() {}

        Date setInterpolationDate() { return  calendar_.advance(today_, 11, Months); }

        void setSpreads(std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates)
        {
            spreads.resize(3);
            spreadDates.resize(3);
            for (Size i = 0; i < 3; i++)
            {
                spreads[i] = Handle<Quote>(boost::make_shared<SimpleQuote>(0.02 + i*pow(-1, i)*0.01));
                spreadDates[i] = calendar_.advance(today_, 8 * (i + 1), Months);
            }
        }

        boost::shared_ptr<ZeroYieldStructure> createSpreadedTermStructure(boost::shared_ptr<YieldTermStructure>& termStructure
                                                                          , std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates)
        {
            Cubic factory = Cubic(CubicInterpolation::Spline, false);
            return boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<Cubic> >(
                Handle<YieldTermStructure>(termStructure),
                spreads, spreadDates, compounding_,
                NoFrequency, dayCount_, factory);
        }

        std::string getName() { return "SetInterpolationFactory"; }
    };

    struct QuoteChanging : Interpol
    {
        QuoteChanging() {}

        Date setInterpolationDate() { return  calendar_.advance(today_, 120, Days); }

        void setSpreads(std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates)
        {
            spreads.resize(2);
            spreadDates.resize(2);
            for (Size i = 0; i < 2; i++)
            {
                spreads[i] = Handle<Quote>(boost::make_shared<SimpleQuote>(0.01*i + 0.02));
                spreadDates[i] = calendar_.advance(today_, 50 * i + 100, Days);
            }
        }

        boost::shared_ptr<ZeroYieldStructure> createSpreadedTermStructure(boost::shared_ptr<YieldTermStructure>& termStructure
                                                                          , std::vector<Handle<Quote> >& spreads, std::vector<Date>& spreadDates)
        {
            return boost::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<BackwardFlat> >(
                Handle<YieldTermStructure>(termStructure),
                spreads, spreadDates);
        }

        std::string getName() { return "QuoteChanging"; }
    };

    template <class InterpolationType>
    struct TestData
    {
        struct Test
        : cl::AdjointTest<Test>
        {
            explicit Test(Size size, TestData* data)
            : AdjointTest()
            , data_(data)
            , size_(size)
            , rates_(size)
            , resRates_(1)
            , spreads_()
            , spreadDates_()
            , interpolationDate_()
            {
                setLogger(&data_->outPerform_);
                std::vector<double> ratesDouble_ = getDoubleVector(size_);

                for (Size i = 0; i < size_; i++)
                {
                    rates_[i] = Real(ratesDouble_[i]);
                }
            }

            Size indepVarNumber() { return size_; }

            Size depVarNumber() { return 1; }

            Size minPerfIteration() { return iterNumFactor; }

            void recordTape()
            {
                cl::Independent(rates_);
                calculateResRates(rates_, resRates_);
                f_ = std::make_unique<cl::tape_function<double>>(rates_, resRates_);
            }

            void calculateResRates(std::vector<Real> rates, std::vector<Real>& resRates)
            {
                resRates.resize(1);
                data_->type_.setSpreads(spreads_, spreadDates_);
                interpolationDate_ = data_->type_.setInterpolationDate();
                boost::shared_ptr<YieldTermStructure> termStructure = createTermStructure(rates);
                boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
                    data_->type_.createSpreadedTermStructure(termStructure, spreads_, spreadDates_);

                Time t = data_->type_.dayCount_.yearFraction(data_->type_.today_, interpolationDate_);

                resRates[0] = spreadedTermStructure->zeroRate(t, data_->type_.compounding_);
            }

            void resetSpreads(std::vector<Real> rates, std::vector<Real>& resRates)
            {
                resRates.resize(1);
                data_->type_.setSpreads(spreads_, spreadDates_);
                interpolationDate_ = data_->type_.setInterpolationDate();
                spreads_[1] = Handle<Quote>(boost::make_shared<SimpleQuote>(0.025));
                boost::shared_ptr<YieldTermStructure> termStructure = createTermStructure(rates);
                boost::shared_ptr<ZeroYieldStructure> spreadedTermStructure =
                    data_->type_.createSpreadedTermStructure(termStructure, spreads_, spreadDates_);

                Time t = data_->type_.dayCount_.yearFraction(data_->type_.today_, interpolationDate_);

                resRates[0] = spreadedTermStructure->zeroRate(t, data_->type_.compounding_);
            }


            boost::shared_ptr<YieldTermStructure> createTermStructure(std::vector<Real>& rates)
            {
                std::vector<Integer> ts = getIntegerVector(size_ - 1);
                std::vector<Date> dates(size_);
                dates[0] = data_->type_.settlementDate_;
                for (Size i = 0; i < size_ - 1; ++i)
                {
                    dates[i + 1] = data_->type_.calendar_.advance(data_->type_.today_, ts[i], Days);
                }
                return boost::make_shared<ZeroCurve>(dates, rates, data_->type_.dayCount_);
            };

            void calcAnalytical()
            {
                double h = 1e-4;
                std::vector<Real> resRight(1), resLeft(1);
                analyticalResults_.resize(size_);
                for (Size i = 0; i < size_; i++)
                {
                    rates_[i] += h;
                    calculateResRates(rates_, resRight);
                    rates_[i] -= 2 * h;
                    calculateResRates(rates_, resLeft);
                    rates_[i] += h;
                    analyticalResults_[i] = (resRight[0] - resLeft[0]) / (2 * h);
                }
            }

            Size size_;
            std::vector<cl::tape_double> rates_;
            std::vector<cl::tape_double> resRates_;
            std::vector<Handle<Quote> > spreads_;
            std::vector<Date> spreadDates_;
            Date interpolationDate_;
            TestData* data_;

        };

        TestData() :
            type_()
            , outPerform_tilte_("Interpolated zero rate (by " + type_.getName() + ") \\ndifferentiation performance with respect to yield rate")
            , outAdjoint_tilte_("Interpolated zero rate (by " + type_.getName() + ") \\nadjoint differentiation performance with respect to yield rate")
            , out_tilte_("Zero Rate dependence on yield rate by" + type_.getName())
            , out_filename_("ZeroRateOnYieldRateBy" + type_.getName())
            , outPerform_(OUTPUT_FOLDER_NAME "//" + type_.getName()
            , { { "filename", "AdjointPerformance" }
        , { "not_clear", "Not" }
        , { "title", outPerform_tilte_ }
        , { "ylabel", "Time (s)" }
        , { "xlabel", "Number of yield rates" }
        , { "smooth", "default" }
        , { "line_box_width", "-5" }
        , { "cleanlog", "true" } })
            , outAdjoint_(OUTPUT_FOLDER_NAME "//" + type_.getName()
            , { { "filename", "Adjoint" }
        , { "not_clear", "Not" }
        , { "title", outAdjoint_tilte_ }
        , { "xlabel", "Number of yield rates" }
        , { "smooth", "default" }
        , { "ylabel", "Time (s)" }
        , { "cleanlog", "false" } })
            , outSize_(OUTPUT_FOLDER_NAME "//" + type_.getName()
            , { { "filename", "TapeSize" }
        , { "not_clear", "Not" }
        , { "title", "Tape size dependence on number of yield rates" }
        , { "xlabel", "Number of yield rates" }
        , { "ylabel", "Memory(MB)" }
        , { "cleanlog", "false" } })
            , out_(OUTPUT_FOLDER_NAME "//" + type_.getName() + "//output"
            , { { "filename", out_filename_ }
        , { "not_clear", "Not" }
        , { "title", out_tilte_ }
        , { "xlabel", "Yield rates" }
        , { "ylabel", "Zero Rate" }
        , { "cleanlog", "false" } })
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
            return std::make_shared<Test>(10 + size, this);
        }

        // Makes plots for sensitivity dependence.
        bool recordDependencePlot()
        {
            std::vector<RateDependence> outData(pointNo);
            auto test = getTest(10);
            std::vector<Real> rates = test->rates_;
            std::vector<Real> result(1);
            for (Size i = 0; i < pointNo; i++)
            {
                rates[5] = (i + 1)*0.01;
                test->calculateResRates(rates, result);
                outData[i] = { rates[5], result[0] };
            }
            out_ << outData;
            return true;
        }

        InterpolationType type_;
        std::string outPerform_tilte_;
        std::string outAdjoint_tilte_;
        std::string out_tilte_;
        std::string out_filename_;
        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output out_;

    };
}
#endif