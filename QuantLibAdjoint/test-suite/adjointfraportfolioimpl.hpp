/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2006 Allen Kuo
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

// Based on FRA (Forward Rate Agreement) project from Quantlib.

#ifndef cl_adjoint_fra_portfolio_impl_hpp
#define cl_adjoint_fra_portfolio_impl_hpp
#pragma once

#include "adjointfraportfoliotest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"
#include <ql/cashflows/coupon.hpp>
#include <ql/time/daycounters/thirty360.hpp>
#include <ql/indexes/ibor/euribor.hpp>
#include <ql/time/schedule.hpp>
#include <boost/make_shared.hpp>
#include <ql/quantlib.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;

#define OUTPUT_FOLDER_NAME "AdjointFRAPortfolio"

namespace
{
    enum
    {
#if defined CL_GRAPH_GEN
        // Number of points for dependency plots.
        pointNo = 50,
        // Number of points for performance plot.
        iterNo = 50,
        // Step for portfolio size for performance testing .
        step = 1,
#else
        // Number of points for dependency plots.
        pointNo = 1,
        // Number of points for performance plot.
        iterNo = 1,
        // Step for portfolio size for performance testing .
        step = 1,
#endif
        // Defines performance accuracy. Its value is a minimum number
        // of calling of O(1) complexity methods per one performance test.
        iterNumFactor = 0,
    };

    struct RateVariation
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "Forward Rate", ""
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type&
            operator << (stream_type& stm, RateVariation& v)
        {
                stm << v.rate_
                    << ";" << v.contractValue_
                    << std::endl;
                return stm;
            }

        Real rate_;
        Real contractValue_;
    };

    struct FRAPortfolioData
    {
        FRAPortfolioData() :
        euriborTermStructure_()
        , euribor_(boost::make_shared<Euribor3M>(euriborTermStructure_))
        , todaysDate_(27, July, 2015)
        , calendar_(euribor_->fixingCalendar())
        , fixingDays_(euribor_->fixingDays())
        , settlementDate_(calendar_.advance(todaysDate_, fixingDays_, Days))
        , fraDayCounter_(euribor_->dayCounter())
        , convention_(euribor_->businessDayConvention())
        , endOfMonth_(euribor_->endOfMonth())
        , termStructureDayCounter_(ActualActual(ActualActual::ISDA))
        , fraFwdType_(Position::Long)
        , fraNotional_(100.0)
        , quoteHandles_()
        , fra_()
        , discountingTermStructure_()
        , fraTermStructure_()
        {
            Settings::instance().evaluationDate() = todaysDate_;
        }

        // This function creates term structure for Forward Rate Agreements (FRA)
        // based on values of forward rates. Each forward rate refers to
        // 3 month term FRA quotes (index refers to months to start) are considered.
        // Input parameters:
        //  - FraRate - vector of forward rate for every FRA in portfolio;
        void createTermStructure(std::vector<Real>& FraRate)
        {
            Size n = FraRate.size();

            //Create relinkable handle to rate for Forward Rate Agreement.
            quoteHandles_.resize(n);
            for (Size i = 0; i < n; ++i)
                quoteHandles_[i] = RelinkableHandle<Quote>(boost::make_shared<SimpleQuote>(FraRate[i]));

            // Rate helpers
            // RateHelpers are built from the above quotes together with
            // other instrument dependent infos.  Quotes are passed in
            // relinkable handles which could be relinked to some other
            // data source later.
            fra_.resize(n);
            for (Size i = 0; i < n; ++i)
                fra_[i] = boost::make_shared<FraRateHelper>(
                quoteHandles_[i], (i + 1), (i + 4), fixingDays_, calendar_, convention_,
                endOfMonth_, fraDayCounter_);

            double tolerance = 1.0e-15;

            // Create Forward Rate Agreement curve.
            fraTermStructure_ = boost::make_shared<PiecewiseYieldCurve<Discount, LogLinear>>(
                settlementDate_, fra_,
                termStructureDayCounter_, tolerance);

            // Term structures used for pricing/discounting.
            discountingTermStructure_.linkTo(fraTermStructure_);
            euriborTermStructure_.linkTo(fraTermStructure_);

            fraValueDate_.resize(n);
            fraMaturityDate_.resize(n);
            for (Size i = 0; i < n; i++)
            {
                fraValueDate_[i] = calendar_.advance(
                    settlementDate_, i + 1, Months,
                    convention_);

                // Set maturity date (The date on which the notional loan matures).
                fraMaturityDate_[i] = calendar_.advance(
                    fraValueDate_[i],  FraTermMonths_, Months,
                    convention_);
            }
        }

        // Adjust term structure with updating fra rate at specified posittion.
        void adjustTermStructure(Real FraRate, Size position)
        {
            quoteHandles_[position] = RelinkableHandle<Quote>(boost::make_shared<SimpleQuote>(FraRate));
            fra_[position] = boost::make_shared<FraRateHelper>(
                quoteHandles_[position], (position + 1), (position + 4), fixingDays_, calendar_, convention_,
                endOfMonth_, fraDayCounter_);

            double tolerance = 1.0e-15;

            // Create Forward Rate Agreement curve.
            fraTermStructure_ = boost::make_shared<PiecewiseYieldCurve<Discount, LogLinear>>(
                settlementDate_, fra_,
                termStructureDayCounter_, tolerance);

            // Term structures used for pricing/discounting.
            discountingTermStructure_.linkTo(fraTermStructure_);
            euriborTermStructure_.linkTo(fraTermStructure_);
        }

        // Create FRA for defined rate.
        boost::shared_ptr<ForwardRateAgreement> createFRA(Real rate, Size index)
        {
            return boost::make_shared<ForwardRateAgreement>(
                fraValueDate_[index], fraMaturityDate_[index],
                fraFwdType_, rate,
                fraNotional_, euribor_,
                discountingTermStructure_);
        }

        // Use finite difference to approximate the derivatives
        // of contract value of each Forward Rate Agreement on forward rate.
        // dy/dx = (y(x) - y(x-h))/(h)
        // Input parameters:
        //  - FraRate - vector of forward rate for every FRA in portfolio;
        ////  - h - step size for finite difference method;
        //  - sf_Finite - calculated derivatives.
        void calculateFinDiff(std::vector<Real>& FraRate, double h, std::vector<Real>& sf_Finite)
        {
            Size n = FraRate.size();
            sf_Finite.resize(n);
            createTermStructure(FraRate);
            for (Size i = 0; i < n; i++)
            {
                sf_Finite[i] += createFRA(FraRate[i], i)->forwardValue() / h;
            }
            for (Size i = 0; i < n; i++)
            {
                FraRate[i] -= h;
                adjustTermStructure(FraRate[i], i);
                sf_Finite[i] -= createFRA(FraRate[i], i)->forwardValue() / h;
                FraRate[i] += h;
            }
        }

        RelinkableHandle<YieldTermStructure> euriborTermStructure_;
        boost::shared_ptr<IborIndex> euribor_;
        Date todaysDate_;

        Calendar calendar_;
        Integer fixingDays_;
        Date settlementDate_;

        DayCounter fraDayCounter_;
        BusinessDayConvention convention_;
        bool endOfMonth_;
        DayCounter termStructureDayCounter_;

        Position::Type fraFwdType_;
        Real fraNotional_;
        const Integer  FraTermMonths_ = 3;
        RelinkableHandle<YieldTermStructure> discountingTermStructure_;
        std::vector<RelinkableHandle<Quote> > quoteHandles_;
        std::vector<boost::shared_ptr<RateHelper>> fra_;

        std::vector<Date> fraValueDate_;
        std::vector<Date> fraMaturityDate_;
        boost::shared_ptr<YieldTermStructure> fraTermStructure_;

    };


    struct TestData
        : public FRAPortfolioData
    {
        struct Test
        : public cl::AdjointTest<Test>
        {
            Test(Size size, TestData* data)
            : size_(size)
            , data_(data)
            , fraRate_(size)
            , totalContractValue_()
            {
                setLogger(&data_->outPerform_);

                for (Size i = 0; i < size_; i++)
                {
                    fraRate_[i] = 0.030 + i*0.0000001;
                }
            }

            Size indepVarNumber() { return size_; }

            Size depVarNumber() { return 1; }

            Size minPerfIteration() { return iterNumFactor; }

            void recordTape()
            {
                cl::Independent(fraRate_);
                calculateTotalContractValue();
                f_ = std::make_unique<cl::tape_function<double>>(fraRate_, totalContractValue_);
            }

            // Calculates price of portfolio and each option.
            void calculateTotalContractValue()
            {
                totalContractValue_.resize(1, 0);
                std::vector<boost::shared_ptr<ForwardRateAgreement>> fraPortfolio(size_, 0);
                data_->createTermStructure(fraRate_);

                for (Size i = 0; i < size_; i++)
                {
                    fraPortfolio[i] = data_->createFRA(fraRate_[i], i);
                    totalContractValue_.front() += fraPortfolio[i]->forwardValue();
                }
            }

            // Calculates derivatives using finite difference method.
            void calcAnalytical()
            {
                double h = 1.0e-6;  // shift for finite diff. method
                analyticalResults_.resize(size_);
                data_->calculateFinDiff(fraRate_, h, analyticalResults_);
            }

            double relativeTol() const { return 1e-5; }

            double absTol() const { return 1e-10; }

            Size size_;
            TestData* data_;
            std::vector<cl::tape_double> fraRate_;
            std::vector<cl::tape_double> totalContractValue_;
        };

        TestData()
            : FRAPortfolioData()
            , outPerform_(OUTPUT_FOLDER_NAME "//",
            { { "filename", "AdjointPerformance" }
        , { "not_clear", "Not" }
        , { "title", "FRA contract value differentiation performance with respect to forward rates" }
        , { "ylabel", "Time (s)" }
        , { "xlabel", "Number of forward rates" }
        , { "line_box_width", "-5" }
        , { "smooth", "default" }
        , { "cleanlog", "true" }
        })
            , outAdjoint_(OUTPUT_FOLDER_NAME "//",
            { { "filename", "Adjoint" }
        , { "not_clear", "Not" }
        , { "title", "FRA contract value adjoint differentiation performance with respect to forward rates" }
        , { "ylabel", "Time (s)" }
        , { "xlabel", "Number of forward rates" }
        , { "smooth", "default" }
        , { "cleanlog", "false" }
        })
            , outSize_(OUTPUT_FOLDER_NAME "//",
            { { "filename", "TapeSize" }
        , { "not_clear", "Not" }
        , { "title", "Tape size dependence on number of forward rates" }
        , { "ylabel", "Memory (MB)" }
        , { "xlabel", "Number of forward rates" }
        , { "smooth", "default" }
        , { "cleanlog", "false" }
        })
            , out_(OUTPUT_FOLDER_NAME "//output",
            { { "filename", "Output" }
        , { "not_clear", "Not" }
        , { "title", "FRA contract value dependence on FRA rate" }
        , { "ylabel", "FRA Contract Value" }
        , { "xlabel", "FRA Rate" }
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
            return std::make_shared<Test>(50 + size, this);
        }

        // Makes plots for strike sensitivity dependence.
        bool recordDependencePlot()
        {
            std::vector<RateVariation> outData(pointNo);
            auto test = getTest(pointNo);
            createTermStructure(test->fraRate_);
            for (Size i = 0; i < pointNo; i += 10)
            {
                boost::shared_ptr<ForwardRateAgreement> fra = createFRA(test->fraRate_[i], i);
                outData[i] = { test->fraRate_[i], fra->forwardValue() };
            }
            out_ << outData;
            return true;
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output out_;
    };


    typedef TestData::Test FRAPortfolioTest;
}

#endif